#!/usr/bin/env perl
# Stage the GV second resolvent in the directory containing completed RbXc vectors.

use strict;
use warnings;

use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Find qw(find);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir tempfile);
use FindBin;
use Getopt::Long qw(GetOptions);
use JSON::PP;
use POSIX qw(isfinite);
use Scalar::Util qw(looks_like_number);

my ($ocean_run_root, $work_dir, $legacy_target, $component, $cv_manifest);
my ($valence_con, $ocean_bin, $photon_energy, $con_bands);
my ($bse2_input, $energy_range, $core_eta);
my $edge_index = 1;
my $dry_run = 0;
my $help = 0;

GetOptions(
  'ocean-run-root=s' => \$ocean_run_root,
  'work-dir=s'       => \$work_dir,
  'target=s'         => \$legacy_target,
  'component=s'      => \$component,
  'edge-index=i'     => \$edge_index,
  'cv-manifest=s'    => \$cv_manifest,
  'valence-con=s'    => \$valence_con,
  'con-bands=s'      => \$con_bands,
  'photon-energy=f'  => \$photon_energy,
  'bse2-input=s'     => \$bse2_input,
  'energy-range=s'   => \$energy_range,
  'core-eta=f'       => \$core_eta,
  'ocean-bin=s'      => \$ocean_bin,
  'dry-run!'         => \$dry_run,
  'help|h'           => \$help,
) or usage(1);
usage(0) if $help;
die "--bse2-input and --energy-range are mutually exclusive\n"
  if defined($bse2_input) && defined($energy_range);
if (defined $bse2_input) {
  $bse2_input = abs_path($bse2_input) or die "Cannot resolve --bse2-input\n";
}
die "--core-eta must be finite and positive\n"
  if defined($core_eta) && (!isfinite($core_eta) || $core_eta <= 0);

for my $required (
  ['--ocean-run-root', $ocean_run_root], ['--component', $component],
  ['--valence-con', $valence_con],
) {
  die "$required->[0] is required\n" unless defined $required->[1];
}

$ocean_run_root = abs_path($ocean_run_root)
  or die "Cannot resolve OCEAN computation directory\n";
$valence_con = abs_path($valence_con)
  or die "Cannot resolve valence con.u2.dat\n";

if (defined $work_dir && defined $legacy_target) {
  my $a = File::Spec->rel2abs($work_dir);
  my $b = File::Spec->rel2abs($legacy_target);
  die "--work-dir and deprecated --target name different directories\n"
    unless $a eq $b;
}
$work_dir = $work_dir // $legacy_target // getcwd();
$work_dir = abs_path($work_dir)
  or die "The RbXc work directory must already exist\n";
die "The RbXc work path is not a directory: $work_dir\n" unless -d $work_dir;
die "Refusing to stage directly in the OCEAN computation root\n"
  if $work_dir eq $ocean_run_root;

$cv_manifest //= File::Spec->catfile($work_dir, 'shg_rbxc_cv.info');
$cv_manifest = abs_path($cv_manifest)
  or die "Cannot resolve CV manifest\n";
die "The CV manifest must be in the RbXc work directory: $work_dir\n"
  unless dirname($cv_manifest) eq $work_dir;

$component = lc($component);
$component =~ s/[\s,]+//g;
die "--component must contain exactly three x/y/z labels\n"
  unless $component =~ /\A[xyz]{3}\z/;
die "--edge-index must be a positive integer\n"
  unless $edge_index =~ /^\d+$/ && $edge_index >= 1;

for my $stage (qw(Common DFT OPF PREP SCREEN)) {
  my $path = File::Spec->catdir($ocean_run_root, $stage);
  die "Missing completed core input stage: $path\n" unless -d $path;
}
my $core_con = File::Spec->catfile($ocean_run_root, 'PREP', 'BSE', 'con.u2.dat');
die "Missing core conduction wavefunctions: $core_con\n" unless -f $core_con;
die "Empty valence conduction wavefunction file\n" unless -s $valence_con;
die "Core and valence con.u2.dat sizes differ\n"
  unless -s($core_con) == -s($valence_con);
my $conduction_bytes = -s $core_con;

my $cv = read_cv_manifest($cv_manifest);
if (defined $energy_range) {
  expand_target_loop($energy_range, scalar(@{$cv->{entries}}));
}
if (defined $bse2_input) {
  validate_bse2(read_text($bse2_input), scalar(@{$cv->{entries}}));
}
die "CV manifest middle/right axes do not match component $component\n"
  unless $cv->{middle_axis} eq substr($component,1,1)
      && $cv->{right_axis} eq substr($component,2,1);
my @active_con = @{$cv->{input_brange}}[2,3];
if (defined $con_bands) {
  die "--con-bands must be FIRST:LAST (positive DFT band indices)\n"
    unless $con_bands =~ /\A([1-9][0-9]*):([1-9][0-9]*)\z/ && $1 <= $2;
  die "--con-bands must match the completed RbXc conduction range @active_con\n"
    unless $1 == $active_con[0] && $2 == $active_con[1];
}
die "RbXc ncon does not match its input conduction range\n"
  unless $cv->{ncon} == $active_con[1] - $active_con[0] + 1;

# PREP's binary con.u2.dat and parcksc files retain the original, wider band
# window. The solver's active window is taken from the completed CV manifest.
my $core_prep = read_json(File::Spec->catfile($ocean_run_root, 'PREP', 'prep.json'));
my $stored = $core_prep->{bse}{brange};
die "Core PREP lacks a four-element stored band range\n"
  unless ref($stored) eq 'ARRAY' && @$stored == 4
      && !grep { !defined($_) || !/\A[1-9][0-9]*\z/ } @$stored;
my @stored_con = @$stored[2,3];
die "Core PREP does not cover the RbXc conduction range @active_con\n"
  unless $active_con[0] >= $stored_con[0] && $active_con[1] <= $stored_con[1];
# OCEAN's core wvfcninfo and dipole readers select the first N conduction
# bands. A shifted lower bound needs additional OCEAN reader changes.
die "The active conduction FIRST must equal core PREP's stored FIRST ($stored_con[0])\n"
  unless $active_con[0] == $stored_con[0];
my $mesh = $core_prep->{bse}{kmesh};
my $xmesh = $core_prep->{bse}{xmesh};
my $nspin = $core_prep->{bse}{nspin};
die "Core PREP lacks valid kmesh, xmesh, or nspin\n"
  unless ref($mesh) eq 'ARRAY' && @$mesh == 3 &&
         ref($xmesh) eq 'ARRAY' && @$xmesh == 3 &&
         defined($nspin) && ($nspin == 1 || $nspin == 2) &&
         !grep { !defined($_) || !/\A[1-9][0-9]*\z/ } (@$mesh,@$xmesh);
my $nk = $mesh->[0] * $mesh->[1] * $mesh->[2];
die "Core PREP kmesh disagrees with RbXc\n"
  unless $nk == $cv->{nk} && join(',',@$mesh) eq join(',',@{$cv->{kmesh}});
my $expected_bytes = 16 * $nk * $nspin * ($stored_con[1]-$stored_con[0]+1)
  * $xmesh->[0] * $xmesh->[1] * $xmesh->[2];
die "Core PREP con.u2.dat size does not match its stored band/x/k/spin dimensions\n"
  unless $conduction_bytes == $expected_bytes;
my $valence_prep_json = File::Spec->catfile(dirname(dirname($valence_con)), 'prep.json');
if (basename(dirname($valence_con)) eq 'BSE' && -f $valence_prep_json) {
  my $valence_prep = read_json($valence_prep_json)->{bse};
  die "Valence PREP stored conduction range differs from core PREP\n"
    unless ref($valence_prep->{brange}) eq 'ARRAY' &&
      join(',',@{$valence_prep->{brange}}[2,3]) eq join(',',@stored_con);
  for my $key (qw(kmesh xmesh)) {
    die "Valence PREP $key differs from core PREP\n"
      unless ref($valence_prep->{$key}) eq 'ARRAY' &&
        join(',',@{$valence_prep->{$key}}) eq join(',',@{$core_prep->{bse}{$key}});
  }
  die "Valence PREP nspin differs from core PREP\n"
    unless ($valence_prep->{nspin} // -1) == $nspin;
}

# This revision deliberately operates where shg_rbxc_cv.x wrote its results.
# Requiring every RHS to be local prevents a supposedly self-contained stage
# from retaining hidden paths into another calculation directory.
my @staged_entries;
my %protected = ($cv_manifest => 1);
for my $entry (@{$cv->{entries}}) {
  my $source = $entry->{file};
  $source = File::Spec->catfile(dirname($cv_manifest), $source)
    unless File::Spec->file_name_is_absolute($source);
  $source = abs_path($source)
    or die "Cannot resolve CV vector $entry->{file}\n";
  die "RbXc vector is not in the current work directory: $source\n"
    unless dirname($source) eq $work_dir;
  die "Empty RbXc vector: $source\n" unless -s $source;
  $protected{$source} = 1;
  push @staged_entries, {
    %{$entry}, file => basename($source),
    target_ev => 2.0*$entry->{omega_ev},
  };
}

my $common_path = File::Spec->catfile(
  $ocean_run_root, 'Common', 'postDefaultsOceanDatafile'
);
my $common = read_json($common_path);
die "Core input must have calc.mode=xas\n"
  unless lc($common->{calc}{mode} // '') eq 'xas';
die "Selected edge does not exist in the core input\n"
  unless ref($common->{calc}{edges}) eq 'ARRAY'
      && $edge_index <= @{$common->{calc}{edges}};
$photon_energy = infer_photon_energy($common) unless defined $photon_energy;
die "--photon-energy must be positive\n" unless $photon_energy > 0.0;

if ($dry_run) {
  print "DRY RUN: in-place GV stage is valid\n";
  print "  OCEAN run root:   $ocean_run_root\n";
  print "  RbXc work dir:    $work_dir\n";
  print "  component:        $component\n";
  print "  final dipole:     d_g^", substr($component,0,1), "\n";
  print "  local CV vectors: ", scalar(@{$cv->{entries}}), "\n";
  print "  conduction bytes: $conduction_bytes\n";
  print "  stored core bands: @stored_con; active CV bands: @active_con\n";
  print "  content hash:     not checked\n";
  exit 0;
}

# stage_core_xas_cv.pl still expects a clean direct child of OCEAN_RUN. Build
# there first, validate it, and only then merge those generated inputs into the
# existing RbXc directory. CLEANUP removes this private staging directory.
my $clean_stage = tempdir(
  'SHG_GV_STAGE_XXXXXXXX', DIR => $ocean_run_root, CLEANUP => 1
);

write_file(
  File::Spec->catfile($clean_stage, 'shg_component'),
  join(' ', split(//, $component)) . "\n",
);

# The only core optical source needed by this pathway is the final d_g^a.
my $left_axis = substr($component,0,1);
my %polarization = (
  x => [1,0,0], y => [0,1,0], z => [0,0,1],
);
my %q_direction = (
  x => [0,1,0], y => [1,0,0], z => [1,0,0],
);
my $e = $polarization{$left_axis};
my $q = $q_direction{$left_axis};
my $photon = "dipole\ncartesian @{$e}\nend\n";
$photon .= "cartesian @{$q}\nend\n";
$photon .= sprintf("%.16g\n", $photon_energy);
write_file(File::Spec->catfile($clean_stage, 'photon1'), $photon);

my $helper = abs_path(File::Spec->catfile($FindBin::Bin, 'stage_core_xas_cv.pl'))
  or die "Cannot locate stage_core_xas_cv.pl\n";
if (defined $ocean_bin) {
  $ocean_bin = abs_path($ocean_bin)
    or die "Cannot resolve --ocean-bin\n";
  die "--ocean-bin does not contain meljtv.x\n"
    unless -x File::Spec->catfile($ocean_bin, 'meljtv.x');
  $ENV{OCEAN_BIN} = $ocean_bin;
}

local $ENV{SHG_EDGE_INDEX} = $edge_index;
local $ENV{SHG_COPY_ONLY} = 1;
my $old_cwd = getcwd();
chdir($clean_stage) or die "Cannot enter $clean_stage\n$!";
my $status = system($^X, $helper);
chdir($old_cwd) or die "Cannot return to $old_cwd\n$!";
die "Clean core-XAS staging failed\n" unless $status == 0;

# Retain the full PREP layout for the real-space Bloch-state reader, but put
# the core Hamiltonian, GMRES vector, and core dipole into the CV band window.
my $stage_brange = File::Spec->catfile($clean_stage, 'brange.ipt');
my @stage_bands = read_brange($stage_brange);
die "Staged core conduction range differs from PREP metadata\n"
  unless $stage_bands[2] == $stored_con[0] && $stage_bands[3] == $stored_con[1];
write_file(File::Spec->catfile($clean_stage, 'file_brange.ipt'),
           "$stored->[0] $stored->[1]\n$stored->[2] $stored->[3]\n");
write_file($stage_brange,
           "$stage_bands[0] $stage_bands[1]\n$active_con[0] $active_con[1]\n");
write_file(File::Spec->catfile($clean_stage, 'nbuse.ipt'),
           ($active_con[1]-$active_con[0]+1)."\n");
my $projectors_sliced = 0;
if ($active_con[1] < $stored_con[1]) {
  my @projectors = glob(File::Spec->catfile($clean_stage, 'parcksc.*'));
  die "A native-stream parcksc.* is needed to subset core projectors\n"
    unless @projectors;
  for my $projector (@projectors) {
    subset_core_projector($projector, $stored_con[1]-$stored_con[0]+1,
                          $active_con[1]-$active_con[0]+1, $nk, $nspin);
  }
  $projectors_sliced = scalar @projectors;
}

my $certificate = File::Spec->catfile($clean_stage, 'cv_core_basis_match.info');
my $verifier = abs_path(File::Spec->catfile($FindBin::Bin, 'verify_cv_core_basis.py'))
  or die "Cannot locate verify_cv_core_basis.py\n";
my $python = $ENV{PYTHON} // 'python3';
$status = system(
  $python, $verifier,
  '--cv-manifest', $cv_manifest,
  '--valence-con', $valence_con,
  '--core-con', File::Spec->catfile($clean_stage, 'con.u2.dat'),
  '--core-brange', File::Spec->catfile($clean_stage, 'brange.ipt'),
  '--core-kmesh', File::Spec->catfile($clean_stage, 'kmesh.ipt'),
  '--out', $certificate,
);
die "Conduction-basis certification failed\n" unless $status == 0;

# bse.in's core header is still used by OCEAN_multiplet. The second
# resolvent uses an independent, full bse2.in grid indexed by CV entry.
my $bse2 = make_bse2(File::Spec->catfile($clean_stage, 'bse.in'),
                     \@staged_entries, $bse2_input, $energy_range, $core_eta);
write_file(File::Spec->catfile($clean_stage, 'bse2.in'), $bse2);
my $targets = validate_bse2($bse2, scalar(@staged_entries));
for my $i (0..$#staged_entries) {
  $staged_entries[$i]{target_ev} = $targets->[$i];
}
my $list = "# index omega_eV target_Omega_eV RbXc_file\n";
for my $entry (@staged_entries) {
  $list .= sprintf("%d %.16g %.16g %s\n", $entry->{index},
                   $entry->{omega_ev}, $entry->{target_ev}, $entry->{file});
}
write_file(File::Spec->catfile($clean_stage, 'cv_vectors.list'), $list);

my $stage_info = {
  format_version => 3,
  second_gmres_input => 'bse2.in',
  index_pairing => 'full_bse2_grid_i_to_manifest_i',
  purpose => 'in-place standalone core-valence second-resolvent stage',
  sector => 'CORE_FROM_VALENCE',
  pathway => 'GV',
  component => [split(//, $component)],
  core_edge_index => 0+$edge_index,
  copied_inputs_only => JSON::PP::true,
  staged_in_existing_rbxc_directory => JSON::PP::true,
  rbxc_vectors_reused_in_place => JSON::PP::true,
  previous_core_response_used => JSON::PP::false,
  conduction_bytes => 0+$conduction_bytes,
  stored_core_conduction_brange => [@stored_con],
  active_core_conduction_brange => [@active_con],
  compact_core_projectors => 0+$projectors_sliced,
  conduction_content_hash_checked => JSON::PP::false,
  basis_check => 'file_size_conduction_brange_kmesh',
  cv_manifest => basename($cv_manifest),
  basis_certificate => basename($certificate),
  vectors => \@staged_entries,
};
write_file(
  File::Spec->catfile($clean_stage, 'shg_gv_stage.json'),
  JSON::PP->new->canonical(1)->pretty(1)->encode($stage_info),
);

my @stage_files = collect_stage_files($clean_stage);
my $inventory_name = 'shg_y1_cv_staged_files.list';
my $inventory = "# Files generated or replaced by prepare_y1_cv.pl\n";
$inventory .= join("\n", sort(@stage_files, $inventory_name)) . "\n";
write_file(File::Spec->catfile($clean_stage, $inventory_name), $inventory);
@stage_files = collect_stage_files($clean_stage);

# No generated core-XAS input may replace the completed CV manifest or an RbXc
# binary. Everything else in this list is an intentional Y1 staging output.
for my $relative (@stage_files) {
  my $destination = File::Spec->catfile($work_dir, split(m{/}, $relative));
  my $destination_abs = File::Spec->rel2abs($destination);
  die "Generated stage would replace protected RbXc data: $destination\n"
    if $protected{$destination_abs};
}
for my $relative (@stage_files) {
  my $source = File::Spec->catfile($clean_stage, split(m{/}, $relative));
  my $destination = File::Spec->catfile($work_dir, split(m{/}, $relative));
  atomic_copy($source, $destination);
}

print "In-place standalone GV stage complete: $work_dir\n";
print "  preserved CV manifest: ", basename($cv_manifest), "\n";
print "  preserved CV vectors:  ", scalar(@staged_entries), "\n";
print "  staged/replaced inputs: ", scalar(@stage_files), "\n";
print "  second GMRES input: bse2.in (full grid; entries are never renumbered)\n";
print "  final core dipole: d_g^$left_axis (generated by shg_y1_cv.x)\n";
print "  previous core response files used: no\n";
print "  core conduction bands: @active_con (stored PREP: @stored_con)\n";

sub read_brange {
  my ($path) = @_;
  open my $fh, '<', $path or die "Cannot read $path: $!\n";
  my @words = split ' ', do { local $/; <$fh> };
  close $fh;
  die "Invalid band range in $path\n"
    unless @words == 4 && !grep { !/\A[1-9][0-9]*\z/ } @words;
  return @words;
}

sub subset_core_projector {
  my ($path, $stored_ncon, $active_ncon, $nk, $nspin) = @_;
  open my $in, '<:raw', $path or die "Cannot open projector $path: $!\n";
  my $header = '';
  die "Short projector header: $path\n" unless read($in,$header,12) == 12;
  my ($nproj,$ntot,$file_spin) = unpack('l3',$header);
  die "Projector dimensions disagree with core PREP: $path\n"
    unless $nproj > 0 && $ntot == $stored_ncon*$nk && $file_spin == $nspin;
  my $band_bytes = 16*$nproj;
  die "Projector size disagrees with its header: $path\n"
    unless -s($path) == 12+$band_bytes*$stored_ncon*$nk*$nspin;
  my $temporary = "$path.subset.$$";
  open my $out, '>:raw', $temporary or die "Cannot create $temporary: $!\n";
  print {$out} pack('l3',$nproj,$active_ncon*$nk,$nspin)
    or die "Cannot write projector header: $!\n";
  my $chunk_bytes = $band_bytes*$active_ncon;
  for my $spin (0..$nspin-1) {
    for my $k (0..$nk-1) {
      my $offset = 12 + ($spin*$nk+$k)*$stored_ncon*$band_bytes;
      seek($in,$offset,0) or die "Cannot seek $path: $!\n";
      my $remaining = $chunk_bytes;
      while ($remaining > 0) {
        my $count = $remaining < 1048576 ? $remaining : 1048576;
        my $got = read($in,my $chunk,$count);
        die "Short projector data in $path\n" unless defined($got) && $got > 0;
        print {$out} $chunk or die "Cannot write $temporary: $!\n";
        $remaining -= $got;
      }
    }
  }
  close $in or die "Cannot close $path: $!\n";
  close $out or die "Cannot close $temporary: $!\n";
  rename $temporary,$path or die "Cannot replace $path: $!\n";
}

sub read_cv_manifest {
  my ($path) = @_;
  open my $fh, '<', $path or die "Cannot open $path\n$!";
  my %data = (entries => []);
  while (my $line = <$fh>) {
    my @w = split ' ', $line;
    next unless @w;
    $data{status} = $w[1] if $w[0] eq 'status';
    $data{sector} = $w[1] if $w[0] eq 'sector';
    $data{completion} = $w[1] if $w[0] eq 'completion';
    $data{middle_axis} = lc($w[1]) if $w[0] eq 'middle_axis';
    $data{right_axis} = lc($w[1]) if $w[0] eq 'right_axis';
    $data{input_brange} = [@w[1..4]] if $w[0] eq 'input_brange' && @w == 5;
    $data{kmesh} = [@w[1..3]] if $w[0] eq 'kmesh' && @w == 4;
    ($data{ncon},$data{nk}) = @w[4,5]
      if $w[0] eq 'dimensions' && @w == 7 && $w[1] eq 'ncon' && $w[2] eq 'nk';
    if ($w[0] eq 'RbXc') {
      die "Malformed RbXc record in $path\n" unless @w == 6;
      push @{$data{entries}}, {
        index => 0+$w[1], omega_ha => 0.0+$w[2],
        omega_ev => 0.0+$w[3], norm => 0.0+$w[4], file => $w[5],
      };
    }
  }
  close $fh or die "Cannot close $path\n$!";
  die "CV manifest is not complete\n"
    unless ($data{status} // '') eq 'RUNNING'
        && ($data{sector} // '') eq 'CORE_FROM_VALENCE'
        && ($data{completion} // '') eq 'RBXC_CV_COMPLETE'
        && @{$data{entries}};
  die "CV manifest lacks valid band/k dimensions\n"
    unless ref($data{input_brange}) eq 'ARRAY' &&
           ref($data{kmesh}) eq 'ARRAY' &&
           defined($data{ncon}) && defined($data{nk}) &&
           !grep { !defined($_) || !/\A[1-9][0-9]*\z/ }
             (@{$data{input_brange}},@{$data{kmesh}},$data{ncon},$data{nk});
  @{$data{entries}} = sort { $a->{index} <=> $b->{index} } @{$data{entries}};
  for my $i (0..$#{$data{entries}}) {
    die "Manifest indices must be contiguous from 1\n" unless $data{entries}[$i]{index} == $i+1;
    die "Nonfinite manifest energy\n"
      unless isfinite($data{entries}[$i]{omega_ev}) && isfinite($data{entries}[$i]{omega_ha});
  }
  return \%data;
}


sub read_text {
  my ($path) = @_;
  open my $fh, '<', $path or die "Cannot read $path: $!\n";
  local $/; my $text = <$fh>; close $fh or die "Cannot close $path: $!\n";
  return $text;
}

sub numeric {
  my ($value) = @_;
  return defined($value) && looks_like_number($value) && isfinite(0+$value);
}

sub expand_target_loop {
  my ($range, $expected_count) = @_;
  my @v = split /:/, $range;
  die "--energy-range must be MIN:MAX:STEP with finite values\n"
    unless @v == 3 && !grep { !numeric($_) } @v;
  my ($first,$last,$step) = map { 0+$_ } @v;
  die "Invalid second-resolvent loop direction/increment\n"
    unless $step != 0 && ($last-$first)/$step >= 0;
  # Match SHG_gmres_setup's OCEAN loop convention exactly.
  my $span = ($last-$first)/$step;
  die "Second-resolvent loop is too large\n" unless isfinite($span) && $span < 2147483645;
  my $count = 1 + POSIX::ceil($span-0.0001);
  die "bse2.in must contain the full manifest grid ($expected_count entries; got $count)\n"
    unless $count == $expected_count;
  return [map { $first+$_*$step } 0..$count-1];
}

sub validate_bse2 {
  my ($text, $expected_count) = @_;
  my @line = split /\n/, $text;
  die "bse2.in needs the full OCEAN bse.in format\n" unless @line >= 6;
  die "bse2.in solver must be inv\n" unless $line[2] =~ /^\s*inv\s*$/;
  my @par = split ' ', $line[3];
  die "Invalid GMRES parameter line in bse2.in\n"
    unless @par >= 4 && !grep { !numeric($_) } @par[0..3];
  die "Invalid GMRES depth/preconditioner/convergence in bse2.in\n"
    unless $par[0] =~ /^\d+$/ && $par[0] > 0 && $par[3] > 0;
  my $style = $line[4]; $style =~ s/^\s+|\s+$//g;
  my $targets;
  if ($style eq 'loop') {
    my @v = split ' ', $line[5];
    die "Invalid bse2.in loop\n" unless @v == 3;
    $targets = expand_target_loop(join(':',@v), $expected_count);
  } elsif ($style eq 'list') {
    die "bse2.in must contain $expected_count targets\n"
      unless $line[5] =~ /^\s*\d+\s*$/ && 0+$line[5] == $expected_count && @line >= 6+$expected_count;
    $targets = [@line[6..5+$expected_count]];
    die "Nonfinite bse2.in list target\n" if grep { !numeric($_) } @$targets;
    @$targets = map { 0+$_ } @$targets;
  } else {
    die "bse2.in energy style must be list or loop\n";
  }
  return $targets;
}

sub make_bse2 {
  my ($core_bse, $entries, $supplied, $range, $eta) = @_;
  my @base = split /\n/, read_text($core_bse);
  die "Core staging did not produce a GMRES bse.in\n"
    unless @base >= 4 && $base[2] =~ /^\s*inv\s*$/;
  my @line;
  if (defined $supplied) {
    @line = split /\n/, read_text($supplied);
    validate_bse2(join("\n",@line)."\n", scalar @$entries);
    # Use the actual staged core header, retain the supplied GMRES settings/grid.
    @line[0,1] = @base[0,1];
  } else {
    @line = @base[0..3];
    my @target = map { 2*$_->{omega_ev} } @$entries;
    if (defined $range) {
      @target = @{expand_target_loop($range, scalar @$entries)};
      push @line, 'loop', join(' ',split /:/,$range);
    } else {
      # A list retains the manifest energies without rounding a reconstructed step.
      push @line, 'list', scalar(@target), map { sprintf('%.16g',$_) } @target;
    }
  }
  if (defined $eta) {
    my @par = split ' ', $line[3]; $par[1] = sprintf('%.16g',$eta);
    $line[3] = join(' ',@par);
  }
  my $text = join("\n",@line)."\n";
  validate_bse2($text, scalar @$entries);
  return $text;
}

sub read_json {
  my ($path) = @_;
  open my $fh, '<', $path or die "Cannot open $path\n$!";
  local $/ = undef;
  my $decoded = JSON::PP->new->decode(<$fh>);
  close $fh or die "Cannot close $path\n$!";
  return $decoded;
}

sub infer_photon_energy {
  my ($data) = @_;
  my $core = $data->{bse}{core} // {};
  my $gmres = $core->{gmres} // {};
  return $gmres->{elist}[0]
    if ref($gmres->{elist}) eq 'ARRAY' && @{$gmres->{elist}};
  return $gmres->{erange}[0]
    if ref($gmres->{erange}) eq 'ARRAY' && @{$gmres->{erange}};
  return $core->{plot}{range}[0]
    if ref($core->{plot}{range}) eq 'ARRAY' && @{$core->{plot}{range}};
  return 1.0;
}

sub write_file {
  my ($path, $content) = @_;
  open my $fh, '>', $path or die "Cannot create $path\n$!";
  print {$fh} $content;
  close $fh or die "Cannot close $path\n$!";
}

sub collect_stage_files {
  my ($root) = @_;
  my @files;
  find({
    no_chdir => 1,
    wanted => sub {
      die "Clean stage unexpectedly contains a symbolic link: $File::Find::name\n"
        if -l $File::Find::name;
      return unless -f $File::Find::name;
      my $relative = File::Spec->abs2rel($File::Find::name, $root);
      $relative =~ s{\\}{/}g;
      die "Invalid clean-stage path: $relative\n" if $relative =~ m{(?:^|/)\.\.(?:/|$)};
      push @files, $relative;
    },
  }, $root);
  return sort @files;
}

sub atomic_copy {
  my ($source, $destination) = @_;
  my $parent = dirname($destination);
  make_path($parent) unless -d $parent;
  my ($tmp_fh, $tmp_path) = tempfile('.prepare_y1_cv.XXXXXX', DIR => $parent);
  close $tmp_fh or die "Cannot close temporary output $tmp_path\n$!";
  copy($source, $tmp_path)
    or die "Failed to copy $source to $tmp_path\n$!";
  my $mode = (stat($source))[2] & 07777;
  chmod($mode, $tmp_path)
    or die "Failed to preserve permissions on $tmp_path\n$!";
  rename($tmp_path, $destination)
    or die "Failed to replace $destination\n$!";
}

sub usage {
  my ($code) = @_;
  print <<'USAGE';
Usage (run in the completed RbXc directory):
  prepare_y1_cv.pl --ocean-run-root OCEAN_RUN --component ABC \
    --valence-con PREP_VAL/BSE/con.u2.dat [options]

The current directory must contain shg_rbxc_cv.info and every RbXc binary
listed by that completed manifest. The script stages/replaces the core-XAS
inputs in this same directory; the manifest and RbXc binaries are preserved.

Options:
  --work-dir DIR      existing RbXc directory (default: current directory)
  --target DIR        deprecated alias for --work-dir
  --cv-manifest FILE  completed manifest (default: WORK_DIR/shg_rbxc_cv.info)
  --edge-index N      selected configured core edge (default: 1)
  --con-bands F:L     active core bands; defaults to the CV manifest's range
                      and must match it exactly (stored core PREP may be wider)
  --photon-energy EV  spectfile value for the newly generated final d_g^a
  --bse2-input FILE   supplied full bse.in-format second-solver input
  --energy-range A:B:S second-solver Omega loop in eV (full manifest length)
                      default: list of twice each manifest omega
  --core-eta EV       override second-solver broadening; default core setting
  --ocean-bin DIR     directory containing the original meljtv.x
  --dry-run           validate without replacing any files
  --help              show this message

OCEAN_RUN is the OCEAN computation directory containing Common, DFT, OPF,
PREP, and SCREEN. A clean core-XAS stage is built temporarily beneath it and
then merged into WORK_DIR. bse2.in is generated for the batch second solver;
bse.in remains for the OCEAN core multiplet header. Existing outputs are not deleted.
Run shg_y1_cv.x --component ABC [--start-index N]; retain the FULL bse2.in grid.
The stored PREP conduction FIRST must equal the CV FIRST; its LAST may be
higher. Native-stream parcksc projectors are subset for the core dipole reader.
USAGE
  exit $code;
}
