#!/usr/bin/env perl
# Create a fresh, copy-only core-XAS directory for the GV second resolvent.

use strict;
use warnings;

use Cwd qw(abs_path getcwd);
use Digest::SHA qw(sha256_hex);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Find qw(find);
use File::Path qw(make_path);
use File::Spec;
use FindBin;
use Getopt::Long qw(GetOptions);
use JSON::PP;

my ($ocean_run_root, $target, $component, $cv_manifest, $valence_con);
my ($ocean_bin, $photon_energy);
my $edge_index = 1;
my $dry_run = 0;
my $help = 0;

GetOptions(
  'ocean-run-root=s' => \$ocean_run_root,
  'target=s'        => \$target,
  'component=s'     => \$component,
  'edge-index=i'    => \$edge_index,
  'cv-manifest=s'   => \$cv_manifest,
  'valence-con=s'   => \$valence_con,
  'photon-energy=f' => \$photon_energy,
  'ocean-bin=s'     => \$ocean_bin,
  'dry-run!'        => \$dry_run,
  'help|h'          => \$help,
) or usage(1);
usage(0) if $help;

for my $required (
  ['--ocean-run-root', $ocean_run_root], ['--target', $target],
  ['--component', $component], ['--cv-manifest', $cv_manifest],
  ['--valence-con', $valence_con],
) {
  die "$required->[0] is required\n" unless defined $required->[1];
}

$ocean_run_root = abs_path($ocean_run_root)
  or die "Cannot resolve OCEAN computation directory\n";
$cv_manifest = abs_path($cv_manifest)
  or die "Cannot resolve CV manifest\n";
$valence_con = abs_path($valence_con)
  or die "Cannot resolve valence con.u2.dat\n";
$target = File::Spec->rel2abs($target);
$component = lc($component);
$component =~ s/[\s,]+//g;
die "--component must contain exactly three x/y/z labels\n"
  unless $component =~ /\A[xyz]{3}\z/;
die "--edge-index must be a positive integer\n"
  unless $edge_index =~ /^\d+$/ && $edge_index >= 1;
die "--target must be a direct child of --ocean-run-root\n"
  unless dirname($target) eq $ocean_run_root;
die "Target already exists; choose a fresh directory: $target\n"
  if -e $target || -l $target;

for my $stage (qw(Common DFT OPF PREP SCREEN)) {
  my $path = File::Spec->catdir($ocean_run_root, $stage);
  die "Missing completed core input stage: $path\n" unless -d $path;
}
my $core_con = File::Spec->catfile($ocean_run_root, 'PREP', 'BSE', 'con.u2.dat');
die "Missing core conduction wavefunctions: $core_con\n" unless -f $core_con;
die "Empty valence conduction wavefunction file\n" unless -s $valence_con;
die "Core and valence con.u2.dat sizes differ\n"
  unless -s($core_con) == -s($valence_con);
my $core_hash = hash_file($core_con);
my $valence_hash = hash_file($valence_con);
die "Core and valence con.u2.dat differ; an explicit basis transformation is required\n"
  unless $core_hash eq $valence_hash;

my $cv = read_cv_manifest($cv_manifest);
die "CV manifest middle/right axes do not match component $component\n"
  unless $cv->{middle_axis} eq substr($component,1,1)
      && $cv->{right_axis} eq substr($component,2,1);

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
  print "DRY RUN: fresh GV stage is valid\n";
  print "  OCEAN run root:  $ocean_run_root\n";
  print "  target:          $target\n";
  print "  component:       $component\n";
  print "  final dipole:    d_g^", substr($component,0,1), "\n";
  print "  CV vectors:      ", scalar(@{$cv->{entries}}), "\n";
  print "  conduction hash: $core_hash\n";
  exit 0;
}

make_path($target) or die "Failed to create $target\n";
write_new(
  File::Spec->catfile($target, 'shg_component'),
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
write_new(File::Spec->catfile($target, 'photon1'), $photon);

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
chdir($target) or die "Cannot enter $target\n$!";
my $status = system($^X, $helper);
chdir($old_cwd) or die "Cannot return to $old_cwd\n$!";
die "Fresh core-XAS staging failed\n" unless $status == 0;

# Copy the complete CV RHS set and its manifest. Runtime uses no valence path.
my $manifest_copy = File::Spec->catfile($target, basename($cv_manifest));
copy_new($cv_manifest, $manifest_copy);
my @staged_entries;
for my $entry (@{$cv->{entries}}) {
  my $source = $entry->{file};
  $source = File::Spec->catfile(dirname($cv_manifest), $source)
    unless File::Spec->file_name_is_absolute($source);
  $source = abs_path($source)
    or die "Cannot resolve CV vector $entry->{file}\n";
  my $destination = File::Spec->catfile($target, basename($source));
  copy_new($source, $destination);
  push @staged_entries, {
    %{$entry}, file => basename($destination),
    target_ev => 2.0*$entry->{omega_ev},
  };
}
my $valence_copy = File::Spec->catfile($target, 'valence_con.u2.dat');
copy_new($valence_con, $valence_copy);

my $certificate = File::Spec->catfile($target, 'cv_core_basis_match.info');
my $verifier = abs_path(File::Spec->catfile($FindBin::Bin, 'verify_cv_core_basis.py'))
  or die "Cannot locate verify_cv_core_basis.py\n";
my $python = $ENV{PYTHON} // 'python3';
$status = system(
  $python, $verifier,
  '--cv-manifest', $manifest_copy,
  '--valence-con', $valence_copy,
  '--core-con', File::Spec->catfile($target, 'con.u2.dat'),
  '--core-brange', File::Spec->catfile($target, 'brange.ipt'),
  '--core-kmesh', File::Spec->catfile($target, 'kmesh.ipt'),
  '--out', $certificate,
);
die "Conduction-basis certification failed\n" unless $status == 0;

my $list = "# index omega_eV target_2omega_eV RbXc_file\n";
for my $entry (@staged_entries) {
  $list .= sprintf("%d %.16g %.16g %s\n", $entry->{index},
                   $entry->{omega_ev}, $entry->{target_ev}, $entry->{file});
}
write_new(File::Spec->catfile($target, 'cv_vectors.list'), $list);

my $stage_info = {
  format_version => 1,
  purpose => 'fresh standalone core-valence second-resolvent stage',
  sector => 'CORE_FROM_VALENCE',
  pathway => 'GV',
  component => [split(//, $component)],
  core_edge_index => 0+$edge_index,
  copied_inputs_only => JSON::PP::true,
  previous_core_response_used => JSON::PP::false,
  conduction_sha256 => $core_hash,
  cv_manifest => basename($manifest_copy),
  basis_certificate => basename($certificate),
  vectors => \@staged_entries,
};
write_new(
  File::Spec->catfile($target, 'shg_gv_stage.json'),
  JSON::PP->new->canonical(1)->pretty(1)->encode($stage_info),
);

my @links;
find(sub { push @links, $File::Find::name if -l $_ }, $target);
die "Fresh stage unexpectedly contains symbolic links: @links\n" if @links;

print "Fresh standalone GV stage complete: $target\n";
print "  copied CV vectors: ", scalar(@staged_entries), "\n";
print "  final core dipole: d_g^$left_axis (generated by shg_y1_cv.x)\n";
print "  previous core response files used: no\n";

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
  return \%data;
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

sub hash_file {
  my ($path) = @_;
  open my $fh, '<:raw', $path or die "Cannot open $path\n$!";
  my $sha = Digest::SHA->new(256);
  $sha->addfile($fh);
  close $fh or die "Cannot close $path\n$!";
  return $sha->hexdigest;
}

sub copy_new {
  my ($source, $destination) = @_;
  die "Refusing to replace $destination\n" if -e $destination || -l $destination;
  copy($source, $destination)
    or die "Failed to copy $source to $destination\n$!";
}

sub write_new {
  my ($path, $content) = @_;
  die "Refusing to replace $path\n" if -e $path || -l $path;
  open my $fh, '>', $path or die "Cannot create $path\n$!";
  print {$fh} $content;
  close $fh or die "Cannot close $path\n$!";
}

sub usage {
  my ($code) = @_;
  print <<'USAGE';
Usage:
  prepare_y1_cv.pl --ocean-run-root OCEAN_RUN --target FRESH_DIR \
    --component ABC --edge-index N --cv-manifest shg_rbxc_cv.info \
    --valence-con PREP_VAL/BSE/con.u2.dat [options]

OCEAN_RUN is the OCEAN computation directory containing Common, DFT, OPF,
PREP, and SCREEN.  The target must not exist and must be its direct child.
All runtime inputs are copied; no symlink or previous core response vector is
used.  Every RbXc dump listed by the completed CV manifest is copied.

Options:
  --photon-energy EV  spectfile value for the newly generated final d_g^a
  --ocean-bin DIR     directory containing the original meljtv.x
  --dry-run           validate without creating the target
  --help              show this message
USAGE
  exit $code;
}
