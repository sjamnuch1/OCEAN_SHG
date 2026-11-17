#!/usr/bin/env perl
# Copyright (C) 2026 OCEAN collaboration
#
# Prepare one core/XAS right-hand side for the standalone SHG exciton-vector
# test. The original OCEAN source and stage outputs are never modified.

use strict;
use warnings;

use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use File::Spec;
use FindBin;
use Getopt::Long qw(GetOptions);
use JSON::PP;

my $run_root;
my $target;
my $component;
my $ocean_bin;
my $photon_energy;
my $edge_index = 1;
my $dry_run = 0;
my $force = 0;
my $help = 0;

GetOptions(
  'root=s'          => \$run_root,
  'target=s'        => \$target,
  'component=s'     => \$component,
  'edge-index=i'    => \$edge_index,
  'photon-energy=f' => \$photon_energy,
  'ocean-bin=s'     => \$ocean_bin,
  'dry-run!'        => \$dry_run,
  'force!'          => \$force,
  'help|h'          => \$help,
) or usage(1);
usage(0) if $help;

die "--component ABC is required\n" unless defined $component;
$component = lc($component);
$component =~ s/[\s,]+//g;
die "--component must contain exactly three x/y/z labels\n"
  unless $component =~ /\A[xyz]{3}\z/;
die "--edge-index must be at least 1\n" unless $edge_index >= 1;
die "--photon-energy must be positive\n"
  if defined($photon_energy) && $photon_energy <= 0.0;

my $cwd = getcwd();
if (!defined $run_root) {
  if (basename($cwd) eq 'SHG_core'
      && -d File::Spec->catdir($cwd, File::Spec->updir(), 'Common')) {
    $run_root = File::Spec->catdir($cwd, File::Spec->updir());
    $target = $cwd unless defined $target;
  } else {
    $run_root = $cwd;
  }
}
$run_root = abs_path($run_root)
  or die "Unable to resolve OCEAN run root '$run_root'\n";
$target = File::Spec->catdir($run_root, 'SHG_core') unless defined $target;
$target = File::Spec->rel2abs($target);

die "The SHG target must be a direct child of the completed run root\n"
  unless dirname($target) eq $run_root;
die "The SHG target must not be the completed run root\n"
  if $target eq $run_root;

my %stage_dir = map {
  $_ => File::Spec->catdir($run_root, $_)
} qw(Common DFT OPF PREP SCREEN);
for my $stage (qw(Common DFT OPF PREP SCREEN)) {
  die "Missing completed-stage directory $stage_dir{$stage}\n"
    unless -d $stage_dir{$stage};
}

my $common_path = File::Spec->catfile($stage_dir{Common}, 'postDefaultsOceanDatafile');
my $dft_path = File::Spec->catfile($stage_dir{DFT}, 'dft.json');
my $opf_path = File::Spec->catfile($stage_dir{OPF}, 'opf.json');
my $prep_path = File::Spec->catfile($stage_dir{PREP}, 'prep.json');
my $screen_path = File::Spec->catfile($stage_dir{SCREEN}, 'screen.json');

my $common = read_json($common_path);
my $dft = read_json($dft_path);
read_json($opf_path);
my $prep = read_json($prep_path);
my $screen = read_json($screen_path);

die "The standalone core exciton test requires calc.mode = xas or rxs\n"
#  unless lc($common->{calc}{mode} // '') eq 'xas';
  unless (lc($common->{calc}{mode} // '') =~ /^(?:xas|rxs)$/) ;
assert_true($dft, [qw(scf complete)], 'BSE wavefunction run');
assert_true($dft, [qw(density complete)], 'scf wavefunction run');
assert_true($dft, [qw(potential complete)], 'Potential postproccessing run');
assert_true($prep, [qw(bse complete)], 'PREP BSE processing');
if (exists $screen->{combine}{complete}) {
  assert_true($screen, [qw(combine complete)], 'SCREEN combination');
} else {
  assert_true($screen, [qw(screen complete)], 'SCREEN calculation');
}

die "PREP/BSE output directory is missing\n"
  unless -d File::Spec->catdir($stage_dir{PREP}, 'BSE');
die "Configured edge $edge_index does not exist\n"
  unless ref($common->{calc}{edges}) eq 'ARRAY'
      && $edge_index <= @{$common->{calc}{edges}};

$photon_energy = infer_photon_energy($common) unless defined $photon_energy;
my $right_axis = substr($component, 2, 1);
my %polarization = (
  x => [1, 0, 0],
  y => [0, 1, 0],
  z => [0, 0, 1],
);
my %q_direction = (
  x => [0, 1, 0],
  y => [1, 0, 0],
  z => [1, 0, 0],
);

my $helper = abs_path(File::Spec->catfile($FindBin::Bin, 'stage_core_xas.pl'))
  or die "Unable to locate stage_core_xas.pl next to prepare_shg.pl\n";

if (defined $ocean_bin) {
  $ocean_bin = abs_path($ocean_bin)
    or die "Unable to resolve --ocean-bin '$ocean_bin'\n";
  die "--ocean-bin does not contain meljtv.x\n"
    unless -x File::Spec->catfile($ocean_bin, 'meljtv.x');
  $ENV{OCEAN_BIN} = $ocean_bin;
}

if ($dry_run) {
  print "DRY RUN: core/XAS standalone staging is valid\n";
  print "  run root:       $run_root\n";
  print "  target:         $target\n";
  print "  component:      $component\n";
  print "  right axis:     $right_axis\n";
  print "  edge index:     $edge_index\n";
  print "  photon energy:  $photon_energy eV\n";
  exit 0;
}

if (-e $target) {
  die "Target already exists; use --force to refresh $target\n" unless $force;
  die "SHG target is not a directory: $target\n" unless -d $target;
} else {
  make_path($target) or die "Failed to create $target\n";
}

write_generated_file(
  File::Spec->catfile($target, 'shg_component'),
  join(' ', split(//, $component)) . "\n",
);

my $e = $polarization{$right_axis};
my $q = $q_direction{$right_axis};
my $photon_text = "dipole\n";
$photon_text .= "cartesian @{$e}\nend\n";
$photon_text .= "cartesian @{$q}\nend\n";
$photon_text .= sprintf("%.16g\n", $photon_energy);
write_generated_file(File::Spec->catfile($target, 'photon1'), $photon_text);

local $ENV{SHG_EDGE_INDEX} = $edge_index;
my $old_cwd = getcwd();
chdir($target) or die "Failed to enter $target\n$!";
my $status = system($^X, $helper);
chdir($old_cwd) or die "Failed to return to $old_cwd\n$!";
die "Core-XAS staging helper failed\n" unless $status == 0;

my $runlist_path = File::Spec->catfile($target, 'runlist');
open my $runlist, '<', $runlist_path or die "Failed to open generated runlist\n$!";
my @runlist_lines = grep { /\S/ } <$runlist>;
close $runlist;
die "Generated runlist is not exactly one XAS run\n"
  unless @runlist_lines == 2
      && $runlist_lines[0] =~ /^\s*1\s*$/
      && $runlist_lines[1] =~ /\bXAS\s*$/;

my $manifest = {
  format_version => 2,
  purpose => 'standalone OCEAN core-SHG exciton-vector input staging',
  sector => 'CORE_XAS',
  component => [split(//, $component)],
  right_axis => $right_axis,
  edge_index => $edge_index,
  photon_energy_ev => 0.0 + $photon_energy,
  run_root => File::Spec->abs2rel($run_root, $target),
};
write_generated_file(
  File::Spec->catfile($target, 'shg_stage.json'),
  JSON::PP->new->canonical(1)->pretty(1)->encode($manifest),
);

print "Standalone core-SHG staging complete: $target\n";
print "  component:     $component\n";
print "  core source:   d_core^$right_axis\n";
print "  selected edge: $edge_index\n";

sub read_json {
  my ($path) = @_;
  die "Missing JSON input $path\n" unless -e $path;
  open my $fh, '<', $path or die "Failed to open $path\n$!";
  local $/ = undef;
  my $text = <$fh>;
  close $fh or die "Failed to close $path\n$!";
  my $decoded = eval { JSON::PP->new->decode($text) };
  die "Failed to parse $path as JSON\n$@" if $@;
  return $decoded;
}

sub assert_true {
  my ($data, $path, $description) = @_;
  my $value = $data;
  for my $key (@{$path}) {
    die "$description is not marked complete (missing " . join('.', @{$path}) . ")\n"
      unless ref($value) eq 'HASH' && exists $value->{$key};
    $value = $value->{$key};
  }
  die "$description is not marked complete\n" unless $value;
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

sub write_generated_file {
  my ($path, $content) = @_;
  if (-e $path || -l $path) {
    die "Refusing to replace directory $path\n" if -d $path && !-l $path;
    unlink($path) or die "Failed to replace $path\n$!";
  }
  open my $out, '>', $path or die "Failed to write $path\n$!";
  print {$out} $content;
  close $out or die "Failed to close $path\n$!";
}

sub usage {
  my ($exit_code) = @_;
  print <<'USAGE';
Usage:
  prepare_shg.pl --component ABC [options]

Run this from the root of a completed core/XAS OCEAN calculation after
DFT, OPF, PREP, and SCREEN. It creates ./SHG_core and stages exactly one XAS edge
and one Cartesian dipole source for the standalone core-exciton GMRES test.

Options:
  --component ABC      SHG tensor component; the third label selects d_core^c
  --edge-index N       1-based index in calc.edges (default: 1)
  --photon-energy EV   spectfile energy (inferred when omitted; dipole source
                       is independent of this value)
  --root DIR           completed OCEAN run root (default: current directory)
  --target DIR         direct child target (default: RUN_ROOT/SHG)
  --ocean-bin DIR      installed original OCEAN bin directory; otherwise the
                       source-tree meljtv.x is used
  --dry-run            validate and print the selected inputs only
  --force              refresh an existing target directory
  --help, -h           show this help
USAGE
  exit $exit_code;
}
