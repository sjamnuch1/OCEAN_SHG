#!/usr/bin/env perl
# Run only the Common/input-processing portion of OCEAN for calc.mode=val.
use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use Getopt::Long qw(GetOptions);
use JSON::PP;

my ($input,$root,$output,$ocean_bin,$help);
GetOptions(
  'input=s'     => \$input,
  'root=s'      => \$root,
  'output=s'    => \$output,
  'ocean-bin=s' => \$ocean_bin,
  'help|h'      => \$help,
) or usage(2);
usage(0) if $help;
$input //= shift @ARGV;
die "Unexpected positional arguments\n" if @ARGV;
die "An OCEAN input file is required (--input FILE)\n" unless defined $input;

$root=abs_path($root // getcwd()) or die "Cannot resolve run root\n";
$input=File::Spec->rel2abs($input,getcwd());
die "Input is not a regular file: $input\n" unless -f $input;
$input=abs_path($input) or die "Cannot resolve input file\n";
$ocean_bin=abs_path($ocean_bin // $ENV{OCEAN_BIN} // '')
  or die "Set OCEAN_BIN or supply --ocean-bin\n";
$output //= File::Spec->catdir($root,'Common_VAL');
$output=File::Spec->rel2abs($output,$root);
die "Output must be a new direct child of the run root\n"
  unless dirname_abs($output) eq $root;
die "Output already exists; refusing to overwrite $output\n"
  if -e $output || -l $output;

my @required=qw(defaults.pl structure.pl edges.pl photon.pl psp_parser.pl);
my $have_json=-f File::Spec->catfile($ocean_bin,'parseJSON.pl');
if($have_json) {
  push @required,qw(parseJSON.pl oparse.json oparse.type.json);
} else {
  push @required,qw(parse.pl oparse.h);
}
push @required,'extractPsp.pl';
for my $file (@required) {
  die "Missing OCEAN Common helper $ocean_bin/$file\n"
    unless -f File::Spec->catfile($ocean_bin,$file);
}

make_path($output) or die "Cannot create $output\n";
my $input_name=basename($input);
copy($input,File::Spec->catfile($output,$input_name))
  or die "Cannot copy input into $output: $!\n";

$ENV{OCEAN_BIN}=$ocean_bin;
$ENV{OCEAN_WORKDIR}=$root;
my $version=File::Spec->catfile($ocean_bin,'Version');
if(open my $vf,'<',$version) {
  chomp($ENV{OCEAN_VERSION}=<$vf> // 'NaN');
  chomp($ENV{OCEAN_VERSION_HASH}=<$vf> // '');
  close $vf;
} else {
  $ENV{OCEAN_VERSION}='NaN';
  $ENV{OCEAN_VERSION_HASH}='';
}

chdir $output or die "Cannot enter $output: $!\n";
if($have_json) {
  run_perl('parseJSON.pl',$input_name,
    File::Spec->catfile($ocean_bin,'oparse.json'),
    File::Spec->catfile($ocean_bin,'oparse.type.json'));
} else {
  run_perl('parse.pl','--all',$input_name,
    File::Spec->catfile($ocean_bin,'oparse.h'));
}
run_perl('extractPsp.pl');
my $psp_status=run_perl_allow_failure('psp_parser.pl');
warn "WARNING: psp_parser.pl failed with exit code $psp_status; continuing as original ocean.pl does\n"
  if $psp_status;
run_perl('defaults.pl');
run_perl('structure.pl');
run_perl('edges.pl');
run_perl('photon.pl');

my $json_path=File::Spec->catfile($output,'postDefaultsOceanDatafile');
die "Common processing did not create postDefaultsOceanDatafile\n" unless -s $json_path;
open my $jf,'<',$json_path or die "Cannot read $json_path: $!\n";
my $data=eval { JSON::PP->new->decode(do { local $/; <$jf> }) };
close $jf;
die "Invalid postDefaultsOceanDatafile: $@\n" if $@ || ref($data) ne 'HASH';
my $mode=lc($data->{calc}{mode} // '');
die "The processed input has calc.mode='$mode', not 'val'. Output left for inspection at $output\n"
  unless $mode eq 'val';
die "Common output is missing bse.kmesh\n"
  unless ref($data->{bse}{kmesh}) eq 'ARRAY' && @{$data->{bse}{kmesh}}==3;
die "Common output is missing calc.photon_q\n"
  unless ref($data->{calc}{photon_q}) eq 'ARRAY' && @{$data->{calc}{photon_q}}==3;
my $q2=0; $q2 += $_*$_ for @{$data->{calc}{photon_q}};
warn "WARNING: photon_q is zero. Normal VAL PREP is possible, but shg_xc_minus_valence.x requires a finite small q.\n"
  unless $q2>0;

open my $sf,'>',File::Spec->catfile($output,'common_val.complete')
  or die "Cannot create completion marker: $!\n";
print {$sf} "COMMON_VAL_COMPLETE\n";
print {$sf} "input $input\n";
print {$sf} "calc.mode val\n";
close $sf or die "Cannot close completion marker: $!\n";
print "Completed VAL Common stage in $output\n";
print "Use with prep_val.pl --common $json_path\n";

sub run_perl {
  my ($script,@args)=@_;
  my $status=run_perl_allow_failure($script,@args);
  die "$script failed with exit code $status\n" if $status;
}
sub run_perl_allow_failure {
  my ($script,@args)=@_;
  my $path=File::Spec->catfile($ocean_bin,$script);
  system {$^X} $^X,$path,@args;
  die "Failed to execute $script: $!\n" if $?==-1;
  die "$script died with signal ".($?&127)."\n" if $?&127;
  return $?>>8;
}
sub dirname_abs {
  my ($path)=@_;
  require File::Basename;
  return abs_path(File::Basename::dirname($path)) // '';
}
sub usage {
  my ($code)=@_;
  print <<'USAGE';
Usage: perl ocean_common_val.pl --input ocean.in [options]
  --root DIR       Run root (default current directory)
  --output DIR     Fresh direct child output (default ROOT/Common_VAL)
  --ocean-bin DIR  Directory containing installed OCEAN Common scripts
  --help           Show this message

Runs parseJSON/parse, extractPsp, psp_parser, defaults, structure, edges,
and photon only. Requires the processed calc.mode to be val. It never starts
OPF, DFT, PREP, SCREEN, or CNBSE and never modifies ROOT/Common.
USAGE
  exit $code;
}

