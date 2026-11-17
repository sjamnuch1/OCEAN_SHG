#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

require JSON::PP;
use JSON::PP;

use Cwd 'abs_path';
use Cwd;
use File::Spec::Functions;
use File::Copy;
use Storable qw(dclone);
use Scalar::Util qw( looks_like_number );
use Time::HiRes qw( gettimeofday tv_interval );

# Standalone VAL-only adaptation: all generated files stay in ROOT/PREP_VAL.
use Getopt::Long qw(GetOptions);
use FindBin;
use File::Path qw(rmtree);
use Text::ParseWords qw(shellwords);
my ($root,$common_file,$core_common_file,$edge_text,$edge_index,$dft_dir,$bin,$dry,$help);
$edge_index = 1;
GetOptions('root=s'=>\$root, 'common=s'=>\$common_file,
           'core-common=s'=>\$core_common_file, 'edge=s'=>\$edge_text,
           'edge-index=i'=>\$edge_index,
           'dft-dir=s'=>\$dft_dir, 'ocean-bin=s'=>\$bin,
           'dry-run'=>\$dry, 'help|h'=>\$help) or die "Use --help\n";
if ($help) {
  print "Usage: perl prep_val.pl [--root RUN] [--common JSON] [--dft-dir DIR]\n";
  print "                        [--core-common JSON] [--edge 'SITE N L']\n";
  print "                        [--edge-index N] [--ocean-bin DIR] [--dry-run]\n";
  print "Requires calc.mode=val and matching completed DFT. Creates fresh RUN/PREP_VAL\n";
  print "with both valence tmels and same-site core projectors for the selected edge.\n";
  exit 0;
}
die "Unexpected positional arguments\n" if @ARGV;
$root=abs_path($root // getcwd()) or die "Cannot resolve run root\n";
$common_file=abs_path($common_file // catfile($root,'Common','postDefaultsOceanDatafile'))
  or die "Cannot resolve Common JSON\n";
$dft_dir=abs_path($dft_dir // catdir($root,'DFT')) or die "Cannot resolve DFT directory\n";
$bin=abs_path($bin // $ENV{OCEAN_BIN} // $FindBin::Bin) or die "Cannot resolve OCEAN bin\n";
$ENV{OCEAN_BIN}=$bin;
my $target=catdir($root,'PREP_VAL');
die "PREP_VAL already exists; refusing to overwrite it\n" if -e $target || -l $target;
for my $file ('ocean_prep.x','specpnt.5','Pquadrature','sphpts') {
  die "Missing $bin/$file\n" unless -f catfile($bin,$file);
}
die "ocean_prep.x is not executable\n" unless -x catfile($bin,'ocean_prep.x');
my $t0 = [gettimeofday];

my $json = JSON::PP->new;
my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);

# Load run info from Common
my $dataFile = $common_file;
die "Failed to find $dataFile\n" unless( -e $dataFile );

my $commonOceanData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $commonOceanData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

# A valence Common file normally has no core edge.  Import one edge from the
# completed core Common metadata, or accept an explicit OCEAN edge triplet.
my $selected_edge;
my $selected_edge_source;
if (defined $edge_text) {
  my @edge = grep { length($_) } split /[\s,]+/, $edge_text;
  die "--edge must be 'SITE N L' with positive SITE/N and 0 <= L < N\n"
    unless @edge == 3 && !grep { !/^\d+$/ } @edge
      && $edge[0] >= 1 && $edge[1] >= 1 && $edge[2] >= 0 && $edge[2] < $edge[1];
  $selected_edge = join(' ', map { 0+$_ } @edge);
  $selected_edge_source = 'command line';
} else {
  die "--edge-index must be a positive 1-based index\n" unless $edge_index >= 1;
  $core_common_file //= catfile($root,'Common','postDefaultsOceanDatafile');
  $core_common_file = abs_path($core_common_file)
    or die "Cannot resolve --core-common; alternatively supply --edge 'SITE N L'\n";
  open my $core_in, '<', $core_common_file
    or die "Cannot read core Common JSON $core_common_file: $!\n";
  local $/ = undef;
  my $core_data = $json->decode(<$core_in>);
  close $core_in;
  my $edges = $core_data->{calc}{edges};
  die "$core_common_file has no calc.edges; supply --edge 'SITE N L'\n"
    unless ref($edges) eq 'ARRAY' && @$edges;
  die "--edge-index $edge_index exceeds the " . scalar(@$edges) . " available edges\n"
    if $edge_index > @$edges;
  $selected_edge = $edges->[$edge_index-1];
  my @edge = grep { length($_) } split /[\s,]+/, $selected_edge;
  die "Selected calc.edges entry is not a valid 'SITE N L' triplet\n"
    unless @edge == 3 && !grep { !/^\d+$/ } @edge
      && $edge[0] >= 1 && $edge[1] >= 1 && $edge[2] >= 0 && $edge[2] < $edge[1];
  $selected_edge = join(' ', map { 0+$_ } @edge);
  $selected_edge_source = $core_common_file . " entry " . $edge_index;
}
$commonOceanData->{calc}{edges} = [ $selected_edge ];

# Load info from OPF
$dataFile = catfile( $root, "OPF", "opf.json" );
my $opfData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $opfData = $json->decode(<$in>);
  close($in);
}
else
{
  $opfData = {};
}

# Load run info from DFT
$dataFile = catfile( $dft_dir, "dft.json" );
die "Failed to find $dataFile\n" unless( -e $dataFile );

my $dftData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $dftData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

my @selected_znl = split ' ', $selected_edge;
die "Selected edge site $selected_znl[0] exceeds the DFT atom count\n"
  unless $selected_znl[0] <= @{$dftData->{structure}{typat} // []};
my $selected_type = $dftData->{structure}{typat}[$selected_znl[0]-1];
my $selected_z = $dftData->{structure}{znucl}[$selected_type-1];
die "Missing OPF completion metadata for Z=$selected_z in $root/OPF/opf.json\n"
  unless defined($selected_z) && ref($opfData->{completed}{$selected_z}) eq 'HASH'
    && defined($opfData->{completed}{$selected_z}{psp_hash})
    && defined($opfData->{completed}{$selected_z}{input_hash});
die "Missing $root/OPF/zpawinfo required for core projectors and mels\n"
  unless -d catdir($root,'OPF','zpawinfo');

# Validate before creating any output; never relabel core metadata as valence.
die "Requires calc.mode=val in the selected Common JSON\n"
  unless ($commonOceanData->{calc}{mode} // '') eq 'val';
die "Selected DFT BSE stage is incomplete\n" unless $dftData->{bse}{complete};
for my $key ('kmesh','kshift') {
  die "Common/DFT mismatch or missing $key\n"
    unless ref($dftData->{bse}{$key}) eq 'ARRAY'
      && $json->encode($commonOceanData->{bse}{$key}) eq $json->encode($dftData->{bse}{$key});
}
for my $key ('photon_q','nonzero_q') {
  die "Common/DFT mismatch or missing $key; run matching valence DFT first\n"
    unless exists $dftData->{bse}{$key}
      && $json->encode($commonOceanData->{calc}{$key}) eq $json->encode($dftData->{bse}{$key});
}
die "Missing four-entry DFT brange\n"
  unless ref($dftData->{bse}{brange}) eq 'ARRAY' && @{$dftData->{bse}{brange}}==4;
my @br=@{$dftData->{bse}{brange}};
my $vfirst=$commonOceanData->{bse}{val_start}>1 ? $commonOceanData->{bse}{val_start} : 1;
die "DFT band window does not cover requested occupied/conduction states\n"
  unless $vfirst >= $br[0] && $vfirst <= $br[1] && $br[2]>0 && $br[3]>=$br[2];
# Probe the same paths that copyLinkQE/copyLinkABI will use.
my $program=$commonOceanData->{dft}{program};
die "Unsupported DFT program $program\n" unless $program eq 'qe' || $program eq 'abi';
my @km=@{$commonOceanData->{bse}{kmesh}};
my @ks=@{$commonOceanData->{bse}{kshift}};
my @q=@{$commonOceanData->{calc}{photon_q}};
die "Invalid k/q dimensions\n" unless @km==3 && @ks==3 && @q==3 && !grep {$_<=0} @km;
my @dirs;
if ($program eq 'qe' && $dftData->{bse}{split}) {
  push @dirs, sprintf('k%i_%i_%iq%.6f_%.6f_%.6f',@km,@ks);
  my @shift;
  for my $i (0..2) {
    my $v=$ks[$i]/$km[$i]-$q[$i];
    $v+=1 while $v<0; $v-=1 while $v>=1;
    push @shift,$v*$km[$i];
  }
  push @dirs, sprintf('k%i_%i_%iq%.6f_%.6f_%.6f',@km,@shift);
} else {
  die "ABINIT split finite-q is unsupported by this original PREP adapter\n"
    if $program eq 'abi' && $dftData->{bse}{split};
  my $prefix=($program eq 'qe' && $commonOceanData->{calc}{nonzero_q}) ? 'ks' : 'k';
  push @dirs, sprintf($prefix.'%i_%i_%iq%.6f_%.6f_%.6f',@km,@ks);
}
for my $dir (@dirs) {
  my $base=catdir($dft_dir,$dir);
  if ($program eq 'qe') {
    my $save=catdir($base,'Out','system.save');
    die "Missing QE wavefunction directory $save\n"
      unless -f catfile($save,'data-file-schema.xml') || -f catfile($save,'data-file.xml');
    die "Missing $base/QE_EIGS.txt\n"
      if -f catfile($save,'data-file-schema.xml') && !-f catfile($base,'QE_EIGS.txt');
  } else {
    die "Missing $base/NSCFx_WFK\n" unless -f catfile($base,'NSCFx_WFK');
  }
}
print "VAL-only PREP: $target\nCommon: $common_file\nDFT: $dft_dir\n";
print "Core edge: $selected_edge (from $selected_edge_source)\n";
print "No DFT, core PREP, or Common data will be modified.\n";
exit 0 if $dry;
mkdir($target) or die "Cannot create $target: $!\n";
chdir($target) or die "Cannot enter $target: $!\n";
copy($common_file,'valence_common.json') or die "Cannot save input snapshot: $!\n";
copy(catfile($dft_dir,'dft.json'),'valence_dft.json') or die "Cannot save DFT snapshot: $!\n";

# Grab previous run info if it exists
my $prepDataFile = "prep.json";
my $prepData = {};
if( -e $prepDataFile && open( my $in, "<", $prepDataFile ) )
{
  local $/ = undef;
  $prepData = $json->decode(<$in>);
  close($in);
}

# Build to-do list
my $newPrepData;


$newPrepData->{'bse'} = {};
$newPrepData->{'bse'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newPrepData->{'bse'}->{'complete'} = JSON::PP::false
    unless( exists $prepData->{'bse'} && exists $prepData->{'bse'}->{'complete'} 
            && $prepData->{'bse'}->{'complete'});

copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'dft'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, [ 'program' ] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, [ 'hash' ] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'scf'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, [ 'fermi', 'nelec' ] );

#copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'screen'}, $prepData->{'bse'},
#                $newPrepData->{'bse'}, [ 'core_offset' ] );

my $fake->{ 'complete' } = JSON::PP::false;
$newPrepData->{'computer'} = {};
copyAndCompare( $newPrepData->{'computer'}, $commonOceanData->{'computer'}, $prepData->{'computer'},
                $fake, [ 'para_prefix' ] );



##### currently needed files
# TODO: make a single input file that can be read all at once by ocean_prep
# avecs
# ntype
# natoms
# znucl
# typat
# taulist
# coord = xred
my @list = ( "avecs", "typat", "znucl", "xred", "elname", "bvecs" );
copyAndCompare( $newPrepData->{'bse'}, $dftData->{'structure'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, \@list );
copyAndCompare( $newPrepData->{'bse'}, $dftData->{'structure'}, $prepData->{'bse'},
                $fake, [ "epsilon"] );
#copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'structure'}, $prepData->{'bse'},
#                $newPrepData->{'bse'}, \@list );

# kmesh.ipt
# nspin
# k0
# qinunitsofbv
# brange
# split
# xmesh
@list = ( "kmesh", "kshift", "xmesh" );
copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, \@list );
copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'bse'}, $prepData->{'bse'},
                $fake, [ "val_start" ] );

copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'calc'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["nonzero_q", "photon_q"] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'general'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["nspin"] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["split"] );
copyAndCompare( $newPrepData->{'bse'}, $dftData->{'bse'}, $prepData->{'bse'},
                $fake, [ "brange", "max_occ_band", "min_unocc_band" ] );

if( $newPrepData->{'bse'}->{'val_start'} > 1 ) {
  printf "Start band %i\n", $newPrepData->{'bse'}->{'val_start'};
  $newPrepData->{'bse'}->{'brange'}[0] = $newPrepData->{'bse'}->{'val_start'};
} else {
  $newPrepData->{'bse'}->{'brange'}[0] = 1;
}

if( $newPrepData->{'bse'}->{'complete'} ) {
  for( my $i = 0; $i < 4; $i++ ) {
    if( $newPrepData->{'bse'}->{'brange'}[$i] != $prepData->{'bse'}->{'brange'}[$i] ) {
      $newPrepData->{'bse'}->{'complete'} = JSON::PP::false;
      printf "Brange doesn't match: %i  %i\n", $newPrepData->{'bse'}->{'brange'}[$i], $prepData->{'bse'}->{'brange'}[$i];
      last;
    }
  }
}




# Setting tmels and cks requires some logic
my $tmels =  JSON::PP::false;
my $cks =  JSON::PP::false;
my $calc = 'val';

# The core-valence pathway needs both ingredients in the same Bloch gauge:
# finite-q valence transition matrix elements and same-site core projectors.
$tmels = JSON::PP::true;
$cks = JSON::PP::true;


$newPrepData->{'bse'}->{'cks'}->{'enabled'} = JSON::PP::false;
$newPrepData->{'bse'}->{'tmels'} = JSON::PP::false;
# TODO, in the future make edges inclusive 
if( $cks )
{
  copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'calc'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["edges"] );

  if( ref $prepData->{'bse'}->{'cks'} eq ref {} ) { 
    unless( $prepData->{'bse'}->{'cks'}->{'enabled'} ) {
      $newPrepData->{'bse'}->{'complete'} = JSON::PP::false;
#      print "  prev enabled failed\n";
    }
  } else {
    $newPrepData->{'bse'}->{'complete'} = JSON::PP::false;
#    print "  ref check failed\n";
  }
  $newPrepData->{'bse'}->{'cks'}->{'enabled'} = JSON::PP::true;

  my $new = 0;
  unless( exists $prepData->{'bse'}->{'cks'} && 
          exists $prepData->{'bse'}->{'cks'}->{'previous'} ) {
    $prepData->{'bse'}->{'cks'}->{'previous'} = {};
    $new = 1;
  }
  $newPrepData->{'bse'}->{'cks'}->{'previous'} = {};

  my %Zset;
  foreach (@{$newPrepData->{'bse'}->{'edges'}}) {
    my @znl = split ' ', $_;
    $Zset{ $newPrepData->{'bse'}->{'znucl'}[$newPrepData->{'bse'}->{'typat'}[$znl[0]-1]-1] } = 1;
  }
  foreach my $z ( keys %Zset ) {
    if( exists $prepData->{'bse'}->{'cks'}->{'previous'}->{$z} ) {
      unless( exists $prepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'psp_hash'} && 
              $prepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'psp_hash'} 
              eq $opfData->{'completed'}->{$z}->{'psp_hash'} ) {
        print ": psp_hash differs\n" unless( $new);
        $new = 1;
      }
      unless( exists $prepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'input_hash'} &&
              $prepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'input_hash'}
              eq $opfData->{'completed'}->{$z}->{'input_hash'} ) {
        print ": input_hash differs\n" unless( $new );
        $new = 1;
      }
    } else {
#      print "  exists previous fail\n";
      $new = 1;
    }
    $newPrepData->{'bse'}->{'cks'}->{'previous'}->{$z} = {};
    $newPrepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'psp_hash'} 
        = $opfData->{'completed'}->{$z}->{'psp_hash'} ;
    $newPrepData->{'bse'}->{'cks'}->{'previous'}->{$z}->{'input_hash'} 
        = $opfData->{'completed'}->{$z}->{'input_hash'} ;
  }
  $newPrepData->{'bse'}->{'complete'} = JSON::PP::false if( $new );
#  print "   new\n" if( $new );
}
if( $tmels )
{
  $newPrepData->{'bse'}->{'complete'} = JSON::PP::false unless( $prepData->{'bse'}->{'tmels'} );
  $newPrepData->{'bse'}->{'tmels'} = JSON::PP::true;
}

$newPrepData->{'bse'}->{'complete'} = JSON::PP::false unless( exists $prepData->{'bse'}->{'cbm'} );
$newPrepData->{'bse'}->{'cbm'} = $prepData->{'bse'}->{'cbm'} if( $newPrepData->{'bse'}->{'complete'} );
$newPrepData->{'time'} = $prepData->{'time'} if( exists $prepData->{'time'} );
$newPrepData->{'dft_time'} = $dftData->{'time'};
$newPrepData->{'calc'}={'mode'=>'val'};
$newPrepData->{'prep_val'}={common_json=>$common_file,dft_dir=>$dft_dir,output_dir=>$target,
                            core_common_json=>($core_common_file // ''),
                            selected_edge=>$selected_edge,
                            selected_edge_source=>$selected_edge_source};

open OUT, ">", "prep.json" or die;
print OUT $json->encode($newPrepData);
close OUT;

unless( $newPrepData->{'bse'}->{'complete'} )
{
######## file nonsense
  mkdir "BSE" or die "Cannot create BSE: $!\n";
  chdir "BSE" or die "Cannot enter BSE: $!\n";

  writeOceanPrepInput( $newPrepData->{'bse'});

  my @extraFiles = ("specpnt.5", "Pquadrature", "sphpts" );
  foreach( @extraFiles )
  {
    copy( catfile( $ENV{'OCEAN_BIN'},$_), $_ ) or die $!;
  }

  if( $newPrepData->{'bse'}->{'program'} eq 'qe' )
  {
    copyLinkQE( $newPrepData->{'bse'});
  } elsif ( $newPrepData->{'bse'}->{'program'} eq 'abi' ) {
    copyLinkABI( $newPrepData->{'bse'});
  } else {
    die "Bad DFT: $newPrepData->{'program'}\n";
  }

#######
  if( -l "zpawinfo" )  # zpawinfo is an existing link
  { 
    unlink "zpawinfo" or die "Problem cleaning old 'zpawinfo' link\n$!";
  }
  elsif(  -d "zpawinfo" ) #or zpawinfo is existing directory
  { 
    rmtree( "zpawinfo" );
  }
  elsif( -e "zpawinfo" ) #or zpawinfo is some other file
  { 
    unlink "zpawinfo";
  }
  if( $newPrepData->{'bse'}->{'cks'}->{'enabled'} )
  {
    print catdir( updir(), updir(), "OPF", "zpawinfo" ) . "\n";
    symlink( catdir( updir(), updir(), "OPF", "zpawinfo" ), "zpawinfo" ) == 1 or die "Failed to link zpawinfo\n$!";
  }


  print "$newPrepData->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1\n";
  # Preserve the launcher arguments without shell interpolation.
  my @command=(shellwords($newPrepData->{'computer'}->{'para_prefix'} // ''),
               catfile($ENV{'OCEAN_BIN'},'ocean_prep.x'));
  my $child=fork();
  die "Cannot fork ocean_prep: $!\n" unless defined $child;
  if ($child==0) {
    open STDOUT,'>','ocean_prep.log' or die $!;
    open STDERR,'>&',STDOUT or die $!;
    exec {$command[0]} @command;
    die "Cannot execute ocean_prep: $!\n";
  }
  waitpid($child,0);
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "ocean_prep died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  ocean_prep exited with value $errorCode\n";
    }
    else {
      printf "ocean_prep exited successfully with value %d\n", $errorCode;
    }
  }

  for my $f ('val.u2.dat','con.u2.dat','ptmels.dat','tmels.info','enkfile','wvfcninfo','wvfvainfo') {
    die "ocean_prep returned success but $f is missing/empty\n" unless -s $f;
  }
  my @parcksv = grep { -s $_ } glob('parcksv.*');
  my @parcksc = grep { -s $_ } glob('parcksc.*');
  die "ocean_prep returned success but no parcksv.* projector was generated\n" unless @parcksv;
  die "ocean_prep returned success but no parcksc.* projector was generated\n" unless @parcksc;
  open IN, '<', 'cbm.out' or die "Failed to open cbm.out\n$!";
  my $eshift = <IN>;
  close IN;

  $newPrepData->{'bse'}->{'cbm'} = $eshift*1;

  chdir updir();
  
  $newPrepData->{'bse'}->{'complete'} = JSON::PP::true;
  $newPrepData->{'time'} = tv_interval( $t0 );

  open OUT, ">", "prep.json" or die;
  print OUT $json->encode($newPrepData);
  close OUT;
} 


exit 0;


sub copyLinkQE
{
  my ($hashRef) = @_;

  foreach my $d ( "Out", "Out_shift" ) {
    if( -l $d )  # Out is an existing link
    {
      unlink $d or die "Problem cleaning old '$d' link\n$!";
    }
    elsif(  -d $d ) #or Out is existing directory
    {
      rmtree( $d );
    }
    elsif( -e $d ) #or Out is some other file
    {
      unlink $d;
    }
  }

  # Three cases here. Either no finite-q, finite-q single directory, or finite-q two directories
  # Do the two-directory versions first
  if( $hashRef->{'split'} ) {
    my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $hashRef->{'kshift'}[0],
                   $hashRef->{'kshift'}[1], $hashRef->{'kshift'}[2];
    my $dirname = catdir( $dft_dir, $rundir, "Out" );
    symlink ($dirname, "Out_shift") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out_shift", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( $dft_dir, $rundir, "QE_EIGS.txt"), "QE_EIGS_shift.txt" );
    }

    my @kshift;
    for( my $i=0; $i<3; $i++ ) {
      my $kactual = $hashRef->{'kshift'}[$i]/$hashRef->{'kmesh'}[$i];
      $kactual -= $hashRef->{'photon_q'}[$i];
      while( $kactual < 0 ) { $kactual += 1; }
      while( $kactual >= 1 ) { $kactual -= 1; }
      $kactual *= $hashRef->{'kmesh'}[$i];
      $kshift[$i] = $kactual;
    } 
    my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $kshift[0],
                   $kshift[1], $kshift[2];
    my $dirname = catdir( $dft_dir, $rundir, "Out" );
    symlink ($dirname, "Out") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( $dft_dir, $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" );
    }
  } else {
    my $rundir;
    if( $hashRef->{'nonzero_q'} ) {
      $rundir = sprintf "ks%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                       $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $hashRef->{'kshift'}[0],
                       $hashRef->{'kshift'}[1], $hashRef->{'kshift'}[2];
    } else {
      $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                       $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $hashRef->{'kshift'}[0],
                       $hashRef->{'kshift'}[1], $hashRef->{'kshift'}[2];
    }
    my $dirname = catdir( $dft_dir, $rundir, "Out" );
    symlink ($dirname, "Out") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( $dft_dir, $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" );
    }
  }

  if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) )
  {
    print "Detected QE62-style DFT run\n";
    open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
    print TMP "qe62\n";
    close TMP;

#    copy( catfile( $dft_dir, $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" );
#    copy("../$rundir/enkfile", "enkfile_raw") or die "Failed to grab enkfile\n$!";
  }
  elsif( -e catfile( "Out", "system.save", "data-file.xml" ) )
  {
    print "Detected QE54-style DFT run\n";
    open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
    print TMP "qe54\n";
    close TMP;
  }
  else
  {
    die "Failed to detect QE style!\n";
  }

  open TMP, ">", "prefix" or die $!;
  print TMP "system\n";
  close TMP;

}

sub copyLinkABI
{
  my ($hashRef) = @_;


  my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $hashRef->{'kshift'}[0],
                   $hashRef->{'kshift'}[1], $hashRef->{'kshift'}[2];

  my $file = "NSCFx_WFK";
  my $targ = "RUN0001_WFK";
  if( -l $targ )  # Out is an existing link
  {
    unlink $targ or die "Problem cleaning old $targ link\n$!";
  }
  elsif( -e $targ ) #or Out is some other file
  {
    unlink $targ;
  }
  my $dirname = catdir( $dft_dir, $rundir, $file );
  symlink ($dirname, $targ) == 1 or die "Failed to link $targ\n$!";

  

  open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
  print TMP "abinit\n";
  close TMP;
}

sub writeOceanPrepInput
{
  my $hashRef = $_[0];
  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "% .16g  % .16g  % .16g\n", $hashRef->{'avecs'}[$i][0],
                                $hashRef->{'avecs'}[$i][1],
                                $hashRef->{'avecs'}[$i][2];

  }
  close OUT;

  open OUT, ">", "ntype" or die "$!";
  printf OUT "%i\n", scalar @{$hashRef->{'znucl'}};
  close OUT;

  open OUT, ">", "typat" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'typat'}}; $i++ )
  {
    printf OUT "%i\n", $hashRef->{'typat'}[$i];
  }
  close OUT;

  open OUT, ">", "natoms" or die "$!";
  printf OUT "%i\n", scalar @{$hashRef->{'typat'}};
  close OUT;

  open OUT, ">", "znucl" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'znucl'}}; $i++ )
  {
    printf OUT "%i\n", $hashRef->{'znucl'}[$i];
  }
  close OUT;

  open OUT, ">", "taulist" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'xred'}}; $i++ )
  {
    printf  OUT "% .16g  % .16g  % .16g\n", $hashRef->{'xred'}[$i][0],
                                $hashRef->{'xred'}[$i][1],
                                $hashRef->{'xred'}[$i][2];
  }
  close OUT;

  open OUT, ">", "coord" or die $!;
  print OUT "xred\n";
  close OUT;


  open OUT, ">", "kmesh.ipt" or die $!;
  printf  OUT "%i  %i  %i\n", $hashRef->{'kmesh'}[0],
                              $hashRef->{'kmesh'}[1],
                              $hashRef->{'kmesh'}[2];
  close OUT;

  open OUT, ">", "k0.ipt" or die $!;
  printf  OUT "% .16g  % .16g  % .16g\n", $hashRef->{'kshift'}[0],
                              $hashRef->{'kshift'}[1],
                              $hashRef->{'kshift'}[2];
  close OUT;

  open OUT, ">", "qinunitsofbvectors.ipt" or die $!;
  printf  OUT "% .16g  % .16g  % .16g\n", $hashRef->{'photon_q'}[0],
                              $hashRef->{'photon_q'}[1],
                              $hashRef->{'photon_q'}[2];
  close OUT;

  open OUT, ">", "xmesh.ipt" or die $!;
  printf  OUT "%i  %i  %i\n", $hashRef->{'xmesh'}[0],
                              $hashRef->{'xmesh'}[1],
                              $hashRef->{'xmesh'}[2];
  close OUT;

  #TODO fix when nbands is less than DFT supply (
  open OUT, ">", "brange.ipt" or die $!;
  printf OUT "%i  %i  %i  %i\n", $hashRef->{'brange'}[0], $hashRef->{'brange'}[1],
                                 $hashRef->{'brange'}[2], $hashRef->{'brange'}[3];
  close OUT;

  open OUT, ">", "nspin" or die $!;
  printf OUT "%i\n", $hashRef->{'nspin'};
  close OUT;

  open OUT, ">", "dft.split" or die $!;
  if( $hashRef->{'split'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;


  open OUT, ">", "prep.tmels" or die $!;
  if( $hashRef->{'tmels'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  open OUT, ">", "prep.cks" or die $!;
  if( $hashRef->{'cks'}->{'enabled'} ) {
    print OUT "1\nNA 1\n";
  } else {
    print OUT "0\n";
  } 
  close OUT;

  if( $hashRef->{'cks'}->{'enabled'} ) 
  {
    open IN, ">", "edges" or die $!;
#TODO: Fix edges, should be array of arrays?
#    for( my $i = 0; $i < scalar @{$hashRef->{'edges'}}; $i++ )
#    {
#      printf IN "%i  %i  %i\n", $hashRef->{'edges'}[$i][0], $hashRef->{'edges'}[$i][1], $hashRef->{'edges'}[$i][2];
#    }
    for( my $i = 0; $i < scalar @{$hashRef->{'edges'}}; $i++ )
    {
      print IN $hashRef->{'edges'}[$i] . "\n";
    }
    close IN;
  }

  open OUT, ">", "efermiinrydberg.ipt" or die $!;
  printf OUT "% .16g\n", $hashRef->{'fermi'}*2;
  close OUT;

#  open OUT, ">", "core_offset" or die $!;
#  print OUT $hashRef->{'core_offset'} . "\n";
#  close OUT;
   
}
sub copyAndCompare
{ 
  my $newRef = $_[0];
  my $commonRef = $_[1];
  my $oldRef = $_[2];
  my $complete = $_[3];
  my @tags = @{$_[4]};
  
  my $comp;# = sub { $_[0] == $_[1] }; 
  
  foreach my $t (@tags)
  { 
    if( ref( $commonRef->{ $t } ) eq '' )
    {
#      print "$commonRef->{ $t } ---\n"; 
      $newRef->{ $t } = $commonRef->{ $t };
    }
    else
    { 
      $newRef->{ $t } = dclone $commonRef->{ $t };
    }
    
    next unless( $complete->{'complete'} );
    unless( exists $oldRef->{ $t } )
    { 
      $complete->{'complete'} = JSON::PP::false;
      next;
    }
    
    recursiveCompare( $newRef->{$t}, $oldRef->{$t}, $complete);
    print "DIFF:   $t\n" unless( $complete->{'complete'} );
  }

}


sub recursiveCompare
{
  my $newRef = $_[0];
  my $oldRef = $_[1];
  my $complete = $_[2];

  return unless( $complete->{'complete'} );


  if( ref( $newRef ) eq 'ARRAY' )
  {
    if( scalar @{ $newRef } != scalar @{ $oldRef } )
    {
      $complete->{'complete'} = JSON::PP::false;
      return;
    }
    for( my $i = 0; $i < scalar @{ $newRef }; $i++ )
    {
      recursiveCompare( @{$newRef}[$i], @{$oldRef}[$i], $complete );
      return unless( $complete->{'complete'} );
    }
  }
  else
  {
#    print "#!  $newRef  $oldRef\n";
    if( looks_like_number( $newRef ) )
    {
      $complete->{'complete'} = JSON::PP::false unless( fcomp($newRef, $oldRef ) );
#      $complete->{'complete'} = JSON::PP::false unless( $newRef == $oldRef );
    }
    else
    {
      $complete->{'complete'} = JSON::PP::false unless( $newRef eq $oldRef );
    }
  }
}

sub fcomp {
  my ($a, $b ) = @_;
  return sprintf("%.12g", $a) eq sprintf("%.12g", $b); 
}
