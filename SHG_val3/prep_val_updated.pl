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
use Getopt::Long qw(GetOptions);

# Core-valence PREP options. With no options, preserve the existing directory
# convention: run from PREP_VAL, read ../Common_VAL and import the first edge
# from ../Common/postDefaultsOceanDatafile.
my ($coreCommonFile, $edgeText, $edgeIndex, $help);
$edgeIndex = 1;
GetOptions( 'core-common=s' => \$coreCommonFile,
            'edge=s'        => \$edgeText,
            'edge-index=i'  => \$edgeIndex,
            'help|h'        => \$help ) or die "Use --help\n";
if( $help ) {
  print "Usage: perl prep_val.pl [--core-common FILE] [--edge-index N]\n";
  print "                        [--edge 'SITE N L']\n";
  print "Generates valence tmels and same-gauge core projectors together.\n";
  exit 0;
}
die "Unexpected positional arguments\n" if @ARGV;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/prep\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
my $t0 = [gettimeofday];

my $json = JSON::PP->new;
my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);

# Load run info from Common
my $dataFile = catfile( updir(), "Common_VAL", "postDefaultsOceanDatafile");
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

# calc.mode=val normally does not retain a core edge. Import one edge from the
# completed core Common metadata, or accept an explicit OCEAN SITE/N/L triplet.
my $selectedEdge;
my $selectedEdgeSource;
if( defined $edgeText ) {
  my @edge = grep { length($_) } split /[\s,]+/, $edgeText;
  die "--edge must be 'SITE N L' with positive SITE/N and 0 <= L < N\n"
    unless @edge == 3 && !grep { !/^\d+$/ } @edge
      && $edge[0] >= 1 && $edge[1] >= 1
      && $edge[2] >= 0 && $edge[2] < $edge[1];
  $selectedEdge = join( ' ', map { 0 + $_ } @edge );
  $selectedEdgeSource = 'command line';
}
else {
  die "--edge-index must be a positive 1-based index\n" unless $edgeIndex >= 1;
  $coreCommonFile //= catfile( updir(), "Common", "postDefaultsOceanDatafile" );
  $coreCommonFile = abs_path( $coreCommonFile )
    or die "Cannot resolve --core-common; alternatively supply --edge 'SITE N L'\n";
  open my $coreIn, '<', $coreCommonFile
    or die "Cannot read core Common JSON $coreCommonFile: $!\n";
  local $/ = undef;
  my $coreData = $json->decode(<$coreIn>);
  close $coreIn;
  my $edges = $coreData->{'calc'}->{'edges'};
  die "$coreCommonFile has no calc.edges; supply --edge 'SITE N L'\n"
    unless ref($edges) eq 'ARRAY' && @{$edges};
  die "--edge-index $edgeIndex exceeds the " . scalar(@{$edges}) . " available edges\n"
    if $edgeIndex > scalar(@{$edges});
  my @edge = grep { length($_) } split /[\s,]+/, $edges->[$edgeIndex-1];
  die "Selected calc.edges entry is not a valid 'SITE N L' triplet\n"
    unless @edge == 3 && !grep { !/^\d+$/ } @edge
      && $edge[0] >= 1 && $edge[1] >= 1
      && $edge[2] >= 0 && $edge[2] < $edge[1];
  $selectedEdge = join( ' ', map { 0 + $_ } @edge );
  $selectedEdgeSource = "$coreCommonFile entry $edgeIndex";
}
$commonOceanData->{'calc'}->{'edges'} = [ $selectedEdge ];
print "Core edge for projector generation: $selectedEdge ($selectedEdgeSource)\n";

# Load info from OPF
$dataFile = catfile( updir(), "OPF", "opf.json" );
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
$dataFile = catfile( updir(), "DFT", "dft.json" );
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

my @selectedZNL = split ' ', $selectedEdge;
die "Selected edge site $selectedZNL[0] exceeds the DFT atom count\n"
  unless $selectedZNL[0] <= scalar(@{$dftData->{'structure'}->{'typat'} // []});
my $selectedType = $dftData->{'structure'}->{'typat'}[$selectedZNL[0]-1];
my $selectedZ = $dftData->{'structure'}->{'znucl'}[$selectedType-1];
die "Missing OPF completion metadata for Z=$selectedZ in ../OPF/opf.json\n"
  unless defined($selectedZ)
    && ref($opfData->{'completed'}->{$selectedZ}) eq 'HASH'
    && defined($opfData->{'completed'}->{$selectedZ}->{'psp_hash'})
    && defined($opfData->{'completed'}->{$selectedZ}->{'input_hash'});
die "Missing ../OPF/zpawinfo required for core-projector generation\n"
  unless -d catdir( updir(), "OPF", "zpawinfo" );

# Grab previous run info if it exists
my $prepDataFile = "prep.json";
my $prepData;
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
      #last;
    }
  }
}

printf "Check if we're doing Wannier\n";
#printf "%s\n",$prepData->{'bse'}->{'wann'};
#printf "%s\n",$newPrepData->{'bse'}->{'wann'};
#die;
#SJ
#Check if we're doing Wannier. If we're doing Wannier then check if the nwann is correct.
if( $prepData->{'bse'}->{'wann'} && $prepData->{'bse'}->{'complete'} ) 
  {
   printf "Detected Wannier run. Checking the band range\n";
   my $wannierFile = catfile( updir(), "DFT", "NWANN");
   open my $fh, '<', $wannierFile or die "Cannot open file: $!";
   chomp(my $NWANN = <$fh>);
   close $fh;
   printf "%i %i %i\n",$prepData->{'bse'}->{'brange'}[1],$prepData->{'bse'}->{'brange'}[0],$NWANN;
   #printf "%s\n",$newPrepData->{'bse'}->{'brange'}[1]-$newPrepData->{'bse'}->{'brange'}[0]+1 == $NWANN;
   if ( $prepData->{'bse'}->{'brange'}[1]-$prepData->{'bse'}->{'brange'}[0]+1 == $NWANN)
     {
      printf "The ran was done with Wannier size of %i\n",$prepData->{'bse'}->{'brange'}[1]-$prepData->{'bse'}->{'brange'}[0]+1;
      printf "Expecting %i\n",$NWANN;
      printf "PREP stage should be completed.";
      $newPrepData->{'bse'}->{'complete'} = JSON::PP::true;
     }
  }



# Setting tmels and cks requires some logic
my $tmels =  JSON::PP::false;
my $cks =  JSON::PP::false;

die "This core-valence PREP script requires calc.mode=val in Common_VAL\n"
  unless( ($commonOceanData->{'calc'}->{'mode'} // '') eq 'val' );

# R_CV requires both quantities to be generated from the same Bloch states:
# finite-q valence transition matrix elements and localized core projectors.
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

#SJ check if we have done wannier to use as valence band instead
##Currently we are in the PREP directory
#Look for the WANN file in DFT.
#If it exists, open it and assign the json file with bse->wannier
my $wannierFile = catfile( updir(), "DFT", "WANN");
if ( -e $wannierFile)
{
	print "Found WANN file in DFT. If the content is true then we will be using Wannier functions for valence bands\n";
open WANN_file, $wannierFile;
<WANN_file> =~m/(\w+)/ or die;
my $WANN_mode= $1;
close WANN_file;
if( $WANN_mode =~ m/true/i )
{
$newPrepData->{'bse'}->{'wann'} = JSON::PP::true ;
}
elsif( $WANN_mode =~ m/false/i )
{
$newPrepData->{'bse'}->{'wann'} = JSON::PP::false ;
}
else
{
	        print "WARNING!!! Wannier file in DFT phase is incorrectly specified. Setting the flag to false\n";
		print "Not using Wannier functions for valence bands and proceed with standard OCEAN way\n";
$newPrepData->{'bse'}->{'wann'} = JSON::PP::false ;		  
}
}
else
{
        print "Proceed with standard OCEAN without wannier\n";
	$newPrepData->{'bse'}->{'wann'} = JSON::PP::false ;
}

#SJ
#Update brange for number of valence to the size of wannier function
#Should have option to get nwann from input
if ($newPrepData->{'bse'}->{'wann'})
{
my $wannierFile = catfile( updir(), "DFT", "NWANN");
open my $fh, '<', $wannierFile or die "Cannot open file: $!";
#<WANN_file> =~m/(\w+)/ or die;
chomp(my $NWANN = <$fh>);
close $fh;
#my $NWANN= $1;
#close WANN_file;
print "$NWANN $newPrepData->{'bse'}->{'brange'}[1] $newPrepData->{'bse'}->{'brange'}[0]\n";
if ($newPrepData->{'bse'}->{'brange'}[1]-$newPrepData->{'bse'}->{'brange'}[0]+1 - $NWANN > 0)
{
	printf "The wannier function built is smaller than the full valence bands\n";
	printf "Adjusting brange to go from %i to %i \n",$newPrepData->{'bse'}->{'brange'}[0],$NWANN;
	#$newPrepData->{'bse'}->{'brange'}[0] = $newPrepData->{'bse'}->{'brange'}[1] - $NWANN + 1;
        $newPrepData->{'bse'}->{'brange'}[1] = $NWANN;
	printf "brange[0:1] %i %i\n",$newPrepData->{'bse'}->{'brange'}[0],$newPrepData->{'bse'}->{'brange'}[1];
} #Comapre nval to nwann
} #Check bse->wann loop

#print "STOP\n";
#die;
#SJ


open OUT, ">", "prep.json" or die;
print OUT $json->encode($newPrepData);
close OUT;


unless( $newPrepData->{'bse'}->{'complete'} )
{
######## file nonsense
  mkdir "BSE";
  chdir "BSE";

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


  #print "$newPrepData->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1\n";
  print "srun -n 1 $ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1\n";
  #system("$newPrepData->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1" );
  system("srun -n 1 $ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1" );
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

  for my $required ( 'val.u2.dat', 'con.u2.dat', 'ptmels.dat', 'tmels.info',
                     'enkfile', 'wvfcninfo', 'wvfvainfo' ) {
    die "ocean_prep returned success but $required is missing or empty\n"
      unless -s $required;
  }
  my @parcksv = grep { -s $_ } glob('parcksv.*');
  my @parcksc = grep { -s $_ } glob('parcksc.*');
  die "ocean_prep returned success but no parcksv.* projector was generated\n"
    unless @parcksv;
  die "ocean_prep returned success but no parcksc.* projector was generated\n"
    unless @parcksc;

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
  # SJ this is the default way.
  if( $hashRef->{'split'} ) { 
    my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $hashRef->{'kshift'}[0],
                   $hashRef->{'kshift'}[1], $hashRef->{'kshift'}[2];
    my $dirname = catdir( updir(), updir(), "DFT", $rundir, "Out" );
    symlink ($dirname, "Out_shift") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out_shift", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( updir(), updir(), "DFT", $rundir, "QE_EIGS.txt"), "QE_EIGS_shift.txt" ); #Here is conduction bands.
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
    if ( $hashRef->{'wann'} ) {
    my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f_WANN", $hashRef->{'kmesh'}[0],
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $kshift[0],
	           $kshift[1], $kshift[2];	
    my $dirname = catdir(updir(),updir(),"DFT",$rundir, "Out");
    symlink ($dirname,"Out") == 1 or die "Failed to link wannier Out\n$!";
    if( -e catfile( "Out","system.save","wann")) {
      copy(catfile(updir(),updir(),"DFT",$rundir,"QE_EIGS_WANN.txt"),"QE_EIGS.txt");
    } }
    else { 
    my $rundir = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],  
                   $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2], $kshift[0],
                   $kshift[1], $kshift[2];
    my $dirname = catdir( updir(), updir(), "DFT", $rundir, "Out" );
    symlink ($dirname, "Out") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( updir(), updir(), "DFT", $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" ); #Here is valence bands.
    } }
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
    my $dirname = catdir( updir(), updir(), "DFT", $rundir, "Out" );
    symlink ($dirname, "Out") == 1 or die "Failed to link Out\n$!";
    if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) ) {
      copy( catfile( updir(), updir(), "DFT", $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" );
    }
  }

  if( -e catfile( "Out", "system.save", "data-file-schema.xml" ) )
  {
    print "Detected QE62-style DFT run\n";
    open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
    print TMP "qe62\n";
    close TMP;

#    copy( catfile( updir(), updir(), "DFT", $rundir, "QE_EIGS.txt"), "QE_EIGS.txt" );
#    copy("../$rundir/enkfile", "enkfile_raw") or die "Failed to grab enkfile\n$!";
  }
  elsif( -e catfile( "Out","system.save","wann"))
  {
    print "Detected Wannier function as valence band built in QE62-style DFT run\n";
    open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
    print TMP "qe62\n";
    close TMP;
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
  my $dirname = catdir( updir(), updir(), "DFT", $rundir, $file );
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

