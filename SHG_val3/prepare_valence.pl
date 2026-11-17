#!/usr/bin/env perl
# Prepare a separate finite-q valence first-resolvent calculation.
# Requires completed Common/DFT/PREP from calc.mode=val. Does not run DFT/BSE.
use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Basename qw(dirname);
use File::Spec;
use File::Path qw(make_path);
use FindBin;
use Getopt::Long qw(GetOptions);
use JSON::PP;
use POSIX qw(isfinite);

my ($root,$target,$axis,$con,$val,$erange,$eta,$mels_energy,$ocean_source,$ocean_bin);
my ($common_file,$dft_dir,$prep_dir);
my ($dry,$help) = (0,0);
GetOptions('root=s'=>\$root, 'target=s'=>\$target, 'axis=s'=>\$axis,
  'common=s'=>\$common_file, 'dft-dir=s'=>\$dft_dir, 'prep-dir=s'=>\$prep_dir,
  'con-bands=s'=>\$con, 'val-bands=s'=>\$val, 'energy-range=s'=>\$erange,
  'eta=f'=>\$eta, 'mels-energy=f'=>\$mels_energy,
  'ocean-source=s'=>\$ocean_source, 'ocean-bin=s'=>\$ocean_bin,
  'dry-run'=>\$dry, 'help|h'=>\$help) or usage(2);
usage(0) if $help;
die "--axis x, y, or z is required\n" unless defined($axis) && $axis =~ /^[xyz]$/;
$root = abs_path($root // getcwd()) or die "Cannot resolve --root\n";
$target //= File::Spec->catdir($root,'SHG_VAL');
$target = File::Spec->rel2abs($target);
my $parent = abs_path(dirname($target)) // '';
die "--target must be a new direct child of --root\n" unless $parent eq $root;
die "Target already exists; choose a fresh directory for this band-window test: $target\n"
  if -e $target || -l $target;
$common_file //= "$root/Common/postDefaultsOceanDatafile";
$dft_dir //= "$root/DFT";
$prep_dir //= "$root/PREP";
$common_file = abs_path($common_file) or die "Cannot resolve --common\n";
$dft_dir = abs_path($dft_dir) or die "Cannot resolve --dft-dir\n";
$prep_dir = abs_path($prep_dir) or die "Cannot resolve --prep-dir\n";
die "--common must identify postDefaultsOceanDatafile\n" unless -f $common_file;
die "--dft-dir must contain dft.json and rhoofr\n"
  unless -f "$dft_dir/dft.json" && -f "$dft_dir/rhoofr";
die "--prep-dir must contain prep.json and BSE\n"
  unless -f "$prep_dir/prep.json" && -d "$prep_dir/BSE";
my $common = read_json($common_file);
my $prep = read_json("$prep_dir/prep.json");
my $dft = read_json("$dft_dir/dft.json");
die "Requires completed DFT and PREP BSE stages\n" unless $prep->{bse}{complete} && $dft->{scf}{complete};
die "Requires a valence PREP calculation (calc.mode=val), not a core/RIXS run\n"
  unless ($common->{calc}{mode} // '') eq 'val';
die "Grid screening is not supported by this OCEAN staging path\n"
  if ($common->{screen}{mode} // '') eq 'grid';
die "PREP lacks same-site core projectors. Rerun prep_val.pl with --core-common or --edge.\n"
  unless $prep->{bse}{cks}{enabled}
    && ref($prep->{bse}{edges}) eq 'ARRAY' && @{$prep->{bse}{edges}};
my $b = $prep->{bse};
die "Missing PREP band range\n" unless ref($b->{brange}) eq 'ARRAY' && @{$b->{brange}} == 4;
my @stored = @{$b->{brange}};
my @active = @stored;
if (($common->{bse}{occupation}{type} // 'none') eq 'none') {
  $active[1] = $b->{max_occ_band} if defined $b->{max_occ_band};
  $active[2] = $b->{min_unocc_band} if defined $b->{min_unocc_band};
}
@active[0,1] = band_pair($val) if defined $val;
@active[2,3] = band_pair($con) if defined $con;
# OCEAN reconstructs the Fermi level from this window. Retain the frontier
# states (and the overlapping occupation region for a metal/spin run).
my $required_vlast = (($common->{bse}{occupation}{type} // 'none') eq 'none')
  ? ($b->{max_occ_band} // $stored[1]) : $stored[1];
my $required_cfirst = (($common->{bse}{occupation}{type} // 'none') eq 'none')
  ? ($b->{min_unocc_band} // $stored[2]) : $stored[2];
die "Keep the valence upper bound $required_vlast and conduction lower bound $required_cfirst for OCEAN Fermi finding\n"
  unless $active[1] == $required_vlast && $active[2] == $required_cfirst;
for my $j (0,2) {
  die "Requested bands are outside PREP coverage; regenerate DFT/PREP with more bands\n"
    unless $active[$j] >= $stored[$j] && $active[$j+1] <= $stored[$j+1]
       && $active[$j] <= $active[$j+1] && $active[$j] >= 1;
}

my @q = @{$b->{photon_q} // []};
my @common_q = @{$common->{calc}{photon_q} // []};
my $bv = $b->{bvecs};
die "Missing PREP q or reciprocal vectors\n" unless @q == 3 && ref($bv) eq 'ARRAY' && @$bv == 3;
die "Selected Common and PREP photon_q values do not match\n"
  unless @common_q == 3
      && !grep { abs($common_q[$_]-$q[$_]) > 1e-12 } 0..2;
# JSON stores the three lattice vectors as rows. OCEAN reads each row
# into one column of sys%bvec. Hence q_cart[mu] = sum_j bvecs[j][mu]*q[j].
my @qc = (0,0,0);
for my $mu (0..2) {
  for my $j (0..2) { $qc[$mu] += $bv->[$j][$mu]*$q[$j]; }
}
my $qlen = sqrt($qc[0]**2+$qc[1]**2+$qc[2]**2);
die "PREP must use a nonzero finite small photon_q\n" unless isfinite($qlen) && $qlen > 0;
my $idir = index('xyz',$axis);
my $qperp2 = 0.0;
for my $mu (0..2) {
  next if $mu == $idir;
  $qperp2 += $qc[$mu]**2;
}
my $qperp = sqrt($qperp2);
my $q_is_oblique = $qperp > 1e-8*$qlen;
if ($q_is_oblique) {
  warn sprintf(
    "WARNING: stored PREP q is not parallel to Cartesian %s.\n" .
    "         q_cart = (% .12e, % .12e, % .12e) bohr^-1\n" .
    "         |q_perp|/|q| = %.6e\n" .
    "         The generated source is sign(q_%s)*(qhat dot d), not pure d_%s.\n",
    $axis, $qc[0], $qc[1], $qc[2], $qperp/$qlen, $axis, $axis
  );
}
my $pb = "$prep_dir/BSE";
for my $f (qw(enkfile wvfcninfo wvfvainfo tmels.info ptmels.dat con.u2.dat val.u2.dat)) {
  die "Missing $pb/$f\n" unless -f "$pb/$f";
}
die "Missing $dft_dir/rhoofr\n" unless -f "$dft_dir/rhoofr";
open my $ti, '<', "$pb/tmels.info" or die $!;
my @tm = split ' ', do { local $/; <$ti> };
close $ti;
die "tmels.info needs nbv, first_conduction, last_conduction, nk, nspin\n" unless @tm == 5;
my $nk = $b->{kmesh}[0]*$b->{kmesh}[1]*$b->{kmesh}[2];
die "tmels.info does not cover the requested band/k window\n"
  unless $tm[0] >= $active[1] && $tm[1] <= $active[2] && $tm[2] >= $active[3] && $tm[3] == $nk;
my $bytes = 16*$tm[0]*($tm[2]-$tm[1]+1)*$nk*$b->{nspin};
die "ptmels.dat size does not match tmels.info and nspin\n" unless -s "$pb/ptmels.dat" == $bytes;

my $gm = $common->{bse}{val}{gmres};
my @energies;
if (defined $erange) {
  @energies = number_list($erange);
  die "--energy-range needs MIN MAX STEP\n" unless @energies == 3;
  die "Energy range requires STEP>0 and MAX>=MIN\n" unless $energies[2] > 0 && $energies[1] >= $energies[0];
}
$eta //= $common->{bse}{val}{broaden};
die "A finite positive --eta (eV) is required\n" unless defined($eta) && isfinite($eta) && $eta > 0;
for my $key (qw(nloop gprc ffff)) {
  die "Missing bse.val.gmres.$key\n" unless defined $gm->{$key};
}
die "GMRES requires nloop>=2, gprc>0, ffff>0\n"
  unless $gm->{nloop} >= 2 && $gm->{gprc} > 0 && $gm->{ffff} > 0;
my $bsein = sprintf("0 0 0\n0 0 0\ninv\n%d %.16g %.16g %.16g\n",
                   $gm->{nloop},$eta,$gm->{gprc},$gm->{ffff});
my $first_energy;
if (@energies || ($gm->{estyle} // '') eq 'range') {
  @energies = @{$gm->{erange} // []} unless @energies;
  die "Invalid bse.val.gmres.erange\n" unless @energies == 3 && $energies[2] > 0 && $energies[1] >= $energies[0];
  $first_energy = $energies[0];
  $bsein .= "loop\n";
  $bsein .= sprintf("%.16g %.16g %.16g\n", @energies);
} elsif (($gm->{estyle} // '') eq 'list') {
  my $list = $gm->{elist};
  die "Empty GMRES energy list\n" unless ref($list) eq 'ARRAY' && @$list;
  $first_energy = $list->[0];
  $bsein .= "list\n".scalar(@$list)."\n".join('',map { sprintf("%.16g\n",$_) } @$list);
} else { die "Specify --energy-range or configure a GMRES list/range\n"; }
$mels_energy //= $first_energy;
die "--mels-energy must be a finite positive energy in eV\n"
  unless defined($mels_energy) && isfinite($mels_energy); #&& $mels_energy > 0;

$ocean_source //= File::Spec->catdir($FindBin::Bin,'..');
$ocean_source = abs_path($ocean_source) or die "Invalid --ocean-source\n";
$ENV{SHG_OCEAN_SOURCE} = $ocean_source;
$ocean_bin //= $ENV{OCEAN_BIN};
$ENV{OCEAN_BIN} = abs_path($ocean_bin) // die "Invalid --ocean-bin\n" if defined $ocean_bin;
for my $f (qw(Pquadrature sphpts meljtv.x)) {
  my @candidates = ("$ocean_source/Common/$f", "$ocean_source/CNBSE/zcain/zexe/$f");
  push @candidates,"$ENV{OCEAN_BIN}/$f" if defined $ENV{OCEAN_BIN};
  die "Cannot find OCEAN support file $f; supply --ocean-source or --ocean-bin\n"
    unless grep { -f $_ } @candidates;
}
print "Target: $target\nVAL bands: $active[0]:$active[1]; CON bands: $active[2]:$active[3]\n";
print "Common: $common_file\nDFT: $dft_dir\nPREP: $prep_dir\n";
printf "q_cart (bohr^-1): %.12g %.12g %.12g; polarization: %s\n",@qc,$axis;
printf "Core-valence mels: photon1=x, photon2=y, photon3=z at %.12g eV\n",$mels_energy;
print "TDA solve; GMRES energies are physical valence excitation energies in eV.\n";
exit 0 if $dry;
make_path($target);
chdir $target or die $!;
# The numbering is fixed so prepare_rbxc_cv.py can select the middle
# polarization without another meljtv.x run.
my @cartesian = ([1,0,0],[0,1,0],[0,0,1]);
for my $i (0..2) {
  my $e = join(' ',@{$cartesian[$i]});
  write_text('photon'.($i+1),
    "dipole\ncartesian $e\nend\ncartesian 0 0 0\nend\n" .
    sprintf("%.16g\n",$mels_energy));
}
$ENV{SHG_VAL_COMMON_FILE} = $common_file;
$ENV{SHG_VAL_DFT_DIR} = $dft_dir;
$ENV{SHG_VAL_PREP_DIR} = $prep_dir;
$ENV{SHG_VAL_BUILD_CV_ASSETS} = 1;
system($^X, "$FindBin::Bin/stage_valence.pl") == 0 or die "Valence staging failed; inspect $target\n";
my @parcksv = grep { -s $_ } glob('parcksv.*');
my @parcksc = grep { -s $_ } glob('parcksc.*');
my @mels = grep { -s $_ } glob('mels.z*n*l*p*');
die "Staging did not produce parcksv.* and parcksc.* projectors\n"
  unless @parcksv && @parcksc;
die "Staging did not produce all three Cartesian mels files\n" unless @mels >= 3;
write_text('brange.ipt',"$active[0] $active[1]\n$active[2] $active[3]\n");
write_text('file_brange.ipt',"$stored[0] $stored[1]\n$stored[2] $stored[3]\n");
write_text('qinunitsofbvectors.ipt',join(' ',map { sprintf('%.16g',$_) } @q)."\n");
write_text('shg_component',"$axis $axis $axis\n");
write_text('bse.in',$bsein);
write_text('cnbse.write_rhs',".false.\n");
write_text('echamp.inp',".false.\n");
write_text('screen.mode',($common->{screen}{mode} // 'hll')."\n");
write_text('nbuse.ipt',($active[1]-$active[0]+1)."\n");
my $manifest = {format_version=>1, purpose=>'finite-q valence Xc first resolvent',
  root=>$root, axis=>$axis, q_reduced=>\@q, q_cart_bohr_inverse=>\@qc,
  q_perpendicular_fraction=>$qperp/$qlen,
  q_parallel_to_requested_axis=>$q_is_oblique ? JSON::PP::false : JSON::PP::true,
  common_file=>$common_file, dft_dir=>$dft_dir, prep_dir=>$prep_dir,
  core_edges=>$prep->{bse}{edges}, mels_energy_ev=>$mels_energy,
  photon_index_map=>{x=>1,y=>2,z=>3},
  generated_projectors=>[map { File::Spec->rel2abs($_) } @parcksv,@parcksc],
  generated_mels=>[map { File::Spec->rel2abs($_) } @mels],
  active_brange=>\@active, stored_brange=>\@stored, eta_ev=>$eta,
  energy_reference=>'valence_transition_energy_no_core_edge_shift',
  source_convention=>'d_code=-i*sign(q_axis)*conjg(M_vn(q))/abs(q)',
  tda=>JSON::PP::true, upstream_common=>$common, upstream_prep=>$prep};
write_text('shg_valence_stage.json',JSON::PP->new->canonical->pretty->encode($manifest));
print "Prepared. Run shg_xc_minus_valence.x from $target\n";

sub read_json {
  open my $f,'<',$_[0] or die "Cannot read $_[0]: $!\n";
  my $data = JSON::PP->new->decode(do { local $/; <$f> });
  close $f;
  return $data;
}
sub write_text {
  open my $f,'>',$_[0] or die "Cannot write $_[0]: $!\n";
  print {$f} $_[1];
  close $f or die $!;
}
sub number_list {
  my @n = split /[\s,:]+/, $_[0];
  for (@n) { die "Invalid numeric input\n" unless /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/ && isfinite(0+$_); }
  return @n;
}
sub band_pair {
  my @p = number_list($_[0]);
  die "Band range must be FIRST:LAST (positive integer DFT indices)\n"
    unless @p == 2 && $p[0] >= 1 && $p[1] >= $p[0] && int($p[0]) == $p[0] && int($p[1]) == $p[1];
  return @p;
}
sub usage {
  print <<'HELP';
Usage: perl prepare_valence.pl --root RUN --axis z [options]
  --target DIR             New direct child of RUN (default RUN/SHG_VAL)
  --common FILE           Valence postDefaultsOceanDatafile
                          (default RUN/Common/postDefaultsOceanDatafile)
  --dft-dir DIR           Matching completed valence DFT directory
                          (default RUN/DFT)
  --prep-dir DIR          Matching completed valence PREP directory
                          (default RUN/PREP)
  --con-bands FIRST:LAST    Active conduction DFT band indices, inclusive
  --val-bands FIRST:LAST    Active occupied DFT band indices, inclusive
  --energy-range 'MIN MAX STEP'  First-resolvent loop grid, eV
  --eta EV                 Positive broadening; default bse.val.broaden
  --mels-energy EV         Photon energy written to Cartesian photon1/2/3;
                           default first GMRES energy
  --ocean-source DIR       OCEAN source root for support tables
  --ocean-bin DIR          Installed OCEAN support files (or OCEAN_BIN)
  --dry-run                Validate without creating the target

Uses completed finite-q VAL PREP outputs containing both tmels and same-site
core projectors. It immediately stages parcksv/parcksc and generates Cartesian
mels files with p01=x, p02=y, and p03=z. It does not generate additional bands
or change the PREP q vector.
Use a new target for each convergence test, or edit only brange.ipt within an
existing stage. Do not edit file_brange.ipt, tmels.info, or binary headers.
Keep the conduction lower and valence upper bounds at their default values;
these frontier states are needed by the original OCEAN Fermi-level finder.
HELP
  exit $_[0];
}

