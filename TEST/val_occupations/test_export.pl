#!/usr/bin/env perl
# Exercise the real preparation routine without launching the CNBSE driver.
use strict;
use warnings;
use Test::More;
use FindBin;
use lib "$FindBin::Bin/../../DFT";
use File::Spec::Functions qw(catfile catdir updir);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Cwd qw(getcwd);
use Scalar::Util qw(looks_like_number);
use POSIX qw(isfinite);
use JSON::PP;

open my $source, '<', "$FindBin::Bin/../../CNBSE/cnbse.pl" or die $!;
my $code = do { local $/; <$source> };
close $source;
$code =~ /(sub prepareValOccupations\s*\{.*?)(?=\nsub writeValAuxFiles)/s or die 'Preparation routine missing';
eval $1;
die $@ if $@;

my $cwd = getcwd();
my $tmp = tempdir(CLEANUP => 1);
make_path("$tmp/CNBSE", "$tmp/PREP/BSE/Out/system.save", "$tmp/PREP/BSE/Out_shift/system.save");
chdir "$tmp/CNBSE" or die $!;
my $common = { dft => { program => 'qe' }, bse => { val => { occupation => {
    type => 'fermi', temperature => 500, mu_mode => 'conserve',
} } } };

for my $spin (1, 2) {
  for my $mesh_case ('zero_q', 'paired_q', 'split_q') {
    my $split = $mesh_case eq 'split_q';
    my $meshes = $mesh_case eq 'paired_q' ? 2 : 1;
    my $prep = { bse => { split => $split, nonzero_q => $mesh_case ne 'zero_q' } };
    my $run = { calc => { mode => 'val' }, bse => {
      val => { gw => { control => 'none' } }, nspin => $spin, kmesh => [1,1,2], brange => [1,1,2,3],
      complete => JSON::PP::true,
    } };
    my $xml = '../PREP/BSE/' . ($split ? 'Out_shift' : 'Out') . '/system.save/data-file-schema.xml';
    open my $fh, '>', $xml or die $!;
    print {$fh} "<input><starting_k_points>\n<k_point weight=\"99\">0 0 0</k_point>\n</starting_k_points></input>\n";
    my @expected;
    for my $k (0 .. 2*$meshes-1) {
      my @e = map { -1 + $_/10 + $k/100 } 0 .. 3;
      push @e, map { $_+0.01 } @e if $spin == 2;
      push @expected, [@e] if $k % $meshes == $meshes-1;
      print {$fh} "<ks_energies>\n<k_point weight=\"0.5\">0 0 0</k_point>\n";
      print {$fh} '<eigenvalues size="' . scalar(@e) . '">', ($k%2 ? "\n" : ' '), join(' ', @e);
      print {$fh} "</eigenvalues>\n</ks_energies>\n";
    }
    close $fh;
    prepareValOccupations($run, $common, $prep, {});
    open $fh, '<', 'val_occupation_energies.dat' or die $!;
    is(scalar <$fh>, "4 2 $spin\n", "$mesh_case spin=$spin header");
    for my $s (0 .. $spin-1) {
      for my $k (0 .. 1) {
        my @values = split ' ', scalar <$fh>;
        is(shift @values, 0.5, 'normalized integration weight');
        is_deeply([map { 0+$_ } @values], [@{$expected[$k]}[4*$s .. 4*$s+3]], 'correct mesh and spin endpoint');
      }
    }
    close $fh;
    ok(!$run->{bse}{complete}, 'new input invalidates old spectrum');
    my $old = decode_json(encode_json($run));
    $run->{bse}{complete} = JSON::PP::true;
    prepareValOccupations($run, $common, $prep, $old);
    ok($run->{bse}{complete}, 'unchanged input retains cache validity');
    $common->{bse}{val}{occupation}{temperature} = 300;
    prepareValOccupations($run, $common, $prep, $old);
    ok(!$run->{bse}{complete}, 'temperature change invalidates cache');
    if ($spin == 1 && $mesh_case eq 'zero_q') {
      unlink $xml or die $!;
      $run->{bse}{complete} = JSON::PP::true;
      $run->{bse}{val}{occupation} = $common->{bse}{val}{occupation};
      $common->{bse}{val}{occupation}{temperature} = 0;
      my $warning = '';
      {
        local $SIG{__WARN__} = sub { $warning .= shift };
        prepareValOccupations($run, $common, $prep, $old);
      }
      like($warning, qr/0 K.*Using legacy valence occupations/, 'zero temperature warns');
      is($run->{bse}{val}{occupation}{type}, 'none', 'auxiliary writer receives legacy mode');
      is($common->{bse}{val}{occupation}{type}, 'fermi', 'requested settings are preserved');
      ok(!$run->{bse}{complete}, 'zero temperature invalidates a previous finite-temperature run');
      ok(!-e 'val_occupation.ipt' && !-e 'val_occupation_energies.dat',
         'zero temperature needs no QE XML and removes stale thermal files');
      for my $invalid (-1, 'NaN', 'Inf') {
        $common->{bse}{val}{occupation}{temperature} = $invalid;
        eval { prepareValOccupations($run, $common, $prep, $old); };
        like($@, qr/must be finite and nonnegative/, "reject invalid temperature $invalid");
      }
    }
    $common->{bse}{val}{occupation}{type} = 'none';
    prepareValOccupations($run, $common, $prep, $old);
    ok(!-e 'val_occupation.ipt' && !-e 'val_occupation_energies.dat', 'legacy mode removes stale thermal files');
    $common->{bse}{val}{occupation} = {type=>'fermi', temperature=>500, mu_mode=>'conserve'};
    $run->{bse}{val}{gw}{control} = 'band';
    eval { prepareValOccupations($run, $common, $prep, {}); };
    like($@, qr/only gw.control = none or cstr/, 'unsupported GW control is not silently ignored');
  }
}
chdir $cwd or die $!;
done_testing();
