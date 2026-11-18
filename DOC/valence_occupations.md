# Finite electronic occupations for valence spectra

Enable this mode in an OCEAN input for a forward valence UV–vis calculation.

```text
calc.mode { val }
bse.val.occupation.type { fermi }
bse.val.occupation.temperature { 500 }
bse.val.occupation.mu_mode { conserve }
```

The temperature is in kelvin. The default `bse.val.occupation.type { none }`
uses the original valence calculation, including its original electron-count
handling. Core occupations retain their separate `bse.occupation` inputs.
With `type { fermi }` and a temperature of exactly 0 K, CNBSE warns and uses
the legacy valence path. Negative or nonfinite temperatures are input errors
in Fermi mode.

The default chemical-potential mode, `conserve`, solves for the electron count
already supplied by DFT, including any fractional part. To supply a chemical
potential directly, set both `bse.val.occupation.mu_mode { fixed }` and
`bse.val.occupation.mu { <value in eV> }`. Use the energy zero of the raw QE
eigenvalues. In fixed mode the log reports the resulting count and its difference
from the DFT count.

CNBSE reads the original QE NSCF XML through the existing `PREP/BSE/Out` or
`Out_shift` link. Keep the XML when reusing completed DFT and PREP runs. For
finite photon momentum the chemical potential is determined on the unshifted
integration mesh. The valence and conduction endpoint energies are read from
`enkfile`. OCEAN checks the conduction bands against the counting spectrum
before using it. The highest included QE band must have occupation at most
`1e-8` at all counting points. If it exceeds this value, increase the NSCF band
count. The CNBSE log records the chemical potential, target and actual electron
counts, upper-band occupation, and active pair count and weight range.

For each forward pair, OCEAN computes `F = f_v - f_c` from the raw DFT
energies. Only positive-weight pairs in the prepared optical band space pass the
binary mask. The optical source is multiplied by `sqrt(F)` before
normalization. The interaction is multiplied on both sides by the diagonal
operator `S = sqrt(F)`, so the Hamiltonian is `D + S K S`. The one-electron
energy difference `D` uses the real scissor/stretch correction specified by
`bse.val.gw.control { cstr }`, while occupations retain the raw DFT energies.
This finite-occupation workflow accepts `gw.control = none` or `cstr` only;
other legacy/manual GW formats are rejected. The existing CNBSE exporter only
writes these two controls, so accepting other values would silently omit the
requested correction.

Use this mode with collinear QE band eigenstates and the uniform OCEAN BSE mesh.
It supports forward Haydock, GMRES and exact diagonalization. Backward-sector, backf,
Wannier-basis, spinor, and imaginary quasiparticle energies are unsupported.
The legacy CNBSE `doping` file is incompatible with this mode. The count solver
takes its target from DFT. The optical space comes from PREP, and the
calculation does not add a Drude term or additional intraband states.
For spin-polarized inputs, one common chemical potential conserves the total
electron count; separate spin populations or magnetization are not constrained.
Self-consistent dielectric restarts reapply the same occupation weights to
each freshly loaded optical source.

To isolate the electronic occupation contribution for BiVO4, run the same
snapshot with `type = none` and then with `type = fermi` at the desired
electronic temperature. Screening and wavefunctions are reused. Structural
temperature effects can still be studied by averaging spectra across the MD
snapshot ensembles. An intrinsic wide-gap snapshot can show very little
occupation-induced change over 300–500 K.

Run the spectral regressions with a compiled OCEAN:

```bash
python3 TEST/val_occupations/test_spectra.py OCEAN2/ocean.x
python3 TEST/val_occupations/test_spectra.py OCEAN2/ocean.x mpirun -n 2
```

These check independent-particle oscillator strengths, exchange-kernel scaling,
dielectric restarts, the cold limit, 300/500 K inputs, spin degeneracy, scissors,
GMRES, exact diagonalization, and input rejection. The fixtures are synthetic;
they do not establish convergence of a BiVO4 spectrum.

The export test needs only Perl and its core modules:

```bash
perl TEST/val_occupations/test_export.pl
```

Build the numerical unit test with the same Fortran compiler used for OCEAN
(replace `gfortran` as appropriate):

```bash
gfortran -I OCEAN2 TEST/val_occupations/test_occupations.f90 \
  OCEAN2/OCEAN_val_occupations.o -o /tmp/test_val_occupations
/tmp/test_val_occupations
```

This checks intrinsic carrier balance when occupation tails underflow,
fractional counts, unequal spin spectra, nonuniform integration weights,
fixed-mu count residuals, cancellation in pair weights, and invalid inputs.
See `DOC/valence_occupations_review.md` for the review findings and validation limits.
