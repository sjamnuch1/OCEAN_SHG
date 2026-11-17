# Finite-occupation implementation review

Reviewed on 2026-11-17, starting from implementation commit `048f0dd` on
`valence-finite-occupations`, based on the user's `shg` commit `6740dc5`.

The review found a spectral correctness bug in the dielectric restart and
fixed it. The checks below pass after the corrections. This is validation
with synthetic fixtures, not a completed BiVO4 production benchmark.

## Findings and corrections

| Finding | Correction | Files |
| --- | --- | --- |
| A dielectric restart reloads bare matrix elements without applying `sqrt(f_v-f_c)`. A fixed-Hamiltonian restart test changed the calculated spectrum. | Apply the same occupation weights before normalization on every reload. The new regression reproduces the failure before the fix and passes afterward. | `OCEAN2/OCEAN_load_data.f90`, `TEST/val_occupations/test_spectra.py` |
| The CNBSE exporter writes `none` for GW controls other than `cstr`; accepting another control could silently omit the requested correction. Some manual GW paths also have legacy assumptions outside this implementation's scope. | Finite-occupation mode explicitly accepts only `none` and `cstr`, in both preparation and the executable. Broadcast rejection before proceeding with MPI collectives. Legacy occupation mode is unchanged. | `CNBSE/cnbse.pl`, `OCEAN2/OCEAN_val_energy.f90` |
| The standalone count solver could report success for a nonfinite fixed chemical potential or invalid weights. Normal input parsing already caught the former, but the numerical routine did not enforce its own contract. | Validate dimensions, spin count, finite inputs, normalized positive weights, and finite results. | `OCEAN2/OCEAN_val_occupations.f90`, `TEST/val_occupations/test_occupations.f90` |
| Nonfinite cutoff or endpoint energies could bypass ordinary comparisons. The last-band convergence check assumes ordered eigenvalues. | Reject nonfinite energies/cutoffs and unordered counting spectra; check coverage of both optical band limits. | `OCEAN2/OCEAN_val_energy.f90`, `OCEAN2/OCEAN_val_occupations.f90` |

No existing variable was renamed or given a new type. The floating-point
`nelec` and separate pair-weight vector remain separate from the original
variables. `bse.val.occupation.type = none` remains the default.

## Validation

Build: GNU Fortran 14, Open MPI 5, BLAS/LAPACK and ScaLAPACK 2.2.2. The changed
Fortran sources were compiled with bounds checking. The legacy FFT routine
uses assumed-size conventions incompatible with those checks and was built
without bounds checking. This was a local Linux build, not a Perlmutter build.

| Check | Result |
| --- | --- |
| Full `ocean.x` build; Perl syntax; Python syntax; `git diff --check` | Pass |
| Analytic independent-particle spectra, fixed mu and conserved fractional count | Pass with one and two MPI ranks |
| Exchange response compared with explicitly scaled source and kernel rows/columns | Pass with one and two ranks |
| Dielectric restart compared with the unchanged-Hamiltonian reference | Pass with one and two ranks |
| Cold-limit equivalence to legacy spectra; finite runs at 300 K and 500 K | Pass with one and two ranks |
| Explicit two-spin degeneracy versus one doubly-degenerate channel | Pass with one and two ranks |
| Scissor correction changes optical poles while retaining raw-DFT occupations | Pass with one and two ranks |
| GMRES and exact-diagonalization exchange spectra versus explicit scaling | Pass with one and two ranks; exact comparison uses a tolerance suitable for its single-precision matrix |
| Invalid temperature/mu, inconsistent spectra, insufficient upper bands, zero source, NaN cutoff, backward mode and unsupported GW controls | Rejected with nonzero exit status, with one and two ranks |
| QE export for zero-q, paired-q and split-q, with one/two spins; cache invalidation and legacy cleanup | 72 Perl assertions pass |
| Low-temperature underflow, fractional counts, unequal spin spectra, weighted counting, fixed-mu residual and pair-weight cancellation | Numerical unit test passes |

Reproduction commands are in `DOC/valence_occupations.md`. The QE export test
executes the actual preparation subroutine against generated QE-style XML;
it does not substitute a second implementation of that subroutine.

## Remaining validation

No real BiVO4 DFT/PREP dataset was available in this checkout. Before using
the implementation for production results on Perlmutter:

1. Build with the intended compiler/MPI stack and run the included tests.
2. Reuse one completed BiVO4 snapshot. Compare legacy and finite-occupation
   spectra at 300 K and 500 K, checking logged electron counts and the upper-band
   occupation check. Repeat with the production MPI rank count.
3. Check convergence with NSCF bands, optical band ranges, k mesh and solver
   settings. Small changes in an insulating snapshot should not be interpreted
   without these convergence checks.

The exchange kernel has a numerical regression. Placement of the weights
around the ladder and ALDA kernels was inspected, but their full numerical
paths have not been validated with a material dataset. Core spectra and SHG
were not rerun. Finite occupations do not add electron-phonon self-energies,
temperature-dependent screening, new intraband transitions or a Drude response.
