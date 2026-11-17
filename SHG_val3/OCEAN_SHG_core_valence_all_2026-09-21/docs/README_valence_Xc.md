# Standalone valence first resolvent for OCEAN SHG

This add-on computes the valence-sector first resolvent

\[
X_V^{(c)}(\omega)=[(\omega+i\eta)I-H_V]^{-1}d_V^{(c)}.
\]

It is a separate executable, analogous to the core `shg_xc_minus` program. It uses OCEAN's valence Hamiltonian and vector machinery and the signed GMRES module from the earlier core add-on. It does not require changes to OCEAN2 sources. It does require a compiled OCEAN tree: “standalone” here means a separate driver, not independence from OCEAN libraries.

This is the first step of the core–valence pathway. It does not yet apply R_CV, solve the second core resolvent, or contract the final susceptibility.

## Compatibility and prerequisites

The interfaces were checked against the supplied `OCEAN_SHG_core_standalone_original_2026-08-19` source tree. The package targets its finite-q VAL path and modern split `val.u2.dat` / `con.u2.dat` storage, with `bloch_selector=3`, `enk_selector=0`, and `tmel_selector=1`.

Start from completed Common, DFT, and PREP stages of a **valence** calculation (`calc.mode=val`). A core or RIXS directory alone is insufficient. The three stages may use their conventional names or be selected explicitly with `--common`, `--dft-dir`, and `--prep-dir`, including the side-by-side `Common_VAL`, `DFT_VAL`, and `PREP_VAL` layout. PREP must contain a nonzero small momentum transfer. A q vector parallel to the requested Cartesian polarization gives a pure Cartesian source. An oblique q is accepted with warnings and produces the directional mixture `sign(q_axis)*(qhat dot d)`, not a pure `d_axis` source. Changing `--axis` does not generate a new shifted mesh. Both signs of q are supported. Converge the result as q tends to zero.

The driver uses the Tamm–Dancoff approximation (one forward electron–hole sector). The preparation script disables the backward and semitda paths. It retains the original valence interaction settings. This staging route does not support grid screening or the legacy combined u2 file.

## Build

Copy the `SHG_VAL` directory alongside `OCEAN2` in the OCEAN source root:

```bash
cp -r SHG_VAL /path/to/OCEAN/
# Build OCEAN and its support modules using your usual procedure first.
make -C /path/to/OCEAN/SHG_VAL
```

The Makefile reuses `Makefile.arch`, OCEAN object files, and module include directories. In particular, it includes `../Modules` for FFT_wrapper.mod. Use the same compiler and MPI implementation as the OCEAN build.

## Prepare and run

Example: valence bands end at global DFT band 8; conduction bands start at 9; PREP stores at least through band 48. Adapt these indices to your material.

```bash
perl /path/to/OCEAN/SHG_VAL/prepare_valence.pl \
  --root /path/to/completed_valence_run \
  --target /path/to/completed_valence_run/SHG_VAL_N40 \
  --ocean-source /path/to/OCEAN \
  --axis z \
  --con-bands 9:48 \
  --energy-range '2480.5 2498.0 0.5' \
  --eta 0.3

cd /path/to/completed_valence_run/SHG_VAL_N40
mpirun -np 4 /path/to/OCEAN/SHG_VAL/shg_xc_minus_valence.x > run.log 2>&1
```

For the side-by-side workflow produced by `ocean_common_val.pl` and
`prep_val.pl`, preserve the ordinary core directories and select the valence
inputs explicitly:

```bash
perl /path/to/OCEAN/SHG_VAL/prepare_valence.pl \
  --root /path/to/ocean_run \
  --target /path/to/ocean_run/SHG_VAL_N40 \
  --common /path/to/ocean_run/Common_VAL/postDefaultsOceanDatafile \
  --dft-dir /path/to/ocean_run/DFT_VAL \
  --prep-dir /path/to/ocean_run/PREP_VAL \
  --ocean-source /path/to/OCEAN \
  --axis z \
  --con-bands 9:48 \
  --energy-range '2480.5 2498.0 0.5' \
  --eta 0.3
```

`--common` names the processed valence `postDefaultsOceanDatafile` itself;
the other two options name directories containing `dft.json`/`rhoofr` and
`prep.json`/`BSE`, respectively.  When omitted, they remain backward
compatible with `ROOT/Common`, `ROOT/DFT`, and `ROOT/PREP`.  The selected
absolute paths are recorded in `shg_valence_stage.json`, and the staging
helper uses those exact locations for wavefunctions, energies, transition
matrices, and density.

The energy range and broadening are in eV. These are **physical fundamental photon energies**, not edge-relative energies and not the second-resolvent 2omega energies. For a final core excitation near 4966 eV, a corresponding fundamental is near 2483 eV. No core-edge shift is added to the valence Hamiltonian.

`--target` must be a fresh direct child of `--root`. The script stages a separate calculation and does not launch DFT or the solver. `--dry-run` validates inputs without creating it. `--help` lists options. Without an energy override it uses the Common valence GMRES energy list/range. The inherited GMRES iteration, preconditioning, and convergence settings are written into `bse.in`; inspect these before running. A range includes points min+n*increment up to max, never an extra point above max.

The selected PREP must contain `cks.enabled=true` and a core edge. During
staging the script copies `parcksv/parcksc` and invokes OCEAN `meljtv.x` for
three generated electric-dipole photon files. Their fixed convention is
`p01=x`, `p02=y`, and `p03=z`; their default photon energy is the first GMRES
energy and can be overridden with `--mels-energy EV`. The created filenames
and convention are recorded in `shg_valence_stage.json`.

## Conduction-band convergence

The active range already has a standard input file: **brange.ipt**.

```text
1 8
9 48
```

Line 1 contains the first and last valence bands; line 2 contains the first and last conduction bands. All bounds are inclusive global DFT band indices. Here nval=8 and ncon=40.

For another test, stage a fresh directory with `--con-bands 9:68`, provided PREP already contains those bands. The script writes the selected `brange.ipt`. It keeps **file_brange.ipt** at the full stored PREP range so file reads remain correctly aligned. Do not change file_brange.ipt to the active range or truncate the wavefunction/overlap files.

The conduction lower bound and valence upper bound are kept at the frontier bounds inferred from PREP, including occupation overlap when applicable. OCEAN needs these states when finding the Fermi level. You can vary the conduction upper bound and, with `--val-bands FIRST:LAST`, the valence lower bound. This is deliberately not an arbitrary interior-band selector.

The script cannot manufacture additional DFT states: prepare a sufficiently large maximum window upstream first. OCEAN also applies its occupation/energy admission mask. In this source version the default upper clipping energy depends on the last retained conduction band across k and spin. Thus enlarging the window can change both vector dimension and admitted transitions. Record the manifest and staging provenance for every test, and converge the final observable as well as band count, q, k mesh, and GMRES residual.

## Source convention and normalization

For an occupied band v and empty band n, the original finite-q overlap is

\[
M_{vn}(\mathbf k,\mathbf q)=\langle u_{v\mathbf k}|u_{n,\mathbf k+\mathbf q}\rangle.
\]

OCEAN's selector-1 reader constructs S=conjg(M)/|q|. For q parallel to axis c, this add-on converts it to the length source

\[
d^{(c)}_{nv\mathbf k}=-i\,M_{vn}^*/q_c
=-i\,\mathrm{sign}(q_c)S_{nv\mathbf k}.
\]

For an oblique q, the implemented phase instead produces

\[
d_{\mathrm{code}}=\mathrm{sign}(q_c)
\sum_\mu \widehat q_\mu d^{(\mu)}.
\]

The preparation script and executable both print warnings in this case. The
staging manifest records `q_perpendicular_fraction` and
`q_parallel_to_requested_axis` so downstream calculations cannot silently
mistake the directional mixture for a pure Cartesian source.

The phase matters in coherent SHG contractions. The original reader applies its masks and normalizes S. The add-on recovers the masked source norm from OCEAN's saved prefactor:

\[
\|d\|=\sqrt{\mathrm{kpref}\,N_k\,\Omega^2/(4\pi)}.
\]

Here Omega is the cell volume in bohr cubed and Nk is the total k-point count. GMRES operates with a unit source; the dumped d and X are restored to physical normalization. Energies are converted internally to Hartree. d has units bohr and X has units bohr/Hartree; the electron charge is not included.

At finite q the Hamiltonian and source use an electron at k+q and a hole at k. This approaches the optical source as q tends to zero. The later R_CV contraction must match the electron band window, momentum mesh, and Bloch gauge; matching array sizes alone is insufficient.

## Outputs and array layout

For z polarization:

- `shg_cache_d_valence_z.bin`: physical masked length source, written once.
- `shg_cache_Gm_valence_w_z_000001.bin`, etc.: physical X at every requested energy.
- `shg_xc_minus_valence_z.info`: dimensions, energy list, q, band ranges, normalization, residuals, iteration counts, and file names.
- `shg_valence_stage.json`: preparation provenance and selected range.

Binary files are native-endian stream-unformatted complex(DP) values (complex128 in the targeted build), without Fortran sequential record markers. The compact Fortran shape is

```fortran
complex(real64) :: X(ncon,nval,nk,nbeta,1)
```

ncon varies fastest, followed by nval, k, valence spin channel beta, and the single TDA sector. No internal OCEAN padding is dumped. A file contains ncon*nval*nk*nbeta complex values; use the manifest rather than guessing nbeta. These are valence transition arrays, not core arrays of shape (ncon,nk,nalpha), and not sequential-unformatted echamp files.

For example, with NumPy:

```python
X = np.fromfile(filename, dtype=np.complex128).reshape(
    (ncon, nval, nk, nbeta, 1), order='F')
```

The driver stops on nonfinite or unconverged residuals before writing that solution. The inherited GMRES tolerance is a squared residual-norm threshold; the manifest also reports residual norms and their physical-source scaling. Earlier completed energy files may remain if a later solve fails.

## Code map

- `shg_xc_minus_valence.f90`: initialize MPI/system, enforce VAL/TDA, prepare d, set up the negative-H resolvent, solve all frequencies, restore normalization, and write compact vectors/manifest.
- `shg_valence_source.f90`: validate q and overlap coverage, call OCEAN_load_data, recover normalization, and rotate the source to the length-gauge phase using public vector operations.
- `shg_gmres.f90`: reused signed arbitrary-RHS GMRES module; the driver selects negative H and disables core spectral source factors.
- `prepare_valence.pl`: validate explicitly selected or conventional upstream data, choose band/energy windows, and create an isolated run directory.
- `stage_valence.pl` and `OCEANcompare.pl`: adapted OCEAN input staging and support routines. The staging helper exits without running a solver.

## Verification and remaining build test

The included synthetic staging tests check two conduction windows, unchanged stored-band metadata and overlaps, exact frequency lists, and rejection of invalid inputs. They passed in the development environment. Perl syntax checks and Fortran parsing also passed. **No Fortran compiler or compiled OCEAN libraries were available here, so this package has not been compiled or run against the actual Hamiltonian.**

```bash
python3 tests/test_staging.py
# After building, run the optional independent-particle end-to-end check:
python3 tests/test_staging.py --exe /path/to/OCEAN/SHG_VAL/shg_xc_minus_valence.x
# MPI variant:
python3 tests/test_staging.py --exe /path/to/OCEAN/SHG_VAL/shg_xc_minus_valence.x --mpi 2
```

The optional test compares dumped components against d=-i*M*/q and X=d/(omega+i*eta-(epsilon_n-epsilon_v)) in a synthetic no-interaction problem. It is supplied for your build environment and was not run here. After it passes, compare a physical valence calculation against OCEAN's response with consistent source phase, normalization, Hamiltonian flags, and broadening before using it for the core–valence SHG contraction.

The included OCEAN-derived routines retain the upstream license in `LICENSE-OCEAN`.
