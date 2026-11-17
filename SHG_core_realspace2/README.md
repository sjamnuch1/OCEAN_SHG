# Standalone OCEAN core-SHG vector tests

2026-10-22 update: both `shg_rbxc_core.x` and
`shg_rbxc_core_realspace.x` accept nonzero upstream `photon_q` while using
OCEAN's base-mesh `con.u2.dat` for the optical core operator. See
[`NONZERO_PHOTON_Q_2026-10-22.md`](NONZERO_PHOTON_Q_2026-10-22.md) for the
changes, assumptions, build commands, and validation limits. This updated
package contains sources; rebuild the executables with your OCEAN compiler.

The first executable computes the core/XAS response vectors

```text
X_core,c^-(omega) = [(omega + i*eta) I - H_BSE^core]^-1 d_core^c.
```

The second executable applies the core-sector exciton--exciton position
operator to every saved first-layer vector,

```text
q_b,c(omega) = R_core^b X_core,c^-(omega).
```

An additional executable, `shg_rbxc_core_realspace.x`, evaluates this same
finite-mesh operator through a Born--von Karman real-space FFT and periodic
sine coordinate. See [`REALSPACE_RX.md`](REALSPACE_RX.md) for its derivation,
array shapes, build/run instructions, and link-versus-real-space comparison.

Both programs build against an unmodified OCEAN source tree. All additions are confined to
`SHG/`; no source file or Makefile under `OCEAN2/`, `CNBSE/`, or another
upstream directory is changed.

## What is implemented

- `shg_core_source.f90` calls original `OCEAN_load_data`. OCEAN's existing
  core/XAS path reads the selected edge and photon from `runlist`, constructs
  `d_core^c` from its CKS coefficients and `mels.*` matrix elements, applies
  the core transition mask, normalizes the source, and initializes the core
  BSE Hamiltonian.
- `shg_gmres.f90` is a private copy/adaptation of OCEAN's recycled GMRES. The
  signed-resolvent, arbitrary-RHS, solver-only, and returned-vector changes
  are compiled only into this standalone executable.
- `shg_xc_minus.f90` restores the physical dipole norm after GMRES, writes the
  core source and each converged response vector, and records diagnostics.
- `shg_core_vertex.f90` reads the staged cell-periodic conduction states from
  `con.u2.dat`, normalizes them, forms the neighboring-k conduction links,
  constructs the localized core-hole phase, and applies one reduced-direction
  contribution to the Cartesian pair-position operator.
- `shg_rbxc_core.f90` reads `shg_xc_minus_core_c.info` and all listed
  `shg_cache_Gm_core_w_c_*.bin` files, accumulates every reduced-direction
  contribution to `R_core^b`, reapplies OCEAN's core admission mask, and writes
  the physical `R_b X_core,c^-(omega)` vectors.
- `shg_core_vertex_realspace.f90` precomputes the full-Bloch phases, periodic
  finite-difference coordinate, and OCEAN FFT workspace, then applies the
  real-space operator without collapsing any core channel.
- `shg_rbxc_core_realspace.f90` mirrors the link driver and writes distinct
  `Rrs` caches plus a real-space manifest.
- `prepare_shg.pl` and `stage_core_xas.pl` stage exactly one existing XAS edge
  and one Cartesian photon into `RUN_ROOT/SHG`. The helper is based on the
  early/staging portion of original `CNBSE/cnbse.pl`; it runs original
  `meljtv.x` but deliberately does not run `ocean.x`.

The two executables now cover the first-layer exciton-vector calculation and
the core-to-core exciton coupling source. They do not yet apply the second
resolvent or evaluate the final second-order susceptibility.

## Required completed calculation

The calculation root must already contain completed `Common/`, `DFT/`, `OPF/`,
`PREP/`, and `SCREEN/` stages, with `calc.mode` set to `xas` and at least one
entry in `calc.edges`. In particular, the staging step consumes:

- `Common/postDefaultsOceanDatafile`;
- `DFT/dft.json` and the completed DFT/BSE state;
- `OPF/opf.json` plus the core radial/projector data;
- `PREP/prep.json` and `PREP/BSE/` wavefunctions, energies, and CKS files;
- `SCREEN/screen.json` and the completed screened core potentials.

The original OCEAN source must also be built first so that its objects and
`meljtv.x` are available.

## Build

Configure and build the unmodified OCEAN tree normally, then build the local
standalone program:

```bash
cp Makefile.arch.example Makefile.arch
# Edit Makefile.arch for the local compiler, MPI, BLAS, and FFTW installation.
make
make -C SHG
```

The executables are `SHG/shg_xc_minus.x`, `SHG/shg_rbxc_core.x`, and
`SHG/shg_rbxc_core_realspace.x`. Their
Makefile links the original OCEAN objects but omits `OCEAN.o`, which contains
the `ocean.x` program entry point.

## Stage one core source

From the root of the completed XAS calculation:

```bash
/path/to/OCEAN/SHG/prepare_shg.pl --component xyz --edge-index 1
```

This creates `./SHG`. For component `xyz`, the rightmost label selects the
source, so the staged source is `d_core^z`. The first two labels are retained
for the later SHG tensor contraction but do not affect this first-layer solve.
The component may contain repeated axes, so `xxx` is valid.

`--edge-index` is the one-based index in `calc.edges`. If original OCEAN
executables/support files are installed outside the source tree, pass
`--ocean-bin /path/to/bin`. Run `prepare_shg.pl --help` for all options.

The generated runlist has the ordinary original-OCEAN core form:

```text
1
Z  n  l  element  corelevel  site-index  photon-index  XAS
```

No `VAL` line is used.

## Run

```bash
cd SHG
mpirun -np 4 /path/to/OCEAN/SHG/shg_xc_minus.x \
  > shg_xc_minus.log 2>&1
```

The frequency grid, GMRES controls, and energy broadening `eta` come from the
ordinary `bse.core.gmres` and `bse.core.broaden` input, staged into `bse.in`.

For a rightmost axis `c`, the outputs are:

```text
shg_cache_d_core_c.bin
shg_cache_Gm_core_w_c_000001.bin
shg_cache_Gm_core_w_c_000002.bin
...
shg_xc_minus_core_c.info
```

The binary cache files are raw `complex(DP)` streams in original
`OCEAN_psi_vtor` order. The `.info` file records the selected edge/site,
energies, `eta`, convergence data, source norm, and vector filenames.

## Apply the core pair-position vertex

After `shg_xc_minus.x` finishes, remain in the same staged `SHG` run directory
and run

```bash
mpirun -np 4 /path/to/OCEAN/SHG/shg_rbxc_core.x \
  shg_xc_minus_core_z.info > shg_rbxc_core.log 2>&1
```

The command-line manifest is optional; when omitted it defaults to
`shg_xc_minus_core_z.info`. The program reads the component stored in that
manifest. For component `xyz`, it applies `R_core^y` to the `c=z` first-layer
vectors.

For `xyz`, the outputs are

```text
shg_cache_R_y_Gm_core_w_z_000001.bin
shg_cache_R_y_Gm_core_w_z_000002.bin
...
shg_rbxc_core_yz.info
```

The implementation is the covariant central difference

```text
R_core^b X(k) = i sum_j w_j^b [
    L_h,j^+ U_c,j^+(k) X(k+dk_j)
  - L_h,j^- U_c,j^-(k) X(k-dk_j) ],
```

where `U_c` is formed directly from the staged cell-periodic conduction
wavefunctions, `L_h,j^+ = exp(+i 2*pi*tau_j/N_j)`, and
`L_h,j^- = conjg(L_h,j^+)`. All three reduced directions are included when
required by a nonorthogonal cell. The Cartesian weights are stored in bohr in
the output manifest.

This stage does not call GMRES or apply the core BSE Hamiltonian. It uses the
original `OCEAN_system`, `OCEAN_psi`, and `OCEAN_energies` modules for the
staged run, vector layout, and final core transition mask.
