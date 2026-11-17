# Standalone OCEAN core-SHG workflow

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

The OCEAN-linked programs build against an unmodified OCEAN source tree. All additions are confined to
`SHG_core_realspace2/`; no source file or Makefile under `OCEAN2/`, `CNBSE/`, or another
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
- `shg_y1_core.f90` reads either RbXc manifest and its selected physical vector,
  solves the second core resolvent at a requested target energy, and writes a
  physical Y1 vector and companion manifest.
- `shg_t1_contract.f90` computes the raw T1 contraction between the physical
  core dipole and Y1 vectors, after checking their manifests.
- `prepare_shg.pl` and `stage_core_xas.pl` stage exactly one existing XAS edge
  and one Cartesian photon into `RUN_ROOT/SHG`. The helper is based on the
  early/staging portion of original `CNBSE/cnbse.pl`; it runs original
  `meljtv.x` but deliberately does not run `ocean.x`.

The five executables cover the first resolvent, either representation of the
core position operator, the second resolvent, and the raw T1 contraction.
The raw T1 term is not the final second-order susceptibility.

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
make -C SHG_core_realspace2
make -C SHG_core_realspace2 check
```

The executables are in `SHG_core_realspace2/`: `shg_xc_minus.x`,
`shg_rbxc_core.x`, `shg_rbxc_core_realspace.x`, `shg_y1_core.x`, and
`shg_t1_contract.x`. The OCEAN-linked targets omit `OCEAN.o`, which contains
the `ocean.x` program entry point. The serial T1 target has no OCEAN link
dependency.

## Stage one core source

From the root of the completed XAS calculation:

```bash
/path/to/OCEAN/SHG_core_realspace2/prepare_shg.pl --component xyz --edge-index 1
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
mpirun -np 4 /path/to/OCEAN/SHG_core_realspace2/shg_xc_minus.x \
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
mpirun -np 4 /path/to/OCEAN/SHG_core_realspace2/shg_rbxc_core.x \
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

## Solve the second core resolvent and contract T1

After either RbXc program finishes, choose one row of its manifest and use
that row's vector filename. For the realspace result at source index 1:

```bash
mpirun -np 4 /path/to/OCEAN/SHG_core_realspace2/shg_y1_core.x \
  shg_rbxc_core_realspace_yz.info \
  shg_cache_Rrs_y_Gm_core_w_z_000001.bin TARGET_ENERGY_EV \
  shg_cache_Y1_core_yz_000001.bin
```

For an SHG T1 point, set `TARGET_ENERGY_EV` to twice the source energy in
the selected RbXc row. The program writes a `.info` file beside the Y1 vector.
The link-based RbXc manifest and its `shg_cache_R_y_...` vector are also
accepted.

Contract Y1 with a physical dipole vector for the left polarization `a`:

```bash
/path/to/OCEAN/SHG_core_realspace2/shg_t1_contract.x \
  shg_xc_minus_core_x.info shg_cache_d_core_x.bin \
  shg_cache_Y1_core_yz_000001.bin.info \
  shg_cache_Y1_core_yz_000001.bin shg_T1_core_xyz_000001.dat
```

The example uses `a=x`, `b=y`, `c=z`. The dipole manifest and Y1 manifest
must describe the same core edge, site, k mesh, and vector dimensions. The
contraction is `sum_i conjg(d_a(i))*Y_bc(i)` with no susceptibility prefactor.
Use the same staged XAS Hamiltonian for the RbXc and Y1 stages.
