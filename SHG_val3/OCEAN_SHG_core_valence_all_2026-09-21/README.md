# OCEAN SHG core-valence standalone workflow

This directory consolidates the complete core-valence pathway into one source
tree. It contains the valence Common/PREP helpers, the first valence GMRES
solve, the core-valence vertex, the second core GMRES solve, and the final
contraction.

## Build

Build OCEAN first, then place this directory directly below the OCEAN source
root and run:

```bash
make -j8
```

If the directory is elsewhere:

```bash
make -j8 OCEAN_ROOT=/path/to/OCEAN
```

The default target creates:

- `shg_xc_minus_valence.x`
- `shg_rbxc_cv.x` (MPI implementation)
- `shg_y1_cv.x`
- `shg_t1_gv_contract.x`

The older serial vertex implementation is retained for validation and can be
built with `make serial-rbxc`.

## Workflow

1. `ocean_common_val.pl`: generate the separate valence `Common_VAL` tree.
2. Run the matching valence DFT stage.
3. `prep_val.pl`: generate `PREP_VAL` with both finite-q valence `tmels` and
   same-gauge core projectors for one selected core edge.
4. `prepare_valence.pl`: stage the first resolvent and immediately generate
   `parcksv/parcksc` plus Cartesian `mels` (`p01=x`, `p02=y`, `p03=z`).
5. Run `shg_xc_minus_valence.x`.
6. `prepare_rbxc_cv.py`: stage the core-valence vertex inputs.
7. Run `shg_rbxc_cv.x` with MPI.
8. `prepare_y1_cv.pl`: stage the core second-resolvent calculation.
9. Run `shg_y1_cv.x`.
10. Run `shg_t1_gv_contract.x`.

See the component-specific documents under `docs/` for complete command-line
options, data layouts, and validation rules.

The core edge can be imported from a completed core Common file:

```bash
perl prep_val.pl --root /run --common /run/Common_VAL/postDefaultsOceanDatafile \
  --dft-dir /run/DFT_VAL --core-common /run/Common/postDefaultsOceanDatafile \
  --edge-index 1 --ocean-bin /path/to/OCEAN/bin
```

Alternatively use `--edge 'SITE N L'`. The selected site must have completed
OPF data under `ROOT/OPF`.

## Oblique finite-q behavior

`prepare_valence.pl` and `shg_valence_source.f90` no longer terminate when the
stored Cartesian q vector is not parallel to the requested axis. They print a
warning and continue. The resulting source is

```text
d_code = sign(q_axis) * (qhat dot d),
```

not a pure `d_axis` source. The preparation manifest records both the
perpendicular-q fraction and whether q is parallel to the requested axis. A
zero q remains fatal.

## Included source files

The shared `shg_gmres.f90` serves both GMRES drivers. The two copies in the
original component packages differed only by a trailing blank line. The MPI
`shg_rbxc_cv_mpi.f90` is the default vertex program; `shg_rbxc_cv.f90` is the
serial reference implementation.
