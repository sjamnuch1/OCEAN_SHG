# Batch core–valence second resolvent — 2026-10-02

This release replaces the single-target `shg_y1_cv.x` command line with a batch
interface. Compile it against the same OCEAN installation as the earlier package:

```bash
make -j8 OCEAN_ROOT=/path/to/OCEAN
```

No compiled executable is supplied. The shared `shg_gmres.f90` changed, so rebuild
its object and both solver executables through the Makefile rather than reusing
old `.o` or `.mod` files.

## Prepare the full grid

Run in the directory containing the completed `shg_rbxc_cv.info` and RHS binaries:

```bash
perl /path/to/prepare_y1_cv.pl \
  --ocean-run-root /path/to/OCEAN_RUN \
  --valence-con /path/to/PREP_VAL/BSE/con.u2.dat \
  --component zzz --core-eta 0.2
```

The script still stages a compatible core Hamiltonian, conduction range, compact
projectors when needed, basis certificate, and freshly generated `mels` using
`meljtv.x`. It now writes `bse2.in` with a full `list` of
`Omega(i) = 2 * omega(i)` from the CV manifest. It takes depth, preconditioner,
and convergence from the staged core GMRES configuration. Broadening is the core
setting, or `--core-eta EV`; valence broadening is retained only as provenance.

Alternatively request a full second-solver loop:

```bash
perl /path/to/prepare_y1_cv.pl \
  --ocean-run-root /path/to/OCEAN_RUN \
  --valence-con /path/to/PREP_VAL/BSE/con.u2.dat \
  --component zzz --energy-range=-461:-459:0.25 --core-eta 0.2
```

Or use `--bse2-input /path/to/custom_bse2.in` for a full existing `list` or `loop`
input. `--energy-range` and `--bse2-input` are mutually exclusive. The supplied
file's first two lines are replaced with the actual staged core header, while
its GMRES settings and grid are retained; `--core-eta` overrides its broadening
if requested. Custom target energies may differ from twice the manifest source
energy for diagnostic calculations. The solver warns and records
`target_is_two_omega F`; it still pairs rows by index, without interpolation.

Every grid must contain exactly the full manifest `energy_count`. Staging sorts
manifest rows by original integer index and requires contiguous indices from 1.
`cv_vectors.list` and `shg_gv_stage.json` record the actual requested targets.

`bse.in` remains a staged core input because unchanged `OCEAN_multiplet` reads its
first two header lines. The second GMRES solver reads its depth, broadening,
preconditioner, tolerance, and energies exclusively from `bse2.in`. The first
valence solver continues to read `bse.in` in its valence stage. As before, in-place
Y staging replaces valence Hamiltonian inputs with core inputs; do not rerun the
first solver in that stage without preparing the valence inputs again.

## Run all entries or begin later

All entries, with defaults:

```bash
mpirun -np 4 /path/to/shg_y1_cv.x --component zzz
```

Entries 5 through the last entry (9 in the supplied example):

```bash
mpirun -np 4 /path/to/shg_y1_cv.x --component zzz --start-index 5
```

The same selection by filename:

```bash
mpirun -np 4 /path/to/shg_y1_cv.x --component zzz \
  --start-file shg_cache_RbXc_cv_zz_000005.bin
```

`--start-index` and `--start-file` are mutually exclusive. `--end-index N` optionally
limits the inclusive last index. These options do not change the grid or renumber
sources or outputs. Entries before the selected range are not read by the solver.
When using an external manifest path, relative RHS paths resolve beside it.

| Original index | Source omega (eV) | Second target Omega (eV) |
|---|---:|---:|
| 1 | -230.500 | -461.00 |
| 5 | -230.000 | -460.00 |
| 9 | -229.500 | -459.00 |

Keep **all nine targets** in `bse2.in` when starting at 5. Entry 5 uses target 5,
not target 1. Values use the solver energy reference; the program does not add
an edge offset automatically.

Additional options:

```text
--manifest FILE       default shg_rbxc_cv.info
--basis FILE          default cv_core_basis_match.info
--output-prefix PATH  default shg_cache_Y1_cv_BC
--dipole-output FILE  default shg_cache_d_core_A.bin
--help
```

`--component ABC` is required; the B and C labels must match the manifest, and A
must match the final dipole polarization prepared by the staging script. A batch
uses one fixed edge, site, basis, and component. Prepare a separate compatible
stage for a different edge or final dipole polarization.

## Outputs and contraction

For `--component zzz`, default outputs are:

```text
shg_cache_d_core_z.bin
shg_cache_d_core_z.bin.info
shg_cache_Y1_cv_zz_000001.bin
shg_cache_Y1_cv_zz_000001.bin.info
...
shg_cache_Y1_cv_zz_000009.bin
shg_cache_Y1_cv_zz_000009.bin.info
shg_cache_Y1_cv_zz.batch.info
```

Each Y retains its original source index and the existing `Y1_GV_COMPLETE`
per-vector manifest format. The physical final dipole is reconstructed/saved
once per batch and referenced by every Y manifest. The contraction executable
is unchanged; for entry 5:

```bash
/path/to/shg_t1_gv_contract.x \
  shg_cache_d_core_z.bin.info shg_cache_d_core_z.bin \
  shg_cache_Y1_cv_zz_000005.bin.info shg_cache_Y1_cv_zz_000005.bin \
  shg_T1_gv_zzz_000005.dat
```

Selected existing Y outputs and the dipole are replaced on rerun. Y files outside
the requested range are retained. The batch manifest describes only the latest
requested range and gains `completion Y1_GV_BATCH_COMPLETE` after its last entry.
After an interruption, completed per-vector products remain available and the
batch completion marker is absent. Restart explicitly using the desired index.
Use a different output prefix/dipole filename to preserve a different setup.

## Computation and memory

For each i:

\[
[(\Omega_i+i\eta_C)I-H_C]Y_{bc,i}=R_{g\leftarrow v}^{b}X_V^c(\omega_i).
\]

The core Hamiltonian and its wavefunction data are initialized once per MPI run.
Each selected RHS is read, masked, normalized, solved at its own target, and its
physical norm restored. The solution is written and released before the next
RHS. Only one returned Y is held at a time; the private solver RHS copy is now
explicitly released after each call. Solver iteration history starts fresh for
each new RHS. The existing GMRES algorithm and admission mask are retained.

## Verification

Local integration fixtures passed for default nine-target lists, custom loops
and lists, entry-5 target preservation, rerun protection of CV inputs, shuffled
manifest rows, invalid/truncated grids, and the existing conduction subset/spin
staging paths. Perl syntax and Fortran structural parsing were checked.
A Fortran compiler and a complete runnable OCEAN dataset were unavailable here,
so compilation, MPI execution, and numerical comparison of Y against the earlier
single-target solver remain to be performed with the OCEAN build on your cluster.
