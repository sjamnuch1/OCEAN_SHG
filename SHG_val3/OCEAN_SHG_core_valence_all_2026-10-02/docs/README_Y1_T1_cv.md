# Separate core-valence `Y(2w)` and `T1_gv` programs

This add-on continues the standalone core-valence pathway after
`shg_rbxc_cv.x`.  It intentionally uses separate executable and manifest
names from the core-core pathway.

The implemented sequence is

\[
 B_{gv}^{bc}=R_{g\leftarrow v}^{b}X_v^c(\omega),
\]

\[
 [(\Omega+i\eta_g)I-H_g]Y_{gv}^{bc}(\Omega;\omega)=B_{gv}^{bc},
\]

and, for SHG at \(\Omega=2\omega\),

\[
 T_{1,gv}^{abc}(\omega)=
 \langle d_g^a\mid Y_{gv}^{bc}(2\omega;\omega)\rangle.
\]

Here `g` denotes the core-exciton sector and `v` the valence-exciton
sector.  The last resolvent and contraction act in the ordinary OCEAN core
XAS space, but their manifests remain tagged `CORE_FROM_VALENCE`/`GV` so
they cannot be confused with the `gg` pathway.

## Contents

- `SHG_CV2/shg_y1_cv.f90`: OCEAN-linked second GMRES solve.
- `SHG_CV2/shg_cv_manifest.f90`: strict CV manifest and basis checks.
- `SHG_CV2/shg_gmres.f90`: private solver copy used by the earlier core add-on.
- `SHG_CV2/shg_core_source.f90`: initializes the staged core Hamiltonian.
- `SHG_CV2/prepare_y1_cv.pl`: stages the GV solve in the completed RbXc directory.
- `SHG_CV2/stage_core_xas_cv.pl`: copy-only adaptation of OCEAN core staging.
- `SHG_CV2/verify_cv_core_basis.py`: checks structural conduction-space compatibility.
- `SHG_CV2/shg_t1_gv_contract.f90`: serial final contraction.
- `SHG_CV2/Makefile`: builds both executables inside an OCEAN source tree.
- `SHG_CV2/tests/test_interfaces.py`: standalone interface tests.

## Why a basis certificate is required

The CV vertex leaves the spectator electron in the conduction Bloch frame
used by the valence calculation.  The next GMRES solve interprets the same
coefficients in the conduction frame of the staged core calculation.
The band order, Bloch gauge, and k frame must physically agree. The lightweight
staging check can establish structural compatibility, but it cannot prove all
of these properties.

The supplied verifier compares the `con.u2.dat` byte counts, conduction-band
ranges, `ncon`, `nk`, and k meshes. It deliberately does not hash or read the
large wavefunction files. Consequently, a successful certificate means the
two spaces are structurally compatible; it does not prove bytewise equality or
an identical Bloch gauge. The user must ensure both files originated from the
same DFT/PREP conduction-wavefunction dataset. A failed structural comparison
still requires correcting the staging or applying an explicit basis
transformation.

For a finite-q valence source, the electron frame may differ from a q=0 core
calculation even when the structural checks pass. The GMRES driver also compares
the manifest's reduced q to
`sys%qinunitsofbvectors`.  Use a consistently staged electron frame or
implement the required overlap transformation.

## Installation and build

Copy `SHG_CV2` into the root of the same OCEAN source tree used for the core
calculation.  Its sibling directories must be `OCEAN2`, `Modules`, `CNBSE`,
and `SCREEN`, and `Makefile.arch` must already exist.  Build normal OCEAN
first, then:

```bash
cd SHG_CV2
make
```

This creates:

```text
shg_y1_cv.x
shg_t1_gv_contract.x
```

`shg_t1_gv_contract.x` is serial and has no OCEAN runtime dependency.  The
GMRES executable must be linked with the same compiler/options as OCEAN.

## 1. Stage the GV solve in the completed RbXc directory

Enter the directory in which `shg_rbxc_cv.x` generated `shg_rbxc_cv.info` and
the listed `RbXc` binaries. `--ocean-run-root` identifies the prepared OCEAN
calculation containing `Common`, `DFT`, `OPF`, `PREP`, and `SCREEN`; it does
not identify a previous core-response run:

```bash
cd /run/R_CV_zz

perl prepare_y1_cv.pl \
  --ocean-run-root /run/OCEAN_XAS \
  --component xzz \
  --edge-index 1 \
  --valence-con /run/PREP_VAL/BSE/con.u2.dat \
  --ocean-bin /path/to/ocean/bin
```

The current directory is the default work directory. `--work-dir DIR` may be
used when invoking the script from elsewhere; the old `--target DIR` spelling
is retained as a deprecated alias. The completed manifest and every listed
`RbXc` binary must already be local to this directory. The script:

1. Reads only the OCEAN computation directory's `Common`, `DFT`, `OPF`,
   `PREP`, and `SCREEN` ingredients needed to reconstruct \(H_g\).
2. Selects one edge and creates the core photon for the left index \(a\),
   corresponding to the final dipole \(d_g^a\).
3. Builds a clean core-XAS stage in a temporary child of the OCEAN calculation,
   validates it, and merges its files into the existing RbXc directory.
4. Preserves `shg_rbxc_cv.info` and all of its completed `RbXc` binaries; they
   are neither copied nor replaced.
5. Uses `--valence-con` to compare the valence/core file sizes, band range, and
   k mesh, then produces `cv_core_basis_match.info`. The large valence file is
   neither hashed nor copied into the Y1 directory.
6. Writes `cv_vectors.list`, including the required \(2\omega\) target for
   each local vector, plus `shg_y1_cv_staged_files.list` recording every file
   generated or replaced by the preparation step.

The OCEAN computation directory is therefore an immutable source of Hamiltonian
inputs, not a source of a previous core GMRES solution. No
`shg_cache_Gm_core_*` response is imported. Re-running the preparation step
replaces its staged core inputs, but it does not delete existing `Y1` or dipole
outputs. The batch solver replaces selected Y/dipole products on rerun; use
new output names to preserve products from a different setup.

The CV input is also self-contained at runtime: after staging, the valence
directory is not read by either Fortran executable.

## 2. Batch core resolvent

The 2026-10-02 release reads `bse2.in` and processes the full CV manifest by
original index after loading the core Hamiltonian once. See
[README_Y1_BATCH.md](README_Y1_BATCH.md) for full commands and input examples.

```bash
mpirun -np 4 /path/to/shg_y1_cv.x --component xzz
# Or start at the fifth original source AND fifth full-grid energy:
mpirun -np 4 /path/to/shg_y1_cv.x --component xzz --start-index 5
```

Keep the entire `bse2.in` grid. Per-vector outputs retain the `Y1_GV_COMPLETE`
manifest format. The default final dipole is `shg_cache_d_core_x.bin`.
Selected existing Y/dipole outputs are replaced; unselected Y files remain.

## 3. Contract with the core dipole

Use the physical core dipole generated by the same `shg_y1_cv.x` invocation:

```bash
./shg_t1_gv_contract.x \
  shg_cache_d_core_x.bin.info \
  shg_cache_d_core_x.bin \
  shg_cache_Y1_cv_zz_000001.bin.info \
  shg_cache_Y1_cv_zz_000001.bin \
  shg_T1_gv_xzz.dat
```

The program computes the compensated complex sum

\[
T_{1,gv}^{abc}=\sum_i d_{g,i}^{a*}Y_{gv,i}^{bc},
\]

with no extra \(1/N_k\), charge factor, symmetrization, or permutation sum.
The result is one raw time-ordering contribution, not the completed
\(\chi^{(2)}\).  It checks status, sectors, axes, dimensions, k mesh, edge,
site, layouts, file sizes, finite values, and recorded norms.

## Compatibility summary

| Interface | Required condition |
|---|---|
| CV `RbXc` to core GMRES | Compact complex128, no padding |
| Spectator conduction electron | Same-size file and matching band/k metadata; common DFT/PREP origin required |
| Band window | Same global conduction bounds |
| k sampling | Same `kmesh`, `nk`, ordering, and electron frame |
| Core sector | Same site, `ZNL`, `nalpha`, and staged XAS model space |
| Final contraction | Core dipole and `Y_gv` in identical OCEAN vtor order |
| Prior core response outputs | Never staged or read |
| Runtime valence-directory access | None after preparation |

The last two mathematical operations are the same as in the core-core
pipeline, but the old binaries are not treated as drop-in tools because the
input provenance and sector checks are different.

## Tests

The interface tests exercise successful certification, intentional basis
mismatch rejection, and the fresh Perl stager's validation-only path:

```bash
python3 tests/test_interfaces.py
```

After compiling, also test the serial complex contraction:

```bash
python3 tests/test_interfaces.py --contract-exe ./shg_t1_gv_contract.x
```

No Fortran compiler was available in the packaging environment, so the
OCEAN-linked executable was not compiled there.  Run these tests and a small
staged OCEAN calculation before production use.
