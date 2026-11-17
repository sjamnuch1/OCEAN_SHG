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
- `SHG_CV2/prepare_y1_cv.pl`: creates a new self-contained GV run directory.
- `SHG_CV2/stage_core_xas_cv.pl`: copy-only adaptation of OCEAN core staging.
- `SHG_CV2/verify_cv_core_basis.py`: certifies an identical conduction basis.
- `SHG_CV2/shg_t1_gv_contract.f90`: serial final contraction.
- `SHG_CV2/Makefile`: builds both executables inside an OCEAN source tree.
- `SHG_CV2/tests/test_interfaces.py`: standalone interface tests.

## Why a basis certificate is required

The CV vertex leaves the spectator electron in the conduction Bloch frame
used by the valence calculation.  The next GMRES solve interprets the same
coefficients in the conduction frame of the staged core calculation.
Equality of `ncon`, `nk`, and `nalpha` is therefore not enough: band order,
Bloch gauge, and the k frame must also agree.

The supplied verifier compares the valence and core `con.u2.dat` files by
byte count and SHA-256.  It also compares the conduction band range and
k mesh.  This is deliberately conservative.  A failed comparison means
that an explicit conduction overlap/basis transformation is needed before
the CV `RbXc` vector can be used; do not bypass the check merely because the
array dimensions agree.

For a finite-q valence source, exact matching to a q=0 core calculation will
normally fail.  The GMRES driver also compares the manifest's reduced q to
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

## 1. Create a fresh standalone GV directory

Run the preparation script from anywhere.  `--ocean-run-root` identifies the
prepared OCEAN computation directory containing `Common`, `DFT`, `OPF`,
`PREP`, and `SCREEN`; it does not identify a previous core-response run.  The
target must be a new direct child of that directory:

```bash
perl prepare_y1_cv.pl \
  --ocean-run-root /run/OCEAN_XAS \
  --target /run/OCEAN_XAS/SHG_GV_xzz \
  --component xzz \
  --edge-index 1 \
  --cv-manifest /run/R_CV_zz/shg_rbxc_cv.info \
  --valence-con /run/PREP_VAL/BSE/con.u2.dat \
  --ocean-bin /path/to/ocean/bin
```

The script deliberately has no `--force` mode.  It refuses an existing target,
so stale core response vectors cannot survive into the new calculation.  It:

1. Reads only the OCEAN computation directory's `Common`, `DFT`, `OPF`,
   `PREP`, and `SCREEN` ingredients needed to reconstruct \(H_g\).
2. Selects one edge and creates the core photon for the left index \(a\),
   corresponding to the final dipole \(d_g^a\).
3. Copies all core inputs, including `con.u2.dat`; no symbolic link remains.
4. Copies the completed CV manifest and every listed `RbXc` vector into the
   target.
5. Copies the valence `con.u2.dat` and produces
   `cv_core_basis_match.info` after exact SHA-256, band-range, and k-mesh
   verification.
6. Writes `cv_vectors.list`, including the required \(2\omega\) target for
   each copied vector.

The OCEAN computation directory is therefore an immutable source of Hamiltonian
inputs, not a source of a previous core GMRES solution.  No `shg_cache_Gm_core_*`,
old `RbXc`, old `Y1`, or old core dipole is imported.

The CV input is also self-contained at runtime: after staging, the valence
directory is not read by either Fortran executable.

## 2. Solve the core resolvent at `2w`

Enter the newly created directory—not the original core directory—and select
one line from `cv_vectors.list`:

```bash
cd /run/OCEAN_XAS/SHG_GV_xzz

mpirun -np 4 /path/to/shg_y1_cv.x \
  shg_rbxc_cv.info \
  shg_cache_RbXc_cv_zz_000001.bin \
  cv_core_basis_match.info \
  xzz \
  TARGET_2W_EV \
  shg_cache_Y1_gv_zz_000001.bin \
  CORE_ETA_EV \
  shg_cache_d_core_x_gv_000001.bin
```

Arguments are mandatory and ordered as:

```text
CV_MANIFEST CV_RBXC_BIN BASIS_MATCH_INFO COMPONENT
TARGET_ENERGY_EV Y_OUTPUT_BIN CORE_ETA_EV DIPOLE_OUTPUT_BIN
```

`COMPONENT=abc`; its `b` and `c` labels must match the CV manifest.  `a` is
recorded for the later core dipole contraction.  `TARGET_ENERGY_EV` uses the
same energy convention as the staged core OCEAN Hamiltonian.  For the SHG
point it must equal twice the source omega recorded in the CV manifest.

`CORE_ETA_EV` is explicit.  The program does not silently reuse the
first, valence-sector broadening.  It uses the GMRES depth, preconditioner,
and convergence tolerance from `bse.in`, but overrides that file's energy
grid with the single command-line target.

Before reading the compact stream, the code requires

```text
OCEAN_psi_size_full = ncon * nk * nalpha
```

and stops if the OCEAN vector has padded storage.  The RHS is then passed
through OCEAN's core admission mask.  GMRES solves a normalized RHS and the
program restores the physical norm to the output.

Outputs are the requested stream and `OUTPUT_BIN.info`.  The manifest status
is `Y1_GV_COMPLETE`, sector is `CORE_FROM_VALENCE`, and pathway is `GV`.

During the same invocation, `shg_y1_cv.x` reconstructs the physical
core-absorption source \(d_g^a\), writes `DIPOLE_OUTPUT_BIN`, and creates its
`XC_MINUS_COMPLETE` companion manifest.  This replaces the earlier requirement
to obtain the final dipole from a separate core `shg_xc_minus` run.  All output
files use `status='new'`; an existing file causes a stop instead of being
silently reused or overwritten.

## 3. Contract with the core dipole

Use the physical core dipole generated by the same `shg_y1_cv.x` invocation:

```bash
./shg_t1_gv_contract.x \
  shg_cache_d_core_x_gv_000001.bin.info \
  shg_cache_d_core_x_gv_000001.bin \
  shg_cache_Y1_gv_zz_000001.bin.info \
  shg_cache_Y1_gv_zz_000001.bin \
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
| Spectator conduction electron | Exact certified `con.u2.dat` match |
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
