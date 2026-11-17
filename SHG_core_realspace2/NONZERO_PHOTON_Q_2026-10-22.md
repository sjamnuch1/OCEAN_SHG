# Reuse conduction states from a nonzero-photon-q PREP — 2026-10-22

Compilation correction: the real-space manifest writer uses the public
read-only `SHG_realspace_get_k0` accessor. Workspace components remain private.
Update both `shg_core_vertex_realspace.f90` and
`shg_rbxc_core_realspace.f90`, then rebuild `shg_rbxc_core_realspace.x`.

`shg_rbxc_core.x` now accepts nonzero upstream `photon_q`. The shared validator
also enables `shg_rbxc_core_realspace.x`. Both continue to evaluate the optical
core operator R = r_e - r_h in the base-mesh conduction basis.

## Mesh convention

In the OCEAN QE/PREP implementation traced for this change, the high-band
conduction states are on the base mesh

    k_c(m) = sum_j [(m_j + k0_j)/N_j] b_j,

while occupied states are on k_v(m) = k_c(m) - q, up to reciprocal-lattice
wrapping. PREP writes the base-mesh conduction functions to `con.u2.dat`, their
energies to `wvfcninfo`, and the matching conduction projectors to `parcksc.*`.
`k0.ipt` remains the ordinary mesh offset. It need not be zero.

`sys%kshift` is derived from `qinunitsofbvectors.ipt`; it indicates nonzero
upstream photon-q metadata, rather than selecting a different part of the
already separated `con.u2.dat` stream. There are no extra shifted/unshifted
blocks to skip inside this conduction file.

## Source changes

- `shg_core_vertex.f90`: removes the blanket rejection of `sys%kshift` in
  `SHG_validate_xc_manifest`. Nonfinite q metadata is rejected instead. For
  nonzero q, the root rank reports that base-mesh conduction data are being
  used. Existing core-run, edge/site, photon/component, spin, mesh, dimension,
  vector-file, wavefunction byte-count, and state-norm checks remain in place.
- `shg_rbxc_core.f90`: adds output provenance for the base conduction mesh,
  active band range, upstream photon q, zero conduction photon-q shift, and
  zero core-operator momentum transfer. The link construction, reciprocal-face
  phases, hole phases, Cartesian weights, admission mask, binary format, and
  command-line interface are retained.
- `shg_core_vertex_realspace.f90`: builds the Bloch phase using
  `(k0 + mesh_index)/kmesh`, without adding upstream photon q.
- `shg_rbxc_core_realspace.f90`: records the same provenance and the actual
  workspace k0 in its output manifest.
- `test_realspace_equivalence.py`: checks link/FFT equivalence, Hermiticity,
  and cancellation of a uniform common-q phase for three phase offsets.
- `REALSPACE_RX.md` and `README.md`: document the base-mesh convention and use.

The global OCEAN `sys%kshift` and `sys%qinunitsofbvectors` values are preserved.
The q=0 core operator is evaluated using conduction data prepared alongside a
finite-q valence calculation. This change does not implement finite-q core SHG.

## Required basis compatibility

Use the same PREP conduction dataset for the incoming core X vector, core
projectors, energies, and `con.u2.dat`. A valence X vector is not a core X vector.
Equal dimensions cannot prove equality of Bloch gauges; the existing metadata
checks do not perform a content hash or recover a mesh from the headerless
wavefunction stream. Replacing the file with states from an independent DFT
calculation is not covered by this update.

The existing loader expects the stored conduction window described by
`brange.ipt` and reads its first `sys%num_bands` states. This release does not
change band-window selection or the wavefunction stream layout.

The staged runtime must still be one XAS core run with the required core
projectors, energies, and SCREEN ingredients. `prepare_shg.pl` and
`stage_core_xas.pl` already pass nonzero photon-q metadata through; they need
no change for such a core stage. Their upstream `calc.mode=xas` requirement
remains. Staging directly from a `calc.mode=val` or `rxs` root requires a
separate staging adaptation.

## Build and run

Install the source directory below your built OCEAN source tree, as with the
working Oct22 package. Build using the same compiler and MPI implementation:

```bash
make -C /path/to/OCEAN/SHG_core_realspace \
  shg_rbxc_core.x shg_rbxc_core_realspace.x
```

This source package excludes the old `.o`, `.mod`, and `.x` files so the old
executables cannot be mistaken for the revised build. If updating an existing
build directory, copy all four changed Fortran files and rebuild both targets.

In the staged core directory containing your completed core X manifest:

```bash
mpirun -np 4 /path/to/OCEAN/SHG_core_realspace/shg_rbxc_core.x \
  shg_xc_minus_core_z.info > shg_rbxc_core.log 2>&1

mpirun -np 4 /path/to/OCEAN/SHG_core_realspace/shg_rbxc_core_realspace.x \
  shg_xc_minus_core_z.info > shg_rbxc_core_realspace.log 2>&1
```

Keep the upstream `qinunitsofbvectors.ipt` value. No new runtime flag is needed.
Output filenames and the existing `RbXc` records retain their previous formats;
the new manifest keys are additional provenance.

## Verification

Run the self-contained numerical regression:

```bash
python3 /path/to/OCEAN/SHG_core_realspace/test_realspace_equivalence.py
```

It uses a nonorthogonal cell, nonzero k0, reciprocal-face transport, and a
3x4x5 mesh. It checks the original zero-phase case and two nonzero common-phase
cases against a relative tolerance of 1e-13. These are algebraic tests with
synthetic states, not a production Fortran/MPI execution.

For production validation, run both representations on identical states and
core X vectors, and compare their caches using `compare_rx_caches.py`. To
isolate this metadata change, compare the revised program with zero and
nonzero photon-q metadata while keeping the same conduction data and X.
Retain separate output directories for the comparison. Differences should be
at floating-point/MPI/FFT tolerance. Comparing independently generated DFT
wavefunctions element by element requires aligning their Bloch gauges first.

No Fortran compiler or complete runnable OCEAN dataset was available in the
editing environment. Compilation and a physical OCEAN/MPI run must therefore
be checked with your installation.
