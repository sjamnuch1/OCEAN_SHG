# Changes — 2026-10-02

- `shg_y1_cv.f90`: batch option interface, `bse2.in`, full-grid count checks,
  stable index pairing, start/end index and start-file selection, one-time
  Hamiltonian/dipole initialization, one RHS/Y at a time, compatible per-vector
  products, aggregate completion manifest, selected-output replacement.
- `shg_gmres.f90`: optional setup filename (default remains `bse.in`),
  parser error/finite-value checks, and explicit private RHS cleanup.
- `prepare_y1_cv.pl`: generate the full `bse2.in` grid, optional custom input
  or energy loop, optional core eta, original-index sorting/checks, actual
  target list and stage metadata. Existing conduction-range/mels staging retained.
- Documentation: updated commands, energy example, restart semantics, output
  behavior, and validation limits.

Other preparation scripts, the vertex program, first-solver driver, and final
contraction program require no changes. `stage_core_xas_cv.pl` still creates
`bse.in` for OCEAN core headers; the wrapper creates the separate `bse2.in`.
