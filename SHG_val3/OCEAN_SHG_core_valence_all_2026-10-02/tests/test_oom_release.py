#!/usr/bin/env python3
"""Static regression checks for the memory-safe valence cache lifecycle."""
from pathlib import Path

source = (Path(__file__).resolve().parents[1] / 'shg_xc_minus_valence.f90').read_text()

clean = source.index('call SHG_gmres_clean(ierr)')
dipole_write = source.index('call write_vector_stream(sys, physical_dc')
loop = source.index('do i = 1, size(xc_minus)', dipole_write)
solution_write = source.index('call write_vector_stream(sys, xc_minus(i)', loop)
solution_kill = source.index('call OCEAN_psi_kill(xc_minus(i), release_ierr)', loop)
manifest = source.index('call write_manifest(', loop)
cleanup_start = source.index('do i = solutions_released + 1, size(xc_minus)')

assert clean < dipole_write
assert solution_write < solution_kill < manifest
assert cleanup_start > manifest
assert 'have_physical_dc = .false.' in source[dipole_write:loop]
assert 'solutions_released = i' in source[solution_kill:manifest]

print('PASS: GMRES workspace and materialized vectors are released in safe order')
