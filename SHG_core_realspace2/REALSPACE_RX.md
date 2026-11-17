# Standalone real-space `R X` implementation

`shg_rbxc_core_realspace.x` computes the same finite-mesh core-exciton
pair-position action as the neighboring-`k` implementation,

\[
q_{b,c}(\omega)=R^b_{\mathrm{core}}X^-_{\mathrm{core},c}(\omega),
\]

but applies the operator in the Born--von Karman real-space representation.
It is intended as an independent validation path. The link implementation is
still the recommended production path for repeated sources.

## Equations implemented

For fixed core channel \(\alpha\) and primitive-cell grid point \(r_a\), the
program first contracts the conduction bands:

\[
F_{k\alpha}(r_a)=\sum_{n=1}^{N_c}u_{nk\sigma_e}(r_a)X_{nk\alpha}.
\]

`con.u2.dat` stores the cell-periodic function \(u\), so the within-cell Bloch
phase is restored:

\[
\Phi_{k\alpha}(r_a)=e^{ik\cdot r_a}F_{k\alpha}(r_a).
\]

For reduced mesh indices \(m_j=0,\ldots,N_j-1\), OCEAN's point is

\[
k=\sum_j\frac{k_{0j}+m_j}{N_j}b_j.
\]

The core operator remains the optical \(q=0\) operator. Nonzero upstream
`photon_q` is allowed: OCEAN QE/PREP places the occupied states at \(k-q\),
while `con.u2.dat` contains the base-mesh states at \(k\). The ordinary
`k0.ipt` mesh offset is retained. See `NONZERO_PHOTON_Q_2026-10-22.md` for
the basis requirements, changes, and verification limits.

The positive-exponent transform to supercell translations is

\[
\Psi_\alpha(r_a+R_n)=\frac{1}{N_k}\sum_k
e^{ik\cdot R_n}\Phi_{k\alpha}(r_a).
\]

OCEAN calls this `OCEAN_BACKWARD`; its FFT wrapper supplies the \(1/N_k\)
normalization. The inverse operation is the unnormalized negative-exponent
`OCEAN_FORWARD` transform, so the pair is exactly inverse.

The finite-mesh periodic relative coordinate is

\[
\rho^b_{\mathrm{FD}}(r_a+R_n)
=\sum_{j=1}^3 2w_j^b
\sin\!\left[
\frac{2\pi}{N_j}(n_j+\xi_{a,j}-\tau_j)
\right],
\]

where

\[
w_j^b=\frac{N_j a_j^b}{4\pi},\qquad
r_a=\sum_j\xi_{a,j}a_j,\qquad
r_h=\sum_j\tau_j a_j.
\]

The program multiplies

\[
Y_\alpha(r_a+R_n)=\rho^b_{\mathrm{FD}}(r_a+R_n)
\Psi_\alpha(r_a+R_n),
\]

transforms \(Y\) back to \(k\), and projects with the conjugate full Bloch
state:

\[
(R^bX)_{nk\alpha}=\sum_a
\left[e^{ik\cdot r_a}u_{nk\sigma_e}(r_a)\right]^*
\widetilde\Phi_{k\alpha}(r_a).
\]

No explicit cell-volume factor is inserted: `SHG_load_conduction_states`
normalizes every staged state with the same discrete grid inner product used
by the final projection.

With identical states and boundary phases, the result is algebraically equal
to

\[
i\sum_jw_j^b(T_j^+-T_j^-)X.
\]

The included link builder restores the reciprocal-lattice phase
\(e^{-iG_j\cdot r_a}\) when a forward neighbor crosses a Brillouin-zone face.
That correction is required for a point-by-point comparison with the
full-Bloch real-space route.

## Array shapes and ordering

| Object | Logical shape | Meaning |
|---|---:|---|
| `X` | `(Ncon,Nk,Nalpha)` | Input core-exciton vector |
| `states` | `(Nx_local,Ncon,Nk,Nspin)` | Distributed `con.u2.dat` states |
| `bloch_phase` | `(Nx_local,Nk)` | \(e^{ik\cdot r_a}\) |
| `rho_fd` | `(Nx_local,Nk)` | One value per primitive point and supercell translation |
| `scratch` | `(Nk)` | One \(k/R\) field passed through the 3-D FFT |
| `local_result` | `(Ncon,Nk,Nalpha)` | Rank-local backprojection before `MPI_ALLREDUCE` |

The product `k` index has `k3` fastest, then `k2`, then `k1`. The newer
`con.u2.dat` stream has primitive-grid `x` fastest, then `y`, then `z`. The FFT
plan is therefore initialized with `(N3,N2,N1)` so its linear storage matches
the OCEAN vector index.

Every \(\alpha\) channel is processed independently. The program never sums
core-hole channels into an electron density and never normalizes an incoming
SHG vector.

## Build

Build the unmodified OCEAN tree first, then from this add-on root run:

```bash
cp Makefile.arch.example Makefile.arch
# Edit Makefile.arch for the local compiler, MPI, BLAS, and FFTW setup.
make
make -C SHG shg_rbxc_core_realspace.x
```

The new executable links the original OCEAN objects and `FFT_wrapper.o`; no
source under `OCEAN2/` is modified.

## Run

In the staged core-XAS `SHG` directory, after `shg_xc_minus.x` has produced its
manifest and raw complex caches:

```bash
mpirun -np 4 /path/to/OCEAN/SHG/shg_rbxc_core_realspace.x \
  shg_xc_minus_core_z.info > shg_rbxc_core_realspace.log 2>&1
```

For component `xyz`, the output is

```text
shg_cache_Rrs_y_Gm_core_w_z_000001.bin
shg_cache_Rrs_y_Gm_core_w_z_000002.bin
...
shg_rbxc_core_realspace_yz.info
```

`Rrs` keeps these files distinct from the link result `R`. Both use raw
`complex(DP)` stream data in `OCEAN_psi_vtor` order.

## Compare with the link route

Run `shg_rbxc_core.x` and the real-space executable on the same staged data,
then compare corresponding caches:

```bash
python3 /path/to/OCEAN/SHG/compare_rx_caches.py \
  shg_cache_R_y_Gm_core_w_z_000001.bin \
  shg_cache_Rrs_y_Gm_core_w_z_000001.bin
```

The relative difference should be near floating-point/FFT tolerance. A larger
difference usually indicates inconsistent `xmesh.ipt`, `kmesh.ipt`,
`con.u2.dat`, band ordering, a missing Brillouin-zone boundary phase, or a
different admission mask.

The self-contained algebraic test exercises the same ordering, phases, and
nonorthogonal Cartesian weights without requiring an OCEAN run:

```bash
python3 SHG/test_realspace_equivalence.py
```

## Cost and intended use

The constant workspace (`bloch_phase`, `rho_fd`, and the FFT plan) is built
once and reused for every energy. The source-dependent work still requires,
for every \((r_a,\alpha)\), two band contractions and two FFTs. Its cost is
therefore substantially higher than applying precomputed small link matrices,
and it offers no accuracy improvement when both representations use the same
finite-mesh operator. Use it for validation, small systems, and development of
genuinely local real-space vertices; use `shg_rbxc_core.x` for production.
