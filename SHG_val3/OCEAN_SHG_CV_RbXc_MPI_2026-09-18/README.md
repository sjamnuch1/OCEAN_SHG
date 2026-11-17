# Core <- valence R_b X_c standalone

This is the next step of the core–valence T1 pathway started with
`shg_xc_minus_valence.x`. It computes **R_CV X_V**, producing a core-hole
vector. It is not the valence–valence operator R_VV. The default Fortran
executable uses MPI to distribute contiguous conduction-band blocks. It uses
standard MPI directly and requires no OCEAN modules, FFT, or BLAS libraries.
The original serial implementation remains available as a reference target.

## Operator and conventions

Use determinants |n g^-1>=a_n^dagger a_g|0> and
|n' v^-1>=a_n'^dagger a_v|0>. The Slater–Condon one-body rule gives

    <n g^-1|R_b|n' v^-1> = -delta(n,n') <v|r_b|g>.

The electron term vanishes because the core and valence holes are orthogonal.
The minus sign is the hole term of r_e-r_h. Thus

    B(n,k,alpha) = (R_CV^b X_V^c)(n,k,alpha)
                 = -sum_v d(v,k,m,h) X(n,v,k,beta(h,e)).

Variables:

| Symbol | Meaning |
|---|---|
| n | Local conduction index, 1..ncon; global bounds from input brange |
| v | Local active valence index, 1..nval |
| k | OCEAN flattened k index, 1..nk |
| g=(m,h) | Localized core orbital and its spin |
| m | Core magnetic index ordered from -l through +l |
| h,e | Core-hole/valence-hole spin and spectator-electron spin, each 1,2 |
| alpha | ((h-1)*(2*l+1)+(m+l))*2+e, OCEAN core ordering |
| beta | 2*(h-1)+e for an explicit four-channel valence vector |
| p | OPF projector index, 1..nproj |
| P(p,v,k,h) | Complex valence projector coefficient in parcksv |
| M(p,m) | Complex core/projector dipole from the selected mels file |
| E(v,k,m,h) | sum_p M(p,m)*P(p,v,k,h), OCEAN's xes_vec |
| d(v,k,m,h) | conjg(E), corresponding to <v|r_b|g> |

`OCEAN_rixs_holder.f90:cksv_read` constructs E by the ordinary complex
product M*P, without conjugating either factor. Its RIXS routine takes
core vectors to valence vectors. Here the direction is reversed, giving
conjg(E), and the hole operator supplies the minus sign. The code uses
SUM(M*P), not complex DOT_PRODUCT, because Fortran DOT_PRODUCT would
conjugate M. `OCEAN_psi.f90:OCEAN_psi_dotter` independently uses conjg(M*P)
for the core absorption source.

The historical RIXS routine multiplies its result by nspin and explicitly
comments that this is a legacy intensity convention. This program does
not copy that factor into the physical one-body operator. It also does
not normalize the source or output, or add electron charge factors.
With length-dipole mels in bohr and physical input X in bohr/Hartree,
the output is in bohr squared/Hartree.

For nbeta=4 the spin channels are explicit. For nbeta=1, the previous
valence driver supplies the scalar response per spin. This program
replicates it into h=e=1 and h=e=2, leaving the two spin-flip channels
zero. It does not divide by sqrt(2): this is per-spin replication, not
an embedding of a separately normalized singlet. A differently normalized
singlet input must first be converted to this convention. With nspin=1
and nbeta=4, the same spatial projectors serve both hole spins.

## Required inputs

1. Completed `shg_xc_minus_valence_*.info` and all listed X dumps.
2. `parcksv.<element><site>` for the selected core site and the **same
   valence Bloch states/gauge** as that solve. This file is produced by
   OCEAN PREP's `ocean_cks` when valence site projectors are generated,
   as in the RIXS preparation path. A plain optical VAL preparation may
   not have generated it. Generate the site projectors from the same
   underlying wavefunctions; the GMRES dump alone cannot determine R_CV.
3. The matching `mels.zZZZnNNlLLpPP` file for Cartesian polarization b.
   Select an electric length-dipole photon definition, not quadrupole or
   a differently normalized operator. The integer PP is OCEAN's photon
   file number; it is not universally x=1,y=2,z=3.

The reader implements the inspected original OCEAN parcksv format:
three native int32 values (nproj, nval_stored*nk, nspin), followed by
native complex128 data in Fortran order
(nproj,nval_stored,nk,nspin). Legacy sequential cksv files and int64
header builds are not supported. The mels file contains two real values
per line, projector index fastest and core m slowest.

You must supply the first global valence band stored in parcksv. Its
header contains a count, not the global lower bound. Use the brange
metadata from the PREP run that created that file. The active slice is
taken from the input X manifest; no conduction-band sum or truncation
is performed. Changing the conduction window in X changes ncon in B.

## Build and run

```bash
make MPIFC=mpifort
python3 prepare_rbxc_cv.py \
  --manifest /run/SHG_VAL_N40/shg_xc_minus_valence_z.info \
  --parcksv /matching_projectors/parcksv.Ti0001 \
  --projector-first-valence 1 \
  --mels /matching_dipoles/mels.z022n01l00p03 \
  --znl 22 1 0 --element Ti --site 1 --axis z \
  --out /run/R_CV_zz_N40

cd /run/R_CV_zz_N40
mpirun -np 8 /path/to/shg_rbxc_cv.x cv.in
```

The example assumes photon file 03 defines z-polarized electric dipole
matrix elements for Ti K. Replace the paths, site, edge, photon number,
and stored band lower bound with your data.

Preparation reads the valence manifest and schedules **every** listed
energy. It validates dimensions, file sizes, source units, band slices,
and finite numerical metadata. It writes cv.in and provenance.json into
a fresh directory. The executable loads the OPF/dipole contraction once,
then reads and contracts one X at a time. Both stages reject overwriting
existing outputs. If a run fails, use a new directory for a clean rerun.

The MPI implementation partitions `ncon` into contiguous blocks. For each
energy it processes one k point at a time. MPI-IO subarray views read only
`X(n_first:n_last,:,k,:)` and write only
`RbXc(n_first:n_last,k,:)`; the complete X and output vectors are never
replicated on a rank. This is intended for slabs and other low-dimensional
systems with relatively few k points but many vacuum-derived conduction
bands. At most `ncon` ranks are active. Supplying more ranks is harmless but
the extras wait in `MPI_Finalize`.

The output byte layout is unchanged from the serial program, so downstream
readers do not need an MPI-specific format. Running the MPI executable without
`mpirun`, or with `mpirun -np 1`, is also supported.

To build the retained serial reference instead:

```bash
make FC=gfortran shg_rbxc_cv_serial.x
```

## Outputs and the next core solve

`shg_cache_RbXc_cv_zz_000001.bin`, etc. are compact stream complex128
arrays with shape **(ncon,nk,4*(2*l+1))**, ncon fastest. There are no
Fortran record markers and no OCEAN padding. `shg_rbxc_cv.info` contains
energy, norm, shape, edge/site, spin convention, and momentum metadata.
Treat the run as complete only when its last completion record is
`RBXC_CV_COMPLETE`; an interrupted manifest starts with status RUNNING.

The electron remains at the **input k+q**, with its original band bounds,
phases, and ordering. This operator changes the hole only. In particular,
these outputs must not be passed as if they were in an independently
generated q=0 core conduction basis. A consistent subsequent core solve
must use the same electron frame, or an explicit conduction overlap
transformation must be applied first. Equal dimensions alone do not
establish that correspondence. The q->0 approximation also requires
convergence; merely deleting q metadata is not a basis transformation.

This output manifest is intentionally distinct from the older core–core
manifest. The previous `shg_y1_core` reader must be adapted to accept
this sector, verify its electron-frame compatibility, and translate its
physical omega to the core solver's energy reference. It is not a drop-in
input to that older driver. No core occupation/admission mask is applied
here; apply the target core model-space mask when loading its RHS.

Site identity, photon polarization, radial normalization, actual k
coordinates, and Bloch gauge are not encoded in the raw parcksv header.
The program cannot certify them from dimensions. Preserve the originating
PREP/photon metadata and choose matching inputs.

## Parallel decomposition and cost

The local dipoles are built once on rank zero at cost
O(nproj*nval*nk*(2*l+1)*nspin), then the smaller `d(nval,nk,nm,nspin)` array
is broadcast. Each X costs O(ncon*nval*nk*nalpha), with small matrix-vector
products over valence bands. A rank owning `ncon_local` bands stores only
`X(ncon_local,nval,nbeta)` and `RbXc(ncon_local,nalpha)` for one k point.

No neighboring-k data, conduction-link matrices, FFT, GMRES call, or
inter-rank reduction of the physical vector is needed. Communication is
limited to setup broadcasts, collective MPI-IO, and one scalar norm reduction
per energy. Parallel efficiency is therefore controlled mainly by the number
of conduction bands and the MPI filesystem performance.

## Verification

```bash
python3 test_cv.py                     # NumPy needed for tests only
python3 test_cv.py --exe ./shg_rbxc_cv.x # singleton MPI execution
python3 test_cv.py --exe ./shg_rbxc_cv.x --mpiexec /path/to/mpirun --np 4
```

Synthetic preparation tests cover multiple energies, active valence
slicing, three spin layouts, overwrite rejection, and reconstruction from
1, 2, 4, and more requested ranks than conduction bands. The optional
compiled test compares every output component to an independently
assembled dense matrix for genuinely complex inputs, which detects
conjugation, sign, index, spin-channel, MPI file-view, and distributed-write
errors. It also checks an adjoint identity. Fortran source parsing and the
preparation/decomposition tests passed in the development environment. No
MPI Fortran compiler was available there: the executable has **not** been
compiled or tested on an actual OCEAN dataset. Run the supplied multi-rank
compiled test and then validate the local dipoles against OCEAN xes_vec
before production SHG use.

Source grounding: inspected `OCEAN_SHG_core_standalone_original_2026-08-19`,
especially OCEAN2/OCEAN_rixs_holder.f90 (cksv_read, rixs_seed),
OCEAN2/OCEAN_psi.f90 (OCEAN_psi_dotter), and PREP/src/ocean_cks.f90
(parallel projector-file header and layout). This is a fresh small
implementation of those contractions, not a linked OCEAN executable.
