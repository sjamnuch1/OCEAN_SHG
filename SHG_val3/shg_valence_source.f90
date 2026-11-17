! Standalone finite-q valence optical source for original OCEAN.
module SHG_valence_source
  use AI_kinds, only : DP
  implicit none
  private
  public :: SHG_read_component, SHG_prepare_valence_source
contains
  subroutine SHG_read_component(component, axis, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, MPI_INTEGER
    character(len=3), intent(out) :: component
    character(len=1), intent(out) :: axis
    integer, intent(inout) :: ierr
    character(len=1) :: labels(3)
    integer :: ios, bcast_ierr, j
    ios = 0
    if (myid == root) then
      open(93, file='shg_component', status='old', iostat=ios)
      if (ios == 0) then
        read(93, *, iostat=ios) labels
        close(93)
      endif
    endif
    call MPI_BCAST(ios, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (ios /= 0 .or. bcast_ierr /= 0) then
      ierr = 7801
      return
    endif
    call MPI_BCAST(labels, 3, MPI_CHARACTER, root, comm, ierr)
    if (ierr /= 0) return
    do j = 1, 3
      if (index('xyz',labels(j)) == 0) then
        ierr = 7802
        return
      endif
      component(j:j) = labels(j)
    enddo
    axis = labels(3)
  end subroutine

  subroutine SHG_prepare_valence_source(sys, axis, normalized_d, source_norm, ierr)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use OCEAN_constants, only : PI_DP
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_INTEGER
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_nrm, OCEAN_psi_full2min, &
                         OCEAN_psi_new, OCEAN_psi_zero_min, OCEAN_psi_axpy, &
                         OCEAN_psi_copy_min, OCEAN_psi_kill
    use OCEAN_system, only : o_system
    type(o_system), intent(inout) :: sys
    character(len=1), intent(in) :: axis
    type(OCEAN_vector), intent(inout) :: normalized_d
    real(DP), intent(out) :: source_norm
    integer, intent(inout) :: ierr
    type(OCEAN_vector) :: rotated
    real(DP) :: qcart(3), qlength, qperp, qsign, nrm2
    integer :: idir, ios, bcast_ierr, tm(4), j, cleanup_ierr
    integer(int64) :: expected_bytes, file_bytes
    logical :: exists

    source_norm = 0.0_DP
    idir = index('xyz',axis)
    if (idir == 0) then
      ierr = 7803
      return
    endif
    qcart = matmul(sys%bvec, sys%qinunitsofbvectors)
    qlength = sqrt(sum(qcart**2))
    if (.not. ieee_is_finite(qlength) .or. qlength <= tiny(1.0_DP)) then
      if (myid == root) write(6,*) 'Finite nonzero PREP photon_q is required.'
      ierr = 7804
      return
    endif
    !do j = 1, 3
    !  if (j == idir) cycle
    !  if (abs(qcart(j)) > 1.0e-8_DP*qlength) then
    !    if (myid == root) write(6,*) 'PREP q must be parallel to the requested Cartesian axis.'
    !    ierr = 7805
    !    return
    !  endif
    !enddo
    !qsign = sign(1.0_DP,qcart(idir))

    ! Magnitude of the component perpendicular to the requested axis.
    qperp = 0.0_DP
    do j = 1, 3
     if (j == idir) cycle
     qperp = qperp + qcart(j)**2
    enddo
    qperp = sqrt(qperp)

    ! Continue with a warning rather than terminating.
    if (qperp > 1.0e-8_DP*qlength) then
     if (myid == root) then
      write(6,'(A)') &
         'WARNING: PREP q is not parallel to the requested Cartesian axis.'
      write(6,'(A,A1)') &
         '         Requested axis: ', axis
      write(6,'(A,3(1X,ES16.8))') &
         '         q_cart [bohr^-1]:', qcart
      write(6,'(A,1X,ES16.8)') &
         '         |q_perp|/|q|:', qperp/qlength
      write(6,'(A)') &
         '         Source is sign(q_axis)*(qhat dot d), not pure d_axis.'
     endif
    endif

    qsign = sign(1.0_DP,qcart(idir))

    ! Validate before the original reader enters collective MPI file reads.
    ! tmels.info: number of stored valence bands (indexed from 1),
    !             first/last stored conduction band, number of k points.
    ios = 0
    if (myid == root) then
      open(93, file='tmels.info', status='old', iostat=ios)
      if (ios == 0) then
        read(93, *, iostat=ios) tm
        close(93)
        if (ios == 0) then
          if (tm(1) < sys%brange(2) .or. tm(2) > sys%brange(3) .or. &
              tm(3) < sys%brange(4) .or. tm(4) /= sys%nkpts) ios = 1
          if (any(sys%brange < 1) .or. sys%brange(1) > sys%brange(2) .or. &
              sys%brange(3) > sys%brange(4)) ios = 1
          if (sys%brange(1) < sys%file_brange(1) .or. sys%brange(2) > sys%file_brange(2) .or. &
              sys%brange(3) < sys%file_brange(3) .or. sys%brange(4) > sys%file_brange(4)) ios = 1
          if (ios == 0) then
            expected_bytes = int(storage_size(cmplx(0.0_DP,0.0_DP,DP))/8,int64) * &
                 int(tm(1),int64) * int(tm(3)-tm(2)+1,int64) * int(tm(4),int64) * int(sys%nspn,int64)
            inquire(file='ptmels.dat', exist=exists, size=file_bytes)
            if (.not. exists) then
              ios = 1
            else if (file_bytes /= expected_bytes) then
              write(6,*) 'ptmels.dat byte count: actual, expected ', file_bytes, expected_bytes
              ios = 1
            endif
          endif
        endif
      endif
    endif
    call MPI_BCAST(ios, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (ios /= 0 .or. bcast_ierr /= 0) then
      if (myid == root) write(6,*) 'Incompatible band window or finite-q source files.'
      ierr = 7806
      return
    endif

    ! tmel_selector=1 reads S=conjg(M)/|q| in original OCEAN.
    ! It then applies the occupation/energy mask, normalizes S, and loads H_V.
    call OCEAN_load_data(sys, normalized_d, ierr)
    if (ierr /= 0) return
    source_norm = sqrt(max(0.0_DP,normalized_d%kpref) * real(sys%nkpts,DP) * &
                       sys%celvol**2 / (4.0_DP*PI_DP))
    if (.not. ieee_is_finite(source_norm) .or. source_norm <= tiny(1.0_DP)) then
      ierr = 7807
      return
    endif
    call OCEAN_psi_nrm(nrm2, normalized_d, ierr)
    if (ierr /= 0) return
    if (.not. ieee_is_finite(nrm2) .or. abs(nrm2-1.0_DP) > 1.0e-8_DP) then
      if (myid == root) write(6,*) 'Invalid normalized valence-source norm squared: ', nrm2
      ierr = 7808
      return
    endif

    ! The phase is essential for coherent SHG contractions:
    ! d^c = -i*conjg(M)/q_c = -i*sign(q_c)*S.
    ! Use public vector operations so OCEAN handles FULL/MIN validity.
    call OCEAN_psi_full2min(normalized_d, ierr)
    if (ierr /= 0) return
    call OCEAN_psi_new(rotated, ierr)
    if (ierr /= 0) return
    rotated%kpref = normalized_d%kpref
    call OCEAN_psi_zero_min(rotated, ierr)
    if (ierr == 0) call OCEAN_psi_axpy(0.0_DP, normalized_d, rotated, ierr, ival=-qsign)
    if (ierr == 0) call OCEAN_psi_copy_min(normalized_d, rotated, ierr)
    cleanup_ierr = 0
    call OCEAN_psi_kill(rotated, cleanup_ierr)
    if (ierr == 0) ierr = cleanup_ierr
  end subroutine
end module SHG_valence_source
