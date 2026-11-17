! Copyright (C) 2026 OCEAN collaboration
!
! Standalone wrappers around the core/XAS source path in original OCEAN.
! No change to OCEAN2 is required: OCEAN_load_data constructs and normalizes
! the core dipole source and initializes the requested core BSE Hamiltonian.
module SHG_core_source
  use AI_kinds, only : DP

  implicit none
  private

  public :: SHG_prepare_core_source, SHG_read_component

contains

  subroutine SHG_read_component(component, right_axis, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, MPI_INTEGER

    character(len=3), intent(out) :: component
    character(len=1), intent(out) :: right_axis
    integer, intent(inout) :: ierr

    character(len=1) :: labels(3)
    integer :: ios, ierr_bcast, j

    labels = '?'
    ios = 0
    if (myid .eq. root) then
      open(unit=93, file='shg_component', form='formatted', status='old', iostat=ios)
      if (ios .eq. 0) then
        read(93, *, iostat=ios) labels
        close(93)
      endif
      if (ios .ne. 0) ierr = 7301
    endif

    call MPI_BCAST(ierr, 1, MPI_INTEGER, root, comm, ierr_bcast)
    if (ierr .ne. 0) return
    if (ierr_bcast .ne. 0) then
      ierr = ierr_bcast
      return
    endif
    call MPI_BCAST(labels, 3, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return

    do j = 1, 3
      labels(j) = SHG_lower_axis(labels(j))
      if (labels(j) .eq. '?') then
        if (myid .eq. root) write(6,'(A)') 'shg_component entries must be x, y, or z.'
        ierr = 7302
        return
      endif
      component(j:j) = labels(j)
    enddo
    right_axis = labels(3)
  end subroutine SHG_read_component


  !> Load the physical core-to-conduction dipole source using unchanged OCEAN,
  !! retain its pre-normalization norm, and return OCEAN's unit solver source.
  subroutine SHG_prepare_core_source(sys, normalized_d, source_norm, ierr)
    use OCEAN_constants, only : PI_DP
    use OCEAN_mpi, only : myid, root
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_nrm
    use OCEAN_system, only : o_system

    type(o_system), intent(inout) :: sys
    type(OCEAN_vector), intent(inout) :: normalized_d
    real(DP), intent(out) :: source_norm
    integer, intent(inout) :: ierr

    real(DP) :: normalized_norm_squared

    source_norm = 0.0_DP
    if (.not. associated(sys%cur_run)) then
      ierr = 7310
      return
    endif
    if (.not. sys%cur_run%have_core .or. sys%cur_run%calc_type .ne. 'XAS') then
      if (myid .eq. root) write(6,'(A)') 'The standalone core source requires one XAS run.'
      ierr = 7311
      return
    endif

    ! Original OCEAN performs all core-specific work here: psi allocation,
    ! OCEAN_psi_dotter, core energy loading, transition masking, unit
    ! normalization, multiplet initialization, and long-range initialization.
    call OCEAN_load_data(sys, normalized_d, ierr)
    if (ierr .ne. 0) return

    normalized_norm_squared = 0.0_DP
    call OCEAN_psi_nrm(normalized_norm_squared, normalized_d, ierr)
    if (ierr .ne. 0) return
    if (normalized_norm_squared .le. tiny(1.0_DP)) then
      if (myid .eq. root) write(6,'(A)') 'The normalized core dipole source has zero norm.'
      ierr = 7312
      return
    endif

    ! core_pnorm stores the physical pre-normalization norm in kpref:
    ! kpref = 4*pi*||d||^2 / (Nk*Omega^2).
    source_norm = sqrt(max(0.0_DP, normalized_d%kpref) &
                      * real(sys%nkpts,DP) * sys%celvol**2 / (4.0_DP*PI_DP))
    if (source_norm .le. tiny(1.0_DP)) then
      if (myid .eq. root) write(6,'(A)') 'The physical core dipole source has zero norm.'
      ierr = 7313
      return
    endif

    if (abs(normalized_norm_squared-1.0_DP) .gt. 1.0e-8_DP) then
      if (myid .eq. root) then
        write(6,'(A,ES16.8)') 'Unexpected normalized core-source norm squared: ', &
                              normalized_norm_squared
      endif
      ierr = 7314
      return
    endif
  end subroutine SHG_prepare_core_source


  pure function SHG_lower_axis(axis) result(lower)
    character(len=1), intent(in) :: axis
    character(len=1) :: lower

    select case(axis)
      case('x','X')
        lower = 'x'
      case('y','Y')
        lower = 'y'
      case('z','Z')
        lower = 'z'
      case default
        lower = '?'
    end select
  end function SHG_lower_axis

end module SHG_core_source

