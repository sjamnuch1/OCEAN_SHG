! Copyright (C) 2026 OCEAN collaboration
!
! Standalone core first-layer SHG test:
!   |X_core,c^-(omega)> = [(omega + i*eta) - H_BSE^core]^{-1} |d_core^c>
!
! This program uses only the public machinery in an unmodified OCEAN checkout.
! Every SHG-specific source and signed-GMRES extension is compiled from SHG/.
program shg_xc_minus
  use AI_kinds, only : DP
  use OCEAN_constants, only : Hartree2eV
  use OCEAN_mpi, only : myid, root, comm, OCEAN_mpi_init, OCEAN_mpi_finalize
  use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_copy_min, OCEAN_psi_kill, &
                        OCEAN_psi_new, OCEAN_psi_nrm, OCEAN_psi_scal
  use OCEAN_system, only : o_system, OCEAN_sys_init
  use OCEAN_timekeeper, only : OCEAN_tk_init, OCEAN_tk_printtimes, &
                               OCEAN_tk_start, OCEAN_tk_stop, tk_total
  use SHG_core_source, only : SHG_prepare_core_source, SHG_read_component
  use SHG_gmres, only : SHG_gmres_clean, SHG_gmres_do_recycle, &
                         SHG_gmres_get_parameters, SHG_gmres_setup

  implicit none

  type(o_system) :: sys
  type(OCEAN_vector) :: normalized_dc, physical_dc
  type(OCEAN_vector), allocatable :: xc_minus(:)
  real(DP), allocatable :: omega_ha(:), residual(:), vector_norm(:)
  integer, allocatable :: iterations(:)
  real(DP) :: eta_ha, preconditioner_ha, convergence, source_norm
  integer :: gmres_depth
  integer :: ierr, cleanup_ierr, mpi_ierr, i
  character(len=1) :: c_axis
  character(len=3) :: component
  character(len=128) :: filename
  logical :: mpi_started, timing_started, have_dc, have_physical_dc
  logical :: gmres_ready, have_solutions

  ierr = 0
  cleanup_ierr = 0
  mpi_ierr = 0
  mpi_started = .false.
  timing_started = .false.
  have_dc = .false.
  have_physical_dc = .false.
  gmres_ready = .false.
  have_solutions = .false.

  call OCEAN_mpi_init(ierr)
  if (ierr .ne. 0) goto 900
  mpi_started = .true.

  if (myid .eq. root) then
    write(6,'(A)') 'Standalone core-SHG first-layer solve'
    write(6,'(A)') '  X_core,c^-(omega) = [(omega+i*eta)-H_BSE^core]^{-1} d_core^c'
  endif

  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900
  call validate_core_run(sys, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_component(component, c_axis, ierr)
  if (ierr .ne. 0) goto 900

  call OCEAN_tk_init()
  call OCEAN_tk_start(tk_total)
  timing_started = .true.

  ! The private SHG solver reads the ordinary core OCEAN bse.in and allocates its
  ! own recycled-GMRES workspace while using original OCEAN_vector objects.
  call SHG_gmres_setup(sys, ierr)
  if (ierr .ne. 0) goto 900
  gmres_ready = .true.

  call SHG_gmres_get_parameters(omega_ha, eta_ha, preconditioner_ha, &
                                 convergence, gmres_depth, ierr)
  if (ierr .ne. 0) goto 900
  if (size(omega_ha) .lt. 1 .or. eta_ha .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'Invalid GMRES frequency grid or broadening.'
    ierr = 7201
    goto 900
  endif

  ! Original OCEAN constructs d_core^c with its XAS projector/matrix-element
  ! machinery, applies the core transition mask, initializes H_BSE^core, and
  ! returns a unit-norm source. The SHG wrapper retains its physical norm.
  call SHG_prepare_core_source(sys, normalized_dc, source_norm, ierr)
  if (ierr .ne. 0) goto 900
  have_dc = .true.
  if (source_norm .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'The physical core d^c source has zero norm.'
    ierr = 7202
    goto 900
  endif

  ! Solve the exact SHG resolvent.  apply_sfact=.false. keeps d^c as the RHS;
  ! solver_only=.true. suppresses the ordinary linear-absorption spectrum.
  call SHG_gmres_do_recycle(sys, normalized_dc, ierr, xc_minus, residual, &
       iterations, energy_override_ha=omega_ha, &
       resolution_override_ha=eta_ha, hamiltonian_sign_override=-1, &
       apply_sfact_override=.false., solver_only_override=.true.)
  if (ierr .ne. 0) goto 900
  have_solutions = .true.

  if (size(xc_minus) .ne. size(omega_ha) .or. &
      size(residual) .ne. size(omega_ha) .or. &
      size(iterations) .ne. size(omega_ha)) then
    if (myid .eq. root) write(6,'(A)') 'GMRES returned inconsistent result-array sizes.'
    ierr = 7203
    goto 900
  endif

  allocate(vector_norm(size(xc_minus)), stat=ierr)
  if (ierr .ne. 0) goto 900
  do i = 1, size(xc_minus)
    ! GMRES solved for the unit-normalized RHS.  Restore the physical source
    ! amplitude so the cache contains G_-(omega)|d^c>, not G_-|d^c>/||d^c||.
    call OCEAN_psi_scal(source_norm, xc_minus(i), ierr)
    if (ierr .ne. 0) goto 900
    call OCEAN_psi_nrm(vector_norm(i), xc_minus(i), ierr)
    if (ierr .ne. 0) goto 900
    vector_norm(i) = sqrt(max(0.0_DP, vector_norm(i)))
  enddo

  call OCEAN_psi_new(physical_dc, ierr)
  if (ierr .ne. 0) goto 900
  have_physical_dc = .true.
  call OCEAN_psi_copy_min(physical_dc, normalized_dc, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_scal(source_norm, physical_dc, ierr)
  if (ierr .ne. 0) goto 900
  physical_dc%kpref = normalized_dc%kpref

  ! Preserve the established OCEAN SHG cache name and complex-stream layout.
  write(filename,'(A,"_cache_d_core_",A1,".bin")') 'shg', c_axis
  call write_vector_stream(physical_dc, trim(filename), ierr)
  if (ierr .ne. 0) goto 900

  do i = 1, size(xc_minus)
    write(filename,'(A,"_cache_Gm_core_w_",A1,"_",I6.6,".bin")') 'shg', c_axis, i
    call write_vector_stream(xc_minus(i), trim(filename), ierr)
    if (ierr .ne. 0) goto 900
  enddo

  call write_manifest(sys, component, c_axis, omega_ha, eta_ha, &
       preconditioner_ha, convergence, gmres_depth, source_norm, residual, &
       iterations, vector_norm, ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A,I0,A)') 'Saved ', size(xc_minus), ' physical core X_c^-(omega) vectors.'
  endif

900 continue
  if (have_physical_dc) call OCEAN_psi_kill(physical_dc, cleanup_ierr)
  if (have_solutions) then
    do i = 1, size(xc_minus)
      call OCEAN_psi_kill(xc_minus(i), cleanup_ierr)
    enddo
  endif
  if (allocated(xc_minus)) deallocate(xc_minus)
  if (have_dc) call OCEAN_psi_kill(normalized_dc, cleanup_ierr)
  if (gmres_ready) call SHG_gmres_clean(cleanup_ierr)

  if (allocated(omega_ha)) deallocate(omega_ha)
  if (allocated(residual)) deallocate(residual)
  if (allocated(iterations)) deallocate(iterations)
  if (allocated(vector_norm)) deallocate(vector_norm)

  if (timing_started) then
    call OCEAN_tk_stop(tk_total)
    call OCEAN_tk_printtimes(myid, comm)
  endif

  if (ierr .eq. 0 .and. cleanup_ierr .ne. 0) ierr = cleanup_ierr
  if (mpi_started) call OCEAN_mpi_finalize(mpi_ierr)
  if (ierr .eq. 0 .and. mpi_ierr .ne. 0) ierr = mpi_ierr

  if (ierr .ne. 0) then
    if (mpi_started) then
      if (myid .eq. root) write(6,'(A,I0)') 'shg_xc_minus failed with ierr = ', ierr
    else
      write(6,'(A,I0)') 'shg_xc_minus failed during MPI initialization, ierr = ', ierr
    endif
    stop 1
  endif

contains

  subroutine validate_core_run(sys, ierr)
    use OCEAN_mpi, only : myid, root

    type(o_system), intent(in) :: sys
    integer, intent(inout) :: ierr

    if (.not. associated(sys%cur_run)) then
      if (myid .eq. root) write(6,'(A)') 'No run was loaded from runlist.'
      ierr = 7210
      return
    endif
    if (sys%nruns .ne. 1) then
      if (myid .eq. root) write(6,'(A)') 'This first-layer test requires exactly one XAS run.'
      ierr = 7211
      return
    endif
    if (.not. sys%cur_run%have_core .or. sys%cur_run%calc_type .ne. 'XAS') then
      if (myid .eq. root) write(6,'(A)') 'This test currently supports core/XAS SHG only.'
      ierr = 7212
      return
    endif
  end subroutine validate_core_run


  subroutine write_vector_stream(vector, path, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_INTEGER
    use OCEAN_psi, only : OCEAN_psi_min2full, OCEAN_psi_size_full, OCEAN_psi_vtor

    type(OCEAN_vector), intent(inout) :: vector
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: buffer(:)
    integer :: n, ios, local_ierr, bcast_ierr

    local_ierr = 0
    call OCEAN_psi_min2full(vector, ierr)
    if (ierr .ne. 0) return

    if (myid .eq. root) then
      n = OCEAN_psi_size_full(vector)
      allocate(buffer(n), stat=ios)
      if (ios .ne. 0) then
        local_ierr = 7220
      else
        call OCEAN_psi_vtor(vector, buffer)
        open(unit=94, file=path, form='unformatted', access='stream', &
             status='replace', iostat=ios)
        if (ios .ne. 0) then
          local_ierr = 7221
        else
          write(94, iostat=ios) buffer
          close(94)
          if (ios .ne. 0) local_ierr = 7222
        endif
        deallocate(buffer)
      endif
    endif

    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    endif
  end subroutine write_vector_stream


  subroutine write_manifest(sys, abc, axis, omega, eta, preconditioner, &
                            tolerance, depth, dnorm, res, niter, xnorm, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_INTEGER

    type(o_system), intent(in) :: sys
    character(len=3), intent(in) :: abc
    character(len=1), intent(in) :: axis
    real(DP), intent(in) :: omega(:), eta, preconditioner, tolerance, dnorm
    real(DP), intent(in) :: res(:), xnorm(:)
    integer, intent(in) :: depth, niter(:)
    integer, intent(inout) :: ierr

    character(len=128) :: path, vector_file
    integer :: j, ios, local_ierr, bcast_ierr

    local_ierr = 0
    write(path,'(A,"_xc_minus_core_",A1,".info")') 'shg', axis
    if (myid .eq. root) then
      open(unit=95, file=trim(path), form='formatted', status='replace', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7223
      else
        write(95,'(A)') 'status XC_MINUS_COMPLETE'
        write(95,'(A)') 'sector CORE_XAS'
        write(95,'(A,1X,A3)') 'component', abc
        write(95,'(A,1X,A1)') 'right_axis', axis
        write(95,'(A)') 'equation Xcore_c_minus=(omega+i*eta-H_BSE_core)^(-1)*dcore_c'
        write(95,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
        write(95,'(A,1X,ES24.16)') 'eta_hartree', eta
        write(95,'(A,1X,ES24.16)') 'preconditioner_hartree', preconditioner
        write(95,'(A,1X,ES24.16)') 'gmres_tolerance', tolerance
        write(95,'(A,1X,I0)') 'gmres_depth', depth
        write(95,'(A,1X,I0)') 'energy_count', size(omega)
        write(95,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(95,'(A,3(1X,I0))') 'core_edge_ZNL', sys%cur_run%ZNL
        write(95,'(A,1X,A2,1X,I0,1X,A2,1X,I0)') 'core_site', sys%cur_run%elname, &
             sys%cur_run%indx, sys%cur_run%corelevel, sys%cur_run%photon
        write(95,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
             sys%cur_run%num_bands, sys%nkpts, sys%nalpha
        write(95,'(A,1X,ES24.16,1X,A)') 'dipole', dnorm, &
             'shg_cache_d_core_'//axis//'.bin'
        do j = 1, size(omega)
          write(vector_file,'(A,"_cache_Gm_core_w_",A1,"_",I6.6,".bin")') 'shg', axis, j
          write(95,'(A,1X,A1,1X,I0,5(1X,ES24.16),1X,I0,1X,A)') &
               'Gminus_omega', axis, j, omega(j), omega(j)*Hartree2eV, &
               res(j), dnorm*res(j), xnorm(j), niter(j), trim(vector_file)
        enddo
        close(95)
      endif
    endif

    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    endif
  end subroutine write_manifest

end program shg_xc_minus
