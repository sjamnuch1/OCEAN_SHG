! Copyright (C) 2026 OCEAN collaboration
!
! Standalone valence first-layer SHG test:
!   |X_valence,c^-(omega)> = [(omega + i*eta) - H_BSE^valence]^{-1} |d_valence^c>
!
! This program uses only the public machinery in an unmodified OCEAN checkout.
! Every SHG-specific source and signed-GMRES extension is compiled from SHG_VAL/.
program shg_xc_minus_valence
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use AI_kinds, only : DP
  use OCEAN_constants, only : Hartree2eV
  use OCEAN_mpi, only : myid, root, comm, OCEAN_mpi_init, OCEAN_mpi_finalize
  use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_copy_min, OCEAN_psi_kill, &
                        OCEAN_psi_new, OCEAN_psi_nrm, OCEAN_psi_scal
  use OCEAN_system, only : o_system, OCEAN_sys_init
  use OCEAN_timekeeper, only : OCEAN_tk_init, OCEAN_tk_printtimes, &
                               OCEAN_tk_start, OCEAN_tk_stop, tk_total
  use SHG_valence_source, only : SHG_prepare_valence_source, SHG_read_component
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
  integer :: ierr, cleanup_ierr, release_ierr, mpi_ierr, i, solutions_released
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
  solutions_released = 0

  call OCEAN_mpi_init(ierr)
  if (ierr .ne. 0) goto 900
  mpi_started = .true.

  if (myid .eq. root) then
    write(6,'(A)') 'Standalone valence-SHG first-layer solve'
    write(6,'(A)') '  X_valence,c^-(omega) = [(omega+i*eta)-H_BSE^valence]^{-1} d_valence^c'
  endif

  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900
  call validate_valence_run(sys, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_component(component, c_axis, ierr)
  if (ierr .ne. 0) goto 900

  call OCEAN_tk_init()
  call OCEAN_tk_start(tk_total)
  timing_started = .true.

  ! The private SHG solver reads the ordinary valence OCEAN bse.in and allocates its
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

  ! Load the finite-q VAL source, mask it, initialize H_V, restore the
  ! length-gauge phase, and retain the physical source norm.
  call SHG_prepare_valence_source(sys, c_axis, normalized_dc, source_norm, ierr)
  if (ierr .ne. 0) goto 900
  have_dc = .true.
  if (source_norm .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'The physical valence d^c source has zero norm.'
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

  if (any(.not. ieee_is_finite(residual)) .or. any(residual**2 > convergence)) then
    if (myid .eq. root) write(6,'(A)') 'GMRES residual failed the configured squared-norm tolerance.'
    ierr = 7204
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

  ! The converged responses are now independent OCEAN_vector objects.  The
  ! Arnoldi/recycling basis is no longer needed, and retaining it while
  ! OCEAN_psi_min2full materializes vectors for output creates the largest
  ! memory peak of the run.  Release the 2*gmres_depth basis vectors first.
  call SHG_gmres_clean(ierr)
  if (ierr .ne. 0) goto 900
  gmres_ready = .false.

  call OCEAN_psi_new(physical_dc, ierr)
  if (ierr .ne. 0) goto 900
  have_physical_dc = .true.
  call OCEAN_psi_copy_min(physical_dc, normalized_dc, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_scal(source_norm, physical_dc, ierr)
  if (ierr .ne. 0) goto 900
  physical_dc%kpref = normalized_dc%kpref

  ! Preserve the established OCEAN SHG cache name and complex-stream layout.
  write(filename,'(A,"_cache_d_valence_",A1,".bin")') 'shg', c_axis
  call write_vector_stream(sys, physical_dc, trim(filename), ierr)
  if (ierr .ne. 0) goto 900

  ! min2full may attach a full, replicated representation to the vector.
  ! The dipole cache is complete, so release both its minimal and full
  ! representations before materializing the first exciton vector.
  release_ierr = 0
  call OCEAN_psi_kill(physical_dc, release_ierr)
  if (release_ierr .ne. 0) then
    ierr = release_ierr
    goto 900
  endif
  have_physical_dc = .false.

  do i = 1, size(xc_minus)
    write(filename,'(A,"_cache_Gm_valence_w_",A1,"_",I6.6,".bin")') 'shg', c_axis, i
    call write_vector_stream(sys, xc_minus(i), trim(filename), ierr)
    if (ierr .ne. 0) goto 900

    ! Do not retain the full representation produced for output.  At most one
    ! materialized response now coexists with the still-minimal responses.
    release_ierr = 0
    call OCEAN_psi_kill(xc_minus(i), release_ierr)
    if (release_ierr .ne. 0) then
      ierr = release_ierr
      goto 900
    endif
    solutions_released = i
  enddo

  call write_manifest(sys, component, c_axis, omega_ha, eta_ha, &
       preconditioner_ha, convergence, gmres_depth, source_norm, residual, &
       iterations, vector_norm, ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A,I0,A)') 'Saved ', size(xc_minus), ' physical valence X_c^-(omega) vectors.'
  endif

900 continue
  if (have_physical_dc) call OCEAN_psi_kill(physical_dc, cleanup_ierr)
  if (have_solutions) then
    do i = solutions_released + 1, size(xc_minus)
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
      if (myid .eq. root) write(6,'(A,I0)') 'shg_xc_minus_valence failed with ierr = ', ierr
    else
      write(6,'(A,I0)') 'shg_xc_minus_valence failed during MPI initialization, ierr = ', ierr
    endif
    stop 1
  endif

contains

  subroutine validate_valence_run(sys, ierr)
    use OCEAN_mpi, only : myid, root

    type(o_system), intent(in) :: sys
    integer, intent(inout) :: ierr

    if (.not. associated(sys%cur_run)) then
      if (myid .eq. root) write(6,'(A)') 'No run was loaded from runlist.'
      ierr = 7210
      return
    endif
    if (sys%nruns .ne. 1) then
      if (myid .eq. root) write(6,'(A)') 'This first-layer test requires exactly one VAL run.'
      ierr = 7211
      return
    endif
    if (.not. sys%cur_run%have_val .or. sys%cur_run%have_core .or. &
        sys%cur_run%calc_type .ne. 'VAL') then
      if (myid .eq. root) write(6,'(A)') 'Requires a pure VAL run, not a RIXS source.'
      ierr = 7212
      return
    endif
    if (sys%tmel_selector .ne. 1 .or. sys%enk_selector .ne. 0 .or. sys%bloch_selector .ne. 3) then
      if (myid .eq. root) write(6,'(A)') 'Use selectors tmel=1, enk=0, bloch=3 from prepare_valence.pl.'
      ierr = 7213
      return
    endif
    if (sys%cur_run%bwflg .or. sys%cur_run%backf .or. sys%nbw .ne. 1) then
      if (myid .eq. root) write(6,'(A)') 'Only the resonant/TDA valence basis is supported.'
      ierr = 7214
      return
    endif
  end subroutine validate_valence_run


  ! Stream only physical entries: never OCEAN padding.
  subroutine write_vector_stream(sys, vector, path, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_INTEGER
    use OCEAN_psi, only : OCEAN_psi_min2full
    type(o_system), intent(in) :: sys
    type(OCEAN_vector), intent(inout) :: vector
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr
    integer :: ik, ibeta, ios, close_ios, local_ierr, bcast_ierr
    local_ierr = 0
    call OCEAN_psi_min2full(vector, ierr)
    if (ierr .ne. 0) return
    if (myid .eq. root) then
      open(unit=94, file=path, form='unformatted', access='stream', status='replace', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7221
      else
        do ibeta = 1, sys%nbeta
          do ik = 1, sys%nkpts
            write(94, iostat=ios) cmplx( &
                vector%valr(1:sys%cur_run%num_bands,1:sys%cur_run%val_bands,ik,ibeta,1), &
                vector%vali(1:sys%cur_run%num_bands,1:sys%cur_run%val_bands,ik,ibeta,1), DP)
            if (ios .ne. 0) then
              local_ierr = 7222
              exit
            endif
          enddo
          if (local_ierr .ne. 0) exit
        enddo
        close(94, iostat=close_ios)
        if (close_ios .ne. 0) local_ierr = 7222
      endif
    endif
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    else
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
    write(path,'(A,"_xc_minus_valence_",A1,".info")') 'shg', axis
    if (myid .eq. root) then
      open(unit=95, file=trim(path), form='formatted', status='replace', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7223
      else
        write(95,'(A)') 'status XC_MINUS_COMPLETE'
        write(95,'(A)') 'sector VALENCE_TDA'
        write(95,'(A,1X,A3)') 'component', abc
        write(95,'(A,1X,A1)') 'right_axis', axis
        write(95,'(A)') 'equation Xvalence_c_minus=(omega+i*eta-H_BSE_valence)^(-1)*dvalence_c'
        write(95,'(A)') 'binary_layout raw_complex_DP_stream_Fortran_ncon_nval_nk_nbeta_nbw_no_padding'
        write(95,'(A,1X,ES24.16)') 'eta_hartree', eta
        write(95,'(A,1X,ES24.16)') 'preconditioner_hartree', preconditioner
        write(95,'(A,1X,ES24.16)') 'gmres_tolerance', tolerance
        write(95,'(A,1X,I0)') 'gmres_depth', depth
        write(95,'(A,1X,I0)') 'energy_count', size(omega)
        write(95,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(95,'(A,5(1X,I0))') 'dimensions ncon nval nk nbeta nbw', &
             sys%cur_run%num_bands, sys%cur_run%val_bands, sys%nkpts, sys%nbeta, sys%nbw
        write(95,'(A,4(1X,I0))') 'brange', sys%brange
        write(95,'(A,4(1X,I0))') 'file_brange', sys%file_brange
        write(95,'(A,3(1X,ES24.16))') 'q_reduced', sys%qinunitsofbvectors
        write(95,'(A,3(1X,ES24.16))') 'q_cart_bohr_inverse', matmul(sys%bvec,sys%qinunitsofbvectors)
        write(95,'(A)') 'source_convention d_code=-i*sign(q_axis)*conjg(M_vn(q))/abs(q)'
        write(95,'(A)') 'nonparallel_q_means d_code=sign(q_axis)*(qhat_dot_d), not pure d_axis'
        write(95,'(A)') 'source_units bohr'
        write(95,'(A)') 'solution_units bohr_per_hartree'
        write(95,'(A)') 'energy_reference valence_transition_energy_no_core_edge_shift'
        write(95,'(A)') 'Gminus_columns axis index omega_Ha omega_eV residual_normalized residual_physical norm iterations file'
        write(95,'(A,1X,ES24.16,1X,A)') 'dipole', dnorm, &
             'shg_cache_d_valence_'//axis//'.bin'
        do j = 1, size(omega)
          write(vector_file,'(A,"_cache_Gm_valence_w_",A1,"_",I6.6,".bin")') 'shg', axis, j
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

end program shg_xc_minus_valence

