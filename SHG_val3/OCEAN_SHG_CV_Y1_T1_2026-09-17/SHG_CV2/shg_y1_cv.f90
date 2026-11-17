program shg_y1_cv
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use AI_kinds, only : DP
  use OCEAN_constants, only : eV2Hartree, Hartree2eV
  use OCEAN_energies, only : OCEAN_energies_allow_full
  use OCEAN_mpi, only : myid, root, comm, OCEAN_mpi_init, OCEAN_mpi_finalize, &
                        MPI_BCAST, MPI_DOUBLE_COMPLEX, MPI_INTEGER
  use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_kill, OCEAN_psi_min2full, &
                        OCEAN_psi_new, OCEAN_psi_nrm, OCEAN_psi_rtov, &
                        OCEAN_psi_scal, OCEAN_psi_size_full, OCEAN_psi_vtor, &
                        OCEAN_psi_zero_full
  use OCEAN_system, only : o_system, OCEAN_sys_init
  use OCEAN_timekeeper, only : OCEAN_tk_init, OCEAN_tk_printtimes, &
                               OCEAN_tk_start, OCEAN_tk_stop, tk_total
  use SHG_core_source, only : SHG_prepare_core_source
  use SHG_gmres, only : SHG_gmres_clean, SHG_gmres_do_recycle, &
                         SHG_gmres_get_parameters, SHG_gmres_setup
  use SHG_rbxc_cv_manifest_io, only : SHG_rbxc_cv_manifest, &
                                    SHG_basis_match, &
                                    SHG_find_rbxc_cv_vector, &
                                    SHG_read_basis_match, &
                                    SHG_read_rbxc_cv_manifest, &
                                    SHG_validate_rbxc_cv_manifest

  implicit none

  integer, parameter :: path_len = 512

  type(o_system) :: sys
  type(SHG_rbxc_cv_manifest) :: info
  type(SHG_basis_match) :: basis
  type(OCEAN_vector) :: initialization_source, rhs
  type(OCEAN_vector), allocatable :: y1(:)
  real(DP), allocatable :: configured_energy_ha(:), target_grid_ha(:)
  real(DP), allocatable :: residual(:)
  integer, allocatable :: iterations(:)
  real(DP) :: target_energy_ev, target_energy_ha
  real(DP) :: eta_ev, eta_ha, configured_eta_ha
  real(DP) :: preconditioner_ha, convergence, initialization_source_norm
  real(DP) :: rhs_norm_squared, rhs_norm, y_norm_squared, y_norm
  real(DP) :: physical_residual, norm_difference
  integer :: gmres_depth, source_index
  integer :: ierr, cleanup_ierr, mpi_ierr
  character(len=path_len) :: input_manifest, input_vector, basis_manifest
  character(len=path_len) :: output_vector, output_manifest
  character(len=path_len) :: dipole_vector, dipole_manifest
  character(len=3) :: component
  logical :: is_shg_target
  logical :: mpi_started, timing_started, gmres_ready
  logical :: have_initialization_source, have_rhs, have_solution

  ierr = 0
  cleanup_ierr = 0
  mpi_ierr = 0
  mpi_started = .false.
  timing_started = .false.
  gmres_ready = .false.
  have_initialization_source = .false.
  have_rhs = .false.
  have_solution = .false.

  call OCEAN_mpi_init(ierr)
  if (ierr .ne. 0) goto 900
  mpi_started = .true.

  call read_command_line(input_manifest, input_vector, basis_manifest, &
                         component, target_energy_ev, output_vector, &
                         eta_ev, dipole_vector, ierr)
  if (ierr .ne. 0) goto 900

  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900

  call SHG_read_rbxc_cv_manifest(trim(input_manifest), info, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_basis_match(trim(basis_manifest), basis, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_validate_rbxc_cv_manifest(sys, info, basis, component, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_find_rbxc_cv_vector(info, trim(input_vector), source_index, ierr)
  if (ierr .ne. 0) goto 900

  target_energy_ha = target_energy_ev*eV2Hartree
  eta_ha = eta_ev*eV2Hartree
  if (.not. ieee_is_finite(target_energy_ha) .or. &
      .not. ieee_is_finite(eta_ha) .or. eta_ha .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'Target energy must be finite and eta must be positive.'
    ierr = 7701
    goto 900
  endif

  if (len_trim(output_vector) .gt. len(output_manifest)-5) then
    if (myid .eq. root) write(6,'(A)') 'Output path is too long for its companion manifest name.'
    ierr = 7705
    goto 900
  endif
  output_manifest = trim(output_vector)//'.info'
  if (len_trim(dipole_vector) .gt. len(dipole_manifest)-5) then
    if (myid .eq. root) write(6,'(A)') 'Dipole path is too long for its companion manifest name.'
    ierr = 7707
    goto 900
  endif
  dipole_manifest = trim(dipole_vector)//'.info'
  if (trim(output_vector) .eq. trim(input_vector) .or. &
      trim(dipole_vector) .eq. trim(input_vector) .or. &
      trim(dipole_vector) .eq. trim(output_vector)) then
    if (myid .eq. root) write(6,'(A)') 'Input RbXc, Y output, and dipole output must be distinct.'
    ierr = 7702
    goto 900
  endif

  is_shg_target = abs(target_energy_ha-2.0_DP*info%omega_hartree(source_index)) &
                  .le. max(1.0e-10_DP, 1.0e-8_DP*max(1.0_DP, &
                  abs(target_energy_ha), abs(2.0_DP*info%omega_hartree(source_index))))

  if (myid .eq. root) then
    write(6,'(A)') 'Standalone core-valence SHG T1 second-resolvent solve'
    write(6,'(A)') '  [(Omega+i*eta)I-H_BSE^core] Y_gv_bc = R_g<-v_b X_v_c^-(omega)'
    write(6,'(A,1X,A)') '  RbXc manifest:', trim(input_manifest)
    write(6,'(A,1X,A)') '  RbXc vector:  ', trim(input_vector)
    write(6,'(A,1X,A)') '  basis proof:  ', trim(basis_manifest)
    write(6,'(A,1X,A3)') '  component:    ', component
    write(6,'(A,F18.8,A)') '  source omega: ', info%omega_ev(source_index), ' eV'
    write(6,'(A,F18.8,A)') '  target Omega: ', target_energy_ev, ' eV'
    write(6,'(A,F18.8,A)') '  broadening:   ', eta_ev, ' eV'
    write(6,'(A,1X,A)') '  new dipole:   ', trim(dipole_vector)
    if (.not. is_shg_target) then
      write(6,'(A)') '  NOTE: requested Omega differs from 2*source omega.'
    endif
  endif

  call OCEAN_tk_init()
  call OCEAN_tk_start(tk_total)
  timing_started = .true.

  ! Read GMRES depth, preconditioner, and convergence from the staged bse.in.
  ! Its frequency list is not used; the command-line target overrides it.
  call SHG_gmres_setup(sys, ierr)
  if (ierr .ne. 0) goto 900
  gmres_ready = .true.
  call SHG_gmres_get_parameters(configured_energy_ha, configured_eta_ha, &
                                 preconditioner_ha, convergence, gmres_depth, ierr)
  if (ierr .ne. 0) goto 900

  ! OCEAN_load_data, reached through this wrapper, initializes the complete
  ! staged core BSE Hamiltonian.  The optical source itself is not the RHS of
  ! this solve; it is retained only long enough to perform initialization.
  call SHG_prepare_core_source(sys, initialization_source, &
                               initialization_source_norm, ierr)
  if (ierr .ne. 0) goto 900
  have_initialization_source = .true.

  call OCEAN_psi_new(rhs, ierr)
  if (ierr .ne. 0) goto 900
  have_rhs = .true.
  if (OCEAN_psi_size_full(rhs) .ne. info%ncon*info%nk*info%nalpha) then
    if (myid .eq. root) then
      write(6,'(A,2(1X,I0))') 'Compact CV vector/OCEAN full-size mismatch:', &
           info%ncon*info%nk*info%nalpha, OCEAN_psi_size_full(rhs)
      write(6,'(A)') 'Padding is not supported by this interface.'
    endif
    ierr = 7706
    goto 900
  endif
  call read_vector_stream(rhs, trim(input_vector), ierr)
  if (ierr .ne. 0) goto 900

  ! Enforce the same admitted core-transition space used by H_BSE^core.
  call OCEAN_energies_allow_full(sys, rhs, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_nrm(rhs_norm_squared, rhs, ierr)
  if (ierr .ne. 0) goto 900
  rhs_norm = sqrt(max(0.0_DP, rhs_norm_squared))
  if (rhs_norm .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'The admitted RbXc right-hand side has zero norm.'
    ierr = 7703
    goto 900
  endif

  norm_difference = abs(rhs_norm-info%rbxc_norm(source_index)) / &
                    max(rhs_norm, info%rbxc_norm(source_index), tiny(1.0_DP))
  if (myid .eq. root .and. norm_difference .gt. 1.0e-8_DP) then
    write(6,'(A,2(1X,ES16.8))') &
      'WARNING: read/masked RbXc norm differs from manifest:', &
      rhs_norm, info%rbxc_norm(source_index)
  endif

  ! The private GMRES implementation assumes a unit-norm arbitrary RHS.
  ! Restore this physical norm on Y after the solve.
  call OCEAN_psi_scal(1.0_DP/rhs_norm, rhs, ierr)
  if (ierr .ne. 0) goto 900
  allocate(target_grid_ha(1), stat=ierr)
  if (ierr .ne. 0) goto 900
  target_grid_ha(1) = target_energy_ha

  call SHG_gmres_do_recycle(sys, rhs, ierr, y1, residual, iterations, &
       energy_override_ha=target_grid_ha, resolution_override_ha=eta_ha, &
       hamiltonian_sign_override=-1, apply_sfact_override=.false., &
       solver_only_override=.true.)
  if (ierr .ne. 0) goto 900
  have_solution = .true.

  if (size(y1) .ne. 1 .or. size(residual) .ne. 1 .or. &
      size(iterations) .ne. 1) then
    if (myid .eq. root) write(6,'(A)') 'GMRES returned inconsistent result-array sizes.'
    ierr = 7704
    goto 900
  endif

  call OCEAN_psi_scal(rhs_norm, y1(1), ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_nrm(y_norm_squared, y1(1), ierr)
  if (ierr .ne. 0) goto 900
  y_norm = sqrt(max(0.0_DP, y_norm_squared))
  physical_residual = rhs_norm*residual(1)

  call write_vector_stream(y1(1), trim(output_vector), ierr)
  if (ierr .ne. 0) goto 900

  ! OCEAN_load_data returned the unit-normalized core dipole used only to
  ! initialize H_BSE^core.  Restore its physical norm and save it here so
  ! the final contraction never depends on an earlier core-SHG calculation.
  call OCEAN_psi_scal(initialization_source_norm, initialization_source, ierr)
  if (ierr .ne. 0) goto 900
  call write_vector_stream(initialization_source, trim(dipole_vector), ierr)
  if (ierr .ne. 0) goto 900
  call write_dipole_manifest(sys, component, initialization_source_norm, &
       trim(dipole_vector), trim(dipole_manifest), ierr)
  if (ierr .ne. 0) goto 900

  call write_output_manifest(sys, info, source_index, trim(input_manifest), &
       trim(input_vector), trim(basis_manifest), component, &
       trim(output_vector), trim(output_manifest), &
       target_energy_ha, eta_ha, configured_eta_ha, preconditioner_ha, &
       convergence, gmres_depth, rhs_norm, y_norm, residual(1), &
       physical_residual, iterations(1), is_shg_target, ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A,1X,A)') 'Saved physical Y_gv_bc vector:', trim(output_vector)
    write(6,'(A,1X,A)') 'Saved companion manifest:', trim(output_manifest)
    write(6,'(A,1X,A)') 'Saved newly generated d_g^a:', trim(dipole_vector)
    write(6,'(A,1X,A)') 'Saved dipole manifest:', trim(dipole_manifest)
  endif

900 continue
  if (have_solution) call OCEAN_psi_kill(y1(1), cleanup_ierr)
  if (allocated(y1)) deallocate(y1)
  if (have_rhs) call OCEAN_psi_kill(rhs, cleanup_ierr)
  if (have_initialization_source) then
    call OCEAN_psi_kill(initialization_source, cleanup_ierr)
  endif
  if (gmres_ready) call SHG_gmres_clean(cleanup_ierr)

  if (allocated(configured_energy_ha)) deallocate(configured_energy_ha)
  if (allocated(target_grid_ha)) deallocate(target_grid_ha)
  if (allocated(residual)) deallocate(residual)
  if (allocated(iterations)) deallocate(iterations)

  if (timing_started) then
    call OCEAN_tk_stop(tk_total)
    call OCEAN_tk_printtimes(myid, comm)
  endif

  if (ierr .eq. 0 .and. cleanup_ierr .ne. 0) ierr = cleanup_ierr
  if (mpi_started) call OCEAN_mpi_finalize(mpi_ierr)
  if (ierr .eq. 0 .and. mpi_ierr .ne. 0) ierr = mpi_ierr

  if (ierr .ne. 0) then
    if (mpi_started) then
      if (myid .eq. root) write(6,'(A,I0)') 'shg_y1_cv failed with ierr = ', ierr
    else
      write(6,'(A,I0)') 'shg_y1_cv failed during MPI initialization, ierr = ', ierr
    endif
    stop 1
  endif

contains

  subroutine read_command_line(manifest, rbxc_file, basis_file, axes, &
                               target_ev, output_file, eta_value_ev, &
                               dipole_file, ierr)
    character(len=*), intent(out) :: manifest, rbxc_file, basis_file
    character(len=*), intent(out) :: output_file, dipole_file
    character(len=3), intent(out) :: axes
    real(DP), intent(out) :: target_ev, eta_value_ev
    integer, intent(inout) :: ierr

    character(len=128) :: number
    integer :: nargs, ios

    manifest = ''
    rbxc_file = ''
    basis_file = ''
    output_file = ''
    dipole_file = ''
    axes = '???'
    target_ev = 0.0_DP
    eta_value_ev = 0.0_DP
    nargs = command_argument_count()

    if (nargs .ne. 8) then
      if (myid .eq. root) call print_usage()
      ierr = 7710
      return
    endif

    call get_path_argument(1, manifest, ierr)
    if (ierr .ne. 0) return
    call get_path_argument(2, rbxc_file, ierr)
    if (ierr .ne. 0) return
    call get_path_argument(3, basis_file, ierr)
    if (ierr .ne. 0) return
    call get_command_argument(4, value=axes)
    call lower_component(axes)
    if (.not. all_valid_axes(axes)) then
      if (myid .eq. root) write(6,'(A)') 'COMPONENT must contain exactly three axes, e.g. xyz.'
      ierr = 7714
      return
    endif
    call get_command_argument(5, value=number)
    read(number,*,iostat=ios) target_ev
    if (ios .ne. 0 .or. .not. ieee_is_finite(target_ev)) then
      if (myid .eq. root) write(6,'(A)') 'TARGET_ENERGY_EV is not a finite number.'
      ierr = 7711
      return
    endif

    call get_path_argument(6, output_file, ierr)
    if (ierr .ne. 0) return
    call get_command_argument(7, value=number)
    read(number,*,iostat=ios) eta_value_ev
    if (ios .ne. 0 .or. .not. ieee_is_finite(eta_value_ev) .or. &
        eta_value_ev .le. 0.0_DP) then
      if (myid .eq. root) write(6,'(A)') 'CORE_ETA_EV must be a positive finite number.'
      ierr = 7712
      return
    endif
    call get_path_argument(8, dipole_file, ierr)
    if (ierr .ne. 0) return
  end subroutine read_command_line


  subroutine get_path_argument(position, value, ierr)
    integer, intent(in) :: position
    character(len=*), intent(out) :: value
    integer, intent(inout) :: ierr

    integer :: arg_length

    value = ''
    call get_command_argument(position, length=arg_length)
    if (arg_length .lt. 1 .or. arg_length .gt. len(value)) then
      if (myid .eq. root) write(6,'(A,I0)') 'Missing or overlong argument ', position
      ierr = 7713
      return
    endif
    call get_command_argument(position, value=value)
  end subroutine get_path_argument


  subroutine print_usage()
    write(6,'(A)') 'Usage:'
    write(6,'(A)') '  shg_y1_cv.x CV_MANIFEST CV_RBXC_BIN BASIS_MATCH_INFO COMPONENT &'
    write(6,'(A)') '                  TARGET_ENERGY_EV Y_OUTPUT_BIN CORE_ETA_EV DIPOLE_OUTPUT_BIN'
    write(6,'(A)') 'For T1_abc, set TARGET_ENERGY_EV = 2*omega_eV of the selected RbXc dump.'
    write(6,'(A)') 'CORE_ETA_EV is explicit: the valence broadening is not reused automatically.'
  end subroutine print_usage


  subroutine lower_component(axes)
    character(len=3), intent(inout) :: axes
    integer :: j
    do j = 1, 3
      if (axes(j:j) .eq. 'X') axes(j:j) = 'x'
      if (axes(j:j) .eq. 'Y') axes(j:j) = 'y'
      if (axes(j:j) .eq. 'Z') axes(j:j) = 'z'
    enddo
  end subroutine lower_component


  pure logical function all_valid_axes(axes)
    character(len=3), intent(in) :: axes
    integer :: j
    all_valid_axes = .true.
    do j = 1, 3
      if (axes(j:j) .ne. 'x' .and. axes(j:j) .ne. 'y' .and. &
          axes(j:j) .ne. 'z') all_valid_axes = .false.
    enddo
  end function all_valid_axes


  subroutine read_vector_stream(vector, path, ierr)
    type(OCEAN_vector), intent(inout) :: vector
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: buffer(:)
    integer :: n, ios, local_ierr, bcast_ierr, nbytes_complex
    integer(8) :: actual_bytes, expected_bytes

    call OCEAN_psi_zero_full(vector, ierr)
    if (ierr .ne. 0) return
    n = OCEAN_psi_size_full(vector)
    allocate(buffer(n), stat=ierr)
    if (ierr .ne. 0) return

    local_ierr = 0
    if (myid .eq. root) then
      nbytes_complex = storage_size(cmplx(0.0_DP,0.0_DP,DP))/8
      expected_bytes = int(n,8)*int(nbytes_complex,8)
      actual_bytes = -1_8
      inquire(file=trim(path), size=actual_bytes, iostat=ios)
      if (ios .ne. 0 .or. actual_bytes .ne. expected_bytes) then
        write(6,'(A,1X,A)') 'Missing or incompatible RbXc cache:', trim(path)
        write(6,'(A,2(1X,I0))') 'Actual/expected bytes:', actual_bytes, expected_bytes
        local_ierr = 7720
      else
        open(unit=92, file=trim(path), form='unformatted', access='stream', &
             status='old', iostat=ios)
        if (ios .ne. 0) then
          local_ierr = 7721
        else
          read(92,iostat=ios) buffer
          close(92)
          if (ios .ne. 0) local_ierr = 7722
        endif
      endif
    endif

    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    else
      call MPI_BCAST(buffer, n, MPI_DOUBLE_COMPLEX, root, comm, ierr)
      if (ierr .eq. 0) call OCEAN_psi_rtov(vector, buffer)
    endif
    deallocate(buffer)
  end subroutine read_vector_stream


  subroutine write_vector_stream(vector, path, ierr)
    type(OCEAN_vector), intent(inout) :: vector
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: buffer(:)
    integer :: n, ios, local_ierr, bcast_ierr

    call OCEAN_psi_min2full(vector, ierr)
    if (ierr .ne. 0) return

    local_ierr = 0
    if (myid .eq. root) then
      n = OCEAN_psi_size_full(vector)
      allocate(buffer(n), stat=ios)
      if (ios .ne. 0) then
        local_ierr = 7730
      else
        call OCEAN_psi_vtor(vector, buffer)
        open(unit=93, file=trim(path), form='unformatted', access='stream', &
             status='new', iostat=ios)
        if (ios .ne. 0) then
          local_ierr = 7731
        else
          write(93,iostat=ios) buffer
          close(93)
          if (ios .ne. 0) local_ierr = 7732
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


  subroutine write_dipole_manifest(sys, axes, dnorm, vector_path, &
                                   manifest_path, ierr)
    type(o_system), intent(in) :: sys
    character(len=3), intent(in) :: axes
    real(DP), intent(in) :: dnorm
    character(len=*), intent(in) :: vector_path, manifest_path
    integer, intent(inout) :: ierr
    integer :: ios, local_ierr, bcast_ierr

    local_ierr = 0
    if (myid .eq. root) then
      open(unit=95, file=trim(manifest_path), form='formatted', &
           status='new', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7735
      else
        write(95,'(A)') 'status XC_MINUS_COMPLETE'
        write(95,'(A)') 'sector CORE_XAS'
        write(95,'(A)') 'purpose final_core_dipole_for_T1_gv'
        write(95,'(A,1X,A3)') 'component', axes
        write(95,'(A,1X,A1)') 'right_axis', axes(1:1)
        write(95,'(A)') 'equation d_g^a=physical_core_absorption_source'
        write(95,'(A)') 'origin generated_in_same_fresh_core_stage_as_Y_gv'
        write(95,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
        write(95,'(A)') 'units bohr'
        write(95,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(95,'(A,3(1X,I0))') 'core_edge_ZNL', sys%cur_run%ZNL
        write(95,'(A,1X,A2,1X,I0,1X,A2,1X,I0)') 'core_site', &
             sys%cur_run%elname, sys%cur_run%indx, &
             sys%cur_run%corelevel, sys%cur_run%photon
        write(95,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
             sys%cur_run%num_bands, sys%nkpts, sys%nalpha
        write(95,'(A,1X,ES24.16,1X,A)') 'dipole', dnorm, trim(vector_path)
        close(95)
      endif
    endif
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    endif
  end subroutine write_dipole_manifest


  subroutine write_output_manifest(sys, input, index, input_manifest_path, &
       input_vector_path, basis_manifest_path, axes, output_vector_path, &
       manifest_path, target_ha, eta, &
       configured_eta, preconditioner, tolerance, depth, rbxc_norm, &
       output_norm, normalized_residual, physical_residual_value, niter, &
       target_is_two_omega, ierr)
    type(o_system), intent(in) :: sys
    type(SHG_rbxc_cv_manifest), intent(in) :: input
    integer, intent(in) :: index, depth, niter
    character(len=*), intent(in) :: input_manifest_path, input_vector_path
    character(len=*), intent(in) :: basis_manifest_path
    character(len=3), intent(in) :: axes
    character(len=*), intent(in) :: output_vector_path, manifest_path
    real(DP), intent(in) :: target_ha, eta, configured_eta, preconditioner
    real(DP), intent(in) :: tolerance, rbxc_norm, output_norm
    real(DP), intent(in) :: normalized_residual, physical_residual_value
    logical, intent(in) :: target_is_two_omega
    integer, intent(inout) :: ierr

    integer :: ios, local_ierr, bcast_ierr

    local_ierr = 0
    if (myid .eq. root) then
      open(unit=94, file=trim(manifest_path), form='formatted', &
           status='new', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7740
      else
        write(94,'(A)') 'status Y1_GV_COMPLETE'
        write(94,'(A)') 'sector CORE_FROM_VALENCE'
        write(94,'(A)') 'pathway GV'
        write(94,'(A,1X,A3)') 'component', axes
        write(94,'(A,1X,A1)') 'left_axis', axes(1:1)
        write(94,'(A,1X,A1)') 'middle_axis', input%middle_axis
        write(94,'(A,1X,A1)') 'right_axis', input%right_axis
        write(94,'(A)') 'equation Y_gv_bc=Gminus_core(Omega)*R_gv_b*X_v_c'
        write(94,'(A)') 'linear_system [(Omega+i*eta_core)I-H_BSE_core]Y_gv_bc=RbXc_cv'
        write(94,'(A)') 't1_contraction T1_gv_abc=<d_core^a|Y_gv_bc_at_Omega=2omega>'
        write(94,'(A)') 'hamiltonian_sign -1'
        write(94,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
        write(94,'(A,1X,A)') 'input_rbxc_operator', trim(input%operator_label)
        write(94,'(A,1X,A)') 'input_rbxc_binary_layout', trim(input%binary_layout)
        write(94,'(A,1X,A)') 'input_manifest', trim(input_manifest_path)
        write(94,'(A,1X,A)') 'input_vector', trim(input_vector_path)
        write(94,'(A,1X,A)') 'basis_match_manifest', trim(basis_manifest_path)
        write(94,'(A,1X,A)') 'output_vector', trim(output_vector_path)
        write(94,'(A,1X,A)') 'paired_dipole_vector', trim(dipole_vector)
        write(94,'(A,1X,A)') 'paired_dipole_manifest', trim(dipole_manifest)
        write(94,'(A,1X,I0)') 'energy_count', 1
        write(94,'(A,1X,I0)') 'source_index', index
        write(94,'(A,1X,ES24.16)') 'source_omega_hartree', input%omega_hartree(index)
        write(94,'(A,1X,ES24.16)') 'source_omega_ev', input%omega_ev(index)
        write(94,'(A,1X,ES24.16)') 'target_energy_hartree', target_ha
        write(94,'(A,1X,ES24.16)') 'target_energy_ev', target_ha*Hartree2eV
        write(94,'(A,1X,L1)') 'target_is_two_omega', target_is_two_omega
        write(94,'(A,1X,ES24.16)') 'eta_hartree', eta
        write(94,'(A,1X,ES24.16)') 'eta_ev', eta*Hartree2eV
        write(94,'(A,1X,ES24.16)') 'input_valence_eta_hartree', &
             input%eta_valence_hartree
        write(94,'(A,3(1X,ES24.16))') 'input_q_reduced', input%q_reduced
        write(94,'(A,1X,ES24.16)') 'bse_in_eta_hartree', configured_eta
        write(94,'(A,1X,ES24.16)') 'preconditioner_hartree', preconditioner
        write(94,'(A,1X,ES24.16)') 'gmres_tolerance_squared_norm', tolerance
        write(94,'(A,1X,I0)') 'gmres_depth', depth
        write(94,'(A,1X,I0)') 'gmres_iterations', niter
        write(94,'(A,1X,ES24.16)') 'normalized_residual_norm', normalized_residual
        write(94,'(A,1X,ES24.16)') 'physical_residual_norm', physical_residual_value
        write(94,'(A,1X,ES24.16)') 'rbxc_norm_manifest', input%rbxc_norm(index)
        write(94,'(A,1X,ES24.16)') 'rbxc_norm_read_after_mask', rbxc_norm
        write(94,'(A,1X,ES24.16)') 'y1_norm', output_norm
        write(94,'(A,3(1X,I0))') 'core_edge_ZNL', sys%cur_run%ZNL
        write(94,'(A,1X,A,1X,I0,1X,A,1X,I0)') 'core_site', &
             trim(sys%cur_run%elname), sys%cur_run%indx, &
             trim(sys%cur_run%corelevel), sys%cur_run%photon
        write(94,'(A,4(1X,I0))') 'input_brange', input%brange
        write(94,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(94,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
             sys%cur_run%num_bands, sys%nkpts, sys%nalpha
        write(94,'(A,1X,A1,1X,A1,1X,I0,7(1X,ES24.16),1X,I0,1X,A)') &
             'Y1_GV', input%middle_axis, input%right_axis, index, &
             input%omega_hartree(index), input%omega_ev(index), target_ha, &
             target_ha*Hartree2eV, eta, rbxc_norm, output_norm, niter, &
             trim(output_vector_path)
        close(94)
      endif
    endif

    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    endif
  end subroutine write_output_manifest

end program shg_y1_cv
