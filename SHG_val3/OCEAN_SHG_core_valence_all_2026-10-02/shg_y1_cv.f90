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
  use SHG_core_source, only : SHG_prepare_core_source, SHG_read_component
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
  real(DP), allocatable :: configured_energy_ha(:), residual(:)
  integer, allocatable :: iterations(:)
  real(DP) :: target_grid_ha(1), target_energy_ha, eta_ha
  real(DP) :: preconditioner_ha, convergence, initialization_source_norm
  real(DP) :: rhs_norm_squared, rhs_norm, y_norm_squared, y_norm
  real(DP) :: physical_residual, norm_difference
  integer :: gmres_depth, source_index, start_index, end_index
  integer :: ierr, cleanup_ierr, mpi_ierr
  character(len=path_len) :: input_manifest, input_vector, basis_manifest
  character(len=path_len) :: output_prefix, output_vector, output_manifest
  character(len=path_len) :: dipole_vector, dipole_manifest, start_file, batch_manifest
  character(len=3) :: component, staged_component
  character(len=1) :: staged_right_axis
  logical :: is_shg_target, show_help
  logical :: mpi_started, timing_started, gmres_ready
  logical :: have_initialization_source, have_rhs, solution_ready

  ierr = 0
  cleanup_ierr = 0
  mpi_ierr = 0
  mpi_started = .false.
  timing_started = .false.
  gmres_ready = .false.
  have_initialization_source = .false.
  have_rhs = .false.
  solution_ready = .false.

  call OCEAN_mpi_init(ierr)
  if (ierr .ne. 0) goto 900
  mpi_started = .true.
  call read_command_line(ierr)
  if (ierr .ne. 0 .or. show_help) goto 900
  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_rbxc_cv_manifest(trim(input_manifest), info, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_basis_match(trim(basis_manifest), basis, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_validate_rbxc_cv_manifest(sys, info, basis, component, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_read_component(staged_component, staged_right_axis, ierr)
  if (ierr .ne. 0) goto 900
  if (component .ne. staged_component) then
    if (myid .eq. root) write(6,'(A)') '--component differs from the staged shg_component; restage the dipole.'
    ierr = 7707
    goto 900
  endif
  if (len_trim(start_file) .gt. 0) then
    call SHG_find_rbxc_cv_vector(info, trim(start_file), start_index, ierr)
    if (ierr .ne. 0) goto 900
  endif
  if (end_index .eq. 0) end_index = info%energy_count
  if (start_index .lt. 1 .or. end_index .gt. info%energy_count .or. &
      start_index .gt. end_index) then
    if (myid .eq. root) write(6,'(A)') 'Invalid start/end index for the full CV manifest.'
    ierr = 7701
    goto 900
  endif

  ! Only the second-resolvent parameters and targets come from bse2.in.
  ! bse.in remains present for unchanged OCEAN core-multiplet header reads.
  call SHG_gmres_setup(sys, ierr, input_filename='bse2.in')
  if (ierr .ne. 0) goto 900
  gmres_ready = .true.
  call SHG_gmres_get_parameters(configured_energy_ha, eta_ha, &
                              preconditioner_ha, convergence, gmres_depth, ierr)
  if (ierr .ne. 0) goto 900
  if (size(configured_energy_ha) .ne. info%energy_count) then
    if (myid .eq. root) then
      write(6,'(A,2(1X,I0))') 'bse2.in/manifest full-grid count mismatch:', &
                             size(configured_energy_ha), info%energy_count
      write(6,'(A)') 'Keep the complete bse2.in grid, including entries before --start-index.'
    endif
    ierr = 7702
    goto 900
  endif
  call preflight_batch(ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A)') 'Batch core-valence SHG second-resolvent solve; GMRES input: bse2.in'
    write(6,'(A)') '  [(Omega+i*eta)I-H_BSE^core] Y_gv_bc = R_g<-v_b X_v_c^-(omega)'
    write(6,'(A,1X,A)') '  CV manifest:', trim(input_manifest)
    write(6,'(A,1X,A3)') '  component:', component
    write(6,'(A,2(1X,I0))') '  manifest indices (inclusive):', start_index, end_index
    write(6,'(A,F18.8,A)') '  core broadening:', eta_ha*Hartree2eV, ' eV'
    write(6,'(A)') '  Grid entry i is paired with manifest entry i; indices are never renumbered.'
  endif
  call OCEAN_tk_init()
  call OCEAN_tk_start(tk_total)
  timing_started = .true.

  ! Expensive wavefunction/Hamiltonian initialization happens once per batch.
  call SHG_prepare_core_source(sys, initialization_source, initialization_source_norm, ierr)
  if (ierr .ne. 0) goto 900
  have_initialization_source = .true.
  call OCEAN_psi_scal(initialization_source_norm, initialization_source, ierr)
  if (ierr .ne. 0) goto 900
  call write_vector_stream(initialization_source, trim(dipole_vector), ierr)
  if (ierr .ne. 0) goto 900
  call write_dipole_manifest(sys, component, initialization_source_norm, &
                            trim(dipole_vector), trim(dipole_manifest), ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_kill(initialization_source, ierr)
  if (ierr .ne. 0) goto 900
  have_initialization_source = .false.

  call OCEAN_psi_new(rhs, ierr)
  if (ierr .ne. 0) goto 900
  have_rhs = .true.
  if (int(OCEAN_psi_size_full(rhs),8) .ne. &
      int(info%ncon,8)*int(info%nk,8)*int(info%nalpha,8)) then
    if (myid .eq. root) write(6,'(A)') 'Compact CV vector/OCEAN full-size mismatch; padding unsupported.'
    ierr = 7706
    goto 900
  endif
  call begin_batch_manifest(ierr)
  if (ierr .ne. 0) goto 900

  do source_index = start_index, end_index
    target_energy_ha = configured_energy_ha(source_index)
    target_grid_ha(1) = target_energy_ha
    is_shg_target = two_omega_match(target_energy_ha, info%omega_hartree(source_index))
    call entry_paths(source_index, input_vector, output_vector, output_manifest)
    if (myid .eq. root) then
      write(6,'(A,I0,A,F18.8,A,F18.8,A)') 'Entry ',source_index, &
        ': omega=',info%omega_ev(source_index),' eV; Omega=',target_energy_ha*Hartree2eV,' eV'
    endif
    call read_vector_stream(rhs, trim(input_vector), ierr)
    if (ierr .ne. 0) goto 900
    call OCEAN_energies_allow_full(sys, rhs, ierr)
    if (ierr .ne. 0) goto 900
    call OCEAN_psi_nrm(rhs_norm_squared, rhs, ierr)
    if (ierr .ne. 0) goto 900
    if (.not. ieee_is_finite(rhs_norm_squared) .or. rhs_norm_squared .le. tiny(1.0_DP)) then
      if (myid .eq. root) write(6,'(A)') 'The admitted RbXc RHS has zero or nonfinite norm.'
      ierr = 7703
      goto 900
    endif
    rhs_norm = sqrt(rhs_norm_squared)
    norm_difference = abs(rhs_norm-info%rbxc_norm(source_index)) / &
                      max(rhs_norm, info%rbxc_norm(source_index), tiny(1.0_DP))
    if (myid .eq. root .and. norm_difference .gt. 1.0e-8_DP) then
      write(6,'(A,2(1X,ES16.8))') 'WARNING: read/masked RHS norm differs from manifest:', &
                                 rhs_norm, info%rbxc_norm(source_index)
    endif
    call OCEAN_psi_scal(1.0_DP/rhs_norm, rhs, ierr)
    if (ierr .ne. 0) goto 900

    ! A one-element override is deliberate: this energy has its OWN RHS.
    ! GMRES's local history starts fresh on each call; H and storage stay loaded.
    call SHG_gmres_do_recycle(sys, rhs, ierr, y1, residual, iterations, &
         energy_override_ha=target_grid_ha, resolution_override_ha=eta_ha, &
         hamiltonian_sign_override=-1, apply_sfact_override=.false., solver_only_override=.true.)
    if (ierr .ne. 0) goto 900
    solution_ready = .true.
    if (size(y1) .ne. 1 .or. size(residual) .ne. 1 .or. size(iterations) .ne. 1) then
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
    call write_output_manifest(sys, info, source_index, trim(input_manifest), &
         trim(input_vector), trim(basis_manifest), component, trim(output_vector), &
         trim(output_manifest), target_energy_ha, eta_ha, eta_ha, preconditioner_ha, &
         convergence, gmres_depth, rhs_norm, y_norm, residual(1), physical_residual, &
         iterations(1), is_shg_target, ierr)
    if (ierr .ne. 0) goto 900
    call append_batch_entry(ierr)
    if (ierr .ne. 0) goto 900
    call release_solution(ierr)
    if (ierr .ne. 0) goto 900
  enddo
  call complete_batch_manifest(ierr)
  if (ierr .ne. 0) goto 900
  if (myid .eq. root) write(6,'(A,1X,A)') 'Batch complete:',trim(batch_manifest)

900 continue
  call release_solution(cleanup_ierr)
  if (have_rhs) call OCEAN_psi_kill(rhs, cleanup_ierr)
  if (have_initialization_source) call OCEAN_psi_kill(initialization_source, cleanup_ierr)
  if (gmres_ready) call SHG_gmres_clean(cleanup_ierr)
  if (allocated(configured_energy_ha)) deallocate(configured_energy_ha)
  if (timing_started) then
    call OCEAN_tk_stop(tk_total)
    call OCEAN_tk_printtimes(myid, comm)
  endif
  if (ierr .eq. 0 .and. cleanup_ierr .ne. 0) ierr = cleanup_ierr
  if (mpi_started) call OCEAN_mpi_finalize(mpi_ierr)
  if (ierr .eq. 0 .and. mpi_ierr .ne. 0) ierr = mpi_ierr
  if (ierr .ne. 0) then
    if (myid .eq. root) write(6,'(A,I0)') 'shg_y1_cv failed with ierr = ',ierr
    stop 1
  endif

contains

  subroutine read_command_line(ierr)
    integer, intent(inout) :: ierr
    character(len=path_len) :: option, value
    integer :: pos, nargs, ios
    logical :: explicit_start
    input_manifest = 'shg_rbxc_cv.info'
    basis_manifest = 'cv_core_basis_match.info'
    component = '???'
    output_prefix = ''
    dipole_vector = ''
    start_file = ''
    start_index = 1
    end_index = 0
    show_help = .false.
    explicit_start = .false.
    pos = 1
    nargs = command_argument_count()
    do while (pos .le. nargs)
      call get_path_argument(pos, option, ierr)
      if (ierr .ne. 0) return
      if (trim(option) .eq. '--help' .or. trim(option) .eq. '-h') then
        show_help = .true.
        if (myid .eq. root) call print_usage()
        return
      endif
      if (pos .eq. nargs) then
        ierr = 7710
        exit
      endif
      call get_path_argument(pos+1, value, ierr)
      if (ierr .ne. 0) return
      select case(trim(option))
      case('--manifest'); input_manifest = value
      case('--basis'); basis_manifest = value
      case('--component')
        if (len_trim(value) .ne. 3) then
          ierr = 7714
          exit
        endif
        component = value
        call lower_component(component)
      case('--output-prefix'); output_prefix = value
      case('--dipole-output'); dipole_vector = value
      case('--start-file'); start_file = value
      case('--start-index')
        read(value,*,iostat=ios) start_index
        if (ios .ne. 0 .or. start_index .lt. 1) then
          ierr = 7711
          exit
        endif
        explicit_start = .true.
      case('--end-index')
        read(value,*,iostat=ios) end_index
        if (ios .ne. 0 .or. end_index .lt. 1) then
          ierr = 7711
          exit
        endif
      case default
        ierr = 7710
        exit
      end select
      pos = pos+2
    enddo
    if (explicit_start .and. len_trim(start_file) .gt. 0) ierr = 7711
    if (.not. all_valid_axes(component)) ierr = 7714
    if (ierr .ne. 0) then
      if (myid .eq. root) call print_usage()
      return
    endif
    if (len_trim(output_prefix) .eq. 0) output_prefix = 'shg_cache_Y1_cv_'//component(2:3)
    if (len_trim(dipole_vector) .eq. 0) dipole_vector = 'shg_cache_d_core_'//component(1:1)//'.bin'
    if (len_trim(output_prefix) .gt. path_len-32 .or. len_trim(dipole_vector) .gt. path_len-5) then
      ierr = 7713
      return
    endif
    dipole_manifest = trim(dipole_vector)//'.info'
    batch_manifest = trim(output_prefix)//'.batch.info'
  end subroutine read_command_line

  subroutine print_usage()
    write(6,'(A)') 'Usage: shg_y1_cv.x --component ABC [options]'
    write(6,'(A)') '  --manifest FILE      default shg_rbxc_cv.info'
    write(6,'(A)') '  --basis FILE         default cv_core_basis_match.info'
    write(6,'(A)') '  --output-prefix PATH default shg_cache_Y1_cv_BC'
    write(6,'(A)') '  --dipole-output FILE default shg_cache_d_core_A.bin'
    write(6,'(A)') '  --start-index N      default 1; retain ORIGINAL grid and manifest indexing'
    write(6,'(A)') '  --start-file FILE    alternative: begin at this manifest-listed RbXc file'
    write(6,'(A)') '  --end-index N        default final manifest index (inclusive)'
    write(6,'(A)') 'bse2.in must contain the FULL list/loop target grid, one target per manifest entry.'
    write(6,'(A)') 'Existing selected Y/dipole outputs are replaced. Unselected Y files are retained.'
  end subroutine print_usage

  subroutine release_solution(ierr)
    integer, intent(inout) :: ierr
    integer :: k, local_ierr
    if (allocated(y1)) then
      if (solution_ready) then
        do k=1,size(y1)
          local_ierr = 0
          call OCEAN_psi_kill(y1(k),local_ierr)
          if (ierr .eq. 0 .and. local_ierr .ne. 0) ierr = local_ierr
        enddo
      endif
      deallocate(y1)
    endif
    solution_ready = .false.
    if (allocated(residual)) deallocate(residual)
    if (allocated(iterations)) deallocate(iterations)
  end subroutine release_solution

  logical function two_omega_match(target,omega)
    real(DP), intent(in) :: target, omega
    two_omega_match = abs(target-2.0_DP*omega) .le. &
      max(1.0e-10_DP,1.0e-8_DP*max(1.0_DP,abs(target),abs(2.0_DP*omega)))
  end function two_omega_match

  subroutine entry_paths(index,rbxc,outbin,outinfo)
    integer, intent(in) :: index
    character(len=*), intent(out) :: rbxc,outbin,outinfo
    character(len=32) :: number
    integer :: slash
    write(number,'(I0)') index
    number = repeat('0',max(0,6-len_trim(number)))//trim(number)
    rbxc = info%rbxc_file(index)
    ! Relative RHS paths are resolved relative to their manifest.
    slash = scan(trim(input_manifest),'/',back=.true.)
    if (rbxc(1:1) .ne. '/' .and. slash .gt. 0) rbxc = input_manifest(:slash)//trim(rbxc)
    outbin = trim(output_prefix)//'_'//trim(number)//'.bin'
    outinfo = trim(outbin)//'.info'
  end subroutine entry_paths

  subroutine preflight_batch(ierr)
    integer, intent(inout) :: ierr
    integer :: k, ios, local_ierr, bcast_ierr
    integer(8) :: nbytes, expected
    character(len=path_len) :: r,o,m
    local_ierr = 0
    if (myid .eq. root) then
      expected = int(info%ncon,8)*int(info%nk,8)*int(info%nalpha,8)* &
                 int(storage_size(cmplx(0.0_DP,0.0_DP,DP))/8,8)
      do k=1,info%energy_count
        if (.not. ieee_is_finite(info%omega_hartree(k)) .or. &
            .not. ieee_is_finite(info%rbxc_norm(k))) then
          local_ierr = 7750
          exit
        endif
        if (k .lt. start_index .or. k .gt. end_index) cycle
        call entry_paths(k,r,o,m)
        nbytes = -1_8
        inquire(file=trim(r),size=nbytes,iostat=ios)
        if (ios .ne. 0 .or. nbytes .ne. expected) then
          write(6,'(A,I0,1X,A)') 'Missing/incompatible selected RbXc entry ',k,trim(r)
          local_ierr = 7751
          exit
        endif
        if (trim(r) .eq. trim(o) .or. trim(r) .eq. trim(dipole_vector) .or. &
            trim(o) .eq. trim(dipole_vector)) then
          local_ierr = 7752
          exit
        endif
        if (.not. two_omega_match(configured_energy_ha(k),info%omega_hartree(k))) then
          write(6,'(A,I0,A)') 'WARNING: bse2.in target ',k, &
            ' differs from twice its source omega; target_is_two_omega will be false.'
        endif
      enddo
    endif
    call MPI_BCAST(local_ierr,1,MPI_INTEGER,root,comm,bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    else
      ierr = local_ierr
    endif
  end subroutine preflight_batch

  subroutine begin_batch_manifest(ierr)
    integer, intent(inout) :: ierr
    integer :: ios, local_ierr, bcast_ierr
    local_ierr = 0
    if (myid .eq. root) then
      open(unit=96,file=trim(batch_manifest),status='replace',action='write',iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7753
      else
        write(96,'(A)') 'status Y1_GV_BATCH_RUNNING'
        write(96,'(A)') 'sector CORE_FROM_VALENCE'
        write(96,'(A)') 'gmres_input bse2.in'
        write(96,'(A,1X,A3)') 'component',component
        write(96,'(A,1X,A)') 'input_manifest',trim(input_manifest)
        write(96,'(A,1X,I0)') 'full_energy_count',info%energy_count
        write(96,'(A,2(1X,I0))') 'selected_index_range',start_index,end_index
        write(96,'(A,1X,ES24.16)') 'eta_hartree',eta_ha
        write(96,'(A,1X,A)') 'paired_dipole_vector',trim(dipole_vector)
        write(96,'(A)') '# Y1_GV index omega_eV Omega_eV norm residual physical_residual iterations binary manifest'
        close(96)
      endif
    endif
    call MPI_BCAST(local_ierr,1,MPI_INTEGER,root,comm,bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    else
      ierr = local_ierr
    endif
  end subroutine begin_batch_manifest

  subroutine append_batch_entry(ierr)
    integer, intent(inout) :: ierr
    integer :: ios, local_ierr, bcast_ierr
    local_ierr = 0
    if (myid .eq. root) then
      open(unit=96,file=trim(batch_manifest),status='old',position='append',action='write',iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7754
      else
        write(96,'(A,1X,I0,5(1X,ES24.16),1X,I0,2(1X,A))',iostat=ios) 'Y1_GV',source_index, &
          info%omega_ev(source_index),target_energy_ha*Hartree2eV,y_norm,residual(1), &
          physical_residual,iterations(1),trim(output_vector),trim(output_manifest)
        close(96)
        if (ios .ne. 0) local_ierr = 7754
      endif
    endif
    call MPI_BCAST(local_ierr,1,MPI_INTEGER,root,comm,bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    else
      ierr = local_ierr
    endif
  end subroutine append_batch_entry

  subroutine complete_batch_manifest(ierr)
    integer, intent(inout) :: ierr
    integer :: ios, local_ierr, bcast_ierr
    local_ierr = 0
    if (myid .eq. root) then
      open(unit=96,file=trim(batch_manifest),status='old',position='append',action='write',iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7755
      else
        write(96,'(A,1X,I0)',iostat=ios) 'completed_count',end_index-start_index+1
        if (ios .eq. 0) write(96,'(A)',iostat=ios) 'completion Y1_GV_BATCH_COMPLETE'
        close(96)
        if (ios .ne. 0) local_ierr = 7755
      endif
    endif
    call MPI_BCAST(local_ierr,1,MPI_INTEGER,root,comm,bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    else
      ierr = local_ierr
    endif
  end subroutine complete_batch_manifest

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
             status='replace', iostat=ios)
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
           status='replace', iostat=ios)
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
           status='replace', iostat=ios)
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
        write(94,'(A)') 'gmres_input bse2.in'
        write(94,'(A,1X,I0)') 'full_energy_count', input%energy_count
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
        write(94,'(A,1X,ES24.16)') 'bse2_in_eta_hartree', configured_eta
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
