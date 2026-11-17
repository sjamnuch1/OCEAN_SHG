! Copyright (C) 2026 OCEAN collaboration
!
! Standalone second-resolvent solve for the first core-SHG time ordering:
!
!   [(Omega+i*eta) I - H_BSE^core] |Y_bc(Omega;omega)>
!       = |R_b X_c^-(omega)>
!
! For the SHG T1 term, Omega=2*omega and
!
!   T1_abc(omega) = <d_core^a | Y_bc(2*omega;omega)>.
!
! The input and output vectors are raw complex(DP) streams in the ordering
! returned by OCEAN_psi_vtor.  The companion RbXc manifest supplies the
! dimensions, core edge/site, tensor axes, source energy, and physical norm.

module SHG_rbxc_manifest_io
  use AI_kinds, only : DP

  implicit none
  private

  integer, parameter :: path_len = 512

  type, public :: SHG_rbxc_manifest
    character(len=32) :: status = ''
    character(len=32) :: sector = ''
    character(len=3) :: component = '???'
    character(len=1) :: middle_axis = '?'
    character(len=1) :: right_axis = '?'
    character(len=2) :: element = '__'
    character(len=2) :: corelevel = '__'
    character(len=96) :: binary_layout = ''
    character(len=128) :: operator_label = ''
    integer :: site_index = 0
    integer :: photon = 0
    integer :: energy_count = 0
    integer :: kmesh(3) = 0
    integer :: znl(3) = 0
    integer :: ncon = 0
    integer :: nk = 0
    integer :: nalpha = 0
    real(DP) :: eta_hartree = 0.0_DP
    real(DP), allocatable :: omega_hartree(:)
    real(DP), allocatable :: omega_ev(:)
    real(DP), allocatable :: source_norm(:)
    real(DP), allocatable :: rbxc_norm(:)
    character(len=path_len), allocatable :: xc_file(:)
    character(len=path_len), allocatable :: rbxc_file(:)
  end type SHG_rbxc_manifest

  public :: SHG_read_rbxc_manifest, SHG_validate_rbxc_manifest
  public :: SHG_find_rbxc_vector

contains

  subroutine SHG_read_rbxc_manifest(path, info, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, &
                          MPI_DOUBLE_PRECISION, MPI_INTEGER

    character(len=*), intent(in) :: path
    type(SHG_rbxc_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr

    integer :: local_ierr, mpi_ierr, n

    local_ierr = 0
    if (myid .eq. root) call read_rbxc_manifest_root(path, info, local_ierr)

    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, mpi_ierr)
    if (mpi_ierr .ne. 0) then
      ierr = mpi_ierr
      return
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
      return
    endif

    call MPI_BCAST(info%status, len(info%status), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%sector, len(info%sector), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%component, len(info%component), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%middle_axis, 1, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%right_axis, 1, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%element, len(info%element), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%corelevel, len(info%corelevel), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%binary_layout, len(info%binary_layout), MPI_CHARACTER, &
                   root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%operator_label, len(info%operator_label), MPI_CHARACTER, &
                   root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%site_index, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%photon, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%energy_count, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%kmesh, 3, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%znl, 3, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%ncon, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%nk, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%nalpha, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%eta_hartree, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return

    n = info%energy_count
    if (myid .ne. root) call allocate_manifest_arrays(info, n, ierr)
    if (ierr .ne. 0) return

    call MPI_BCAST(info%omega_hartree, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%omega_ev, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%source_norm, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%rbxc_norm, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%xc_file, path_len*n, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%rbxc_file, path_len*n, MPI_CHARACTER, root, comm, ierr)
  end subroutine SHG_read_rbxc_manifest


  subroutine read_rbxc_manifest_root(path, info, ierr)
    character(len=*), intent(in) :: path
    type(SHG_rbxc_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr

    character(len=2048) :: line
    character(len=128) :: key, word2, word3, word4
    character(len=path_len) :: temp_xc_file, temp_rbxc_file
    character(len=1) :: b_axis, c_axis
    logical, allocatable :: seen(:)
    integer :: ios, idx, count
    real(DP) :: temp_omega_ha, temp_omega_ev
    real(DP) :: temp_source_norm, temp_rbxc_norm

    open(unit=91, file=trim(path), form='formatted', status='old', iostat=ios)
    if (ios .ne. 0) then
      write(6,'(A,1X,A)') 'Failed to open RbXc manifest:', trim(path)
      ierr = 7601
      return
    endif

    count = 0
    do
      read(91,'(A)',iostat=ios) line
      if (ios .lt. 0) exit
      if (ios .ne. 0) then
        ierr = 7602
        exit
      endif
      if (len_trim(line) .eq. 0) cycle

      key = ''
      read(line,*,iostat=ios) key
      if (ios .ne. 0) cycle

      select case(trim(key))
        case('status')
          read(line,*,iostat=ios) key, info%status
        case('sector')
          read(line,*,iostat=ios) key, info%sector
        case('component')
          read(line,*,iostat=ios) key, info%component
        case('middle_axis')
          read(line,*,iostat=ios) key, info%middle_axis
        case('right_axis')
          read(line,*,iostat=ios) key, info%right_axis
        case('operator')
          read(line,*,iostat=ios) key, info%operator_label
        case('binary_layout')
          read(line,*,iostat=ios) key, info%binary_layout
        case('eta_hartree')
          read(line,*,iostat=ios) key, info%eta_hartree
        case('energy_count')
          read(line,*,iostat=ios) key, info%energy_count
          if (ios .eq. 0) then
            if (info%energy_count .lt. 1) then
              ios = 1
            else
              call allocate_manifest_arrays(info, info%energy_count, ierr)
              if (ierr .ne. 0) exit
              allocate(seen(info%energy_count), stat=ierr)
              if (ierr .ne. 0) exit
              seen = .false.
            endif
          endif
        case('kmesh')
          read(line,*,iostat=ios) key, info%kmesh
        case('core_edge_ZNL')
          read(line,*,iostat=ios) key, info%znl
        case('core_site')
          read(line,*,iostat=ios) key, info%element, info%site_index, &
                                   info%corelevel, info%photon
        case('dimensions')
          read(line,*,iostat=ios) key, word2, word3, word4, &
                                   info%ncon, info%nk, info%nalpha
        case('RbXc')
          if (.not. allocated(info%omega_hartree) .or. .not. allocated(seen)) then
            ios = 1
          else
            read(line,*,iostat=ios) key, b_axis, c_axis, idx, temp_omega_ha, &
                 temp_omega_ev, temp_source_norm, temp_rbxc_norm, &
                 temp_xc_file, temp_rbxc_file
            if (ios .eq. 0) then
              if (idx .lt. 1 .or. idx .gt. info%energy_count) then
                ios = 1
              elseif (b_axis .ne. info%middle_axis .or. &
                      c_axis .ne. info%right_axis .or. seen(idx)) then
                ios = 1
              else
                info%omega_hartree(idx) = temp_omega_ha
                info%omega_ev(idx) = temp_omega_ev
                info%source_norm(idx) = temp_source_norm
                info%rbxc_norm(idx) = temp_rbxc_norm
                info%xc_file(idx) = temp_xc_file
                info%rbxc_file(idx) = temp_rbxc_file
                seen(idx) = .true.
                count = count + 1
              endif
            endif
          endif
        case default
          ! Real-space and link-based RbXc manifests have extra diagnostics.
          ios = 0
      end select

      if (ios .ne. 0) then
        write(6,'(A,1X,A)') 'Malformed line in RbXc manifest:', trim(line)
        ierr = 7603
        exit
      endif
    enddo
    close(91)

    if (allocated(seen)) deallocate(seen)
    if (ierr .ne. 0) return
    if (trim(info%status) .ne. 'RBXC_COMPLETE' .or. &
        trim(info%sector) .ne. 'CORE_XAS') then
      write(6,'(A)') 'The manifest is not a completed core RbXc calculation.'
      ierr = 7604
      return
    endif
    if (count .ne. info%energy_count) then
      write(6,'(A,2(1X,I0))') 'RbXc entry count mismatch:', count, info%energy_count
      ierr = 7605
    endif
  end subroutine read_rbxc_manifest_root


  subroutine allocate_manifest_arrays(info, n, ierr)
    type(SHG_rbxc_manifest), intent(inout) :: info
    integer, intent(in) :: n
    integer, intent(inout) :: ierr

    if (allocated(info%omega_hartree)) deallocate(info%omega_hartree)
    if (allocated(info%omega_ev)) deallocate(info%omega_ev)
    if (allocated(info%source_norm)) deallocate(info%source_norm)
    if (allocated(info%rbxc_norm)) deallocate(info%rbxc_norm)
    if (allocated(info%xc_file)) deallocate(info%xc_file)
    if (allocated(info%rbxc_file)) deallocate(info%rbxc_file)

    allocate(info%omega_hartree(n), info%omega_ev(n), info%source_norm(n), &
             info%rbxc_norm(n), info%xc_file(n), info%rbxc_file(n), stat=ierr)
    if (ierr .ne. 0) return

    info%omega_hartree = 0.0_DP
    info%omega_ev = 0.0_DP
    info%source_norm = 0.0_DP
    info%rbxc_norm = 0.0_DP
    info%xc_file = ''
    info%rbxc_file = ''
  end subroutine allocate_manifest_arrays


  subroutine SHG_validate_rbxc_manifest(sys, info, ierr)
    use OCEAN_mpi, only : myid, root
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    type(SHG_rbxc_manifest), intent(in) :: info
    integer, intent(inout) :: ierr

    integer :: j

    if (.not. associated(sys%cur_run)) then
      ierr = 7610
    elseif (sys%nruns .ne. 1) then
      ierr = 7611
    elseif (.not. sys%cur_run%have_core .or. sys%cur_run%calc_type .ne. 'XAS') then
      ierr = 7612
    elseif (any(info%kmesh .ne. sys%kmesh)) then
      ierr = 7613
    elseif (info%nk .ne. sys%nkpts .or. &
            info%ncon .ne. sys%cur_run%num_bands .or. &
            info%nalpha .ne. sys%nalpha) then
      ierr = 7614
    elseif (any(info%znl .ne. sys%cur_run%ZNL)) then
      ierr = 7615
    elseif (trim(info%element) .ne. trim(sys%cur_run%elname) .or. &
            info%site_index .ne. sys%cur_run%indx .or. &
            trim(info%corelevel) .ne. trim(sys%cur_run%corelevel) .or. &
            info%photon .ne. sys%cur_run%photon) then
      ierr = 7616
    elseif (info%middle_axis .ne. info%component(2:2) .or. &
            info%right_axis .ne. info%component(3:3)) then
      ierr = 7617
    elseif (trim(info%binary_layout) .ne. &
            'raw_complex_DP_stream_in_OCEAN_psi_vtor_order') then
      ierr = 7618
    elseif (info%eta_hartree .le. tiny(1.0_DP)) then
      ierr = 7619
    endif

    do j = 1, 3
      if (.not. valid_axis(info%component(j:j))) ierr = 7620
    enddo

    if (ierr .eq. 0) then
      if (.not. allocated(info%rbxc_file) .or. &
          size(info%rbxc_file) .ne. info%energy_count) then
        ierr = 7621
      endif
    endif

    if (ierr .ne. 0 .and. myid .eq. root) then
      write(6,'(A,I0)') 'RbXc manifest does not match the staged XAS run, ierr = ', ierr
    endif
  end subroutine SHG_validate_rbxc_manifest


  subroutine SHG_find_rbxc_vector(info, path, index, ierr)
    use OCEAN_mpi, only : myid, root

    type(SHG_rbxc_manifest), intent(in) :: info
    character(len=*), intent(in) :: path
    integer, intent(out) :: index
    integer, intent(inout) :: ierr

    character(len=path_len) :: requested_base, listed_base
    integer :: i, matches

    index = 0
    do i = 1, info%energy_count
      if (trim(path) .eq. trim(info%rbxc_file(i))) then
        index = i
        return
      endif
    enddo

    requested_base = basename(path)
    matches = 0
    do i = 1, info%energy_count
      listed_base = basename(info%rbxc_file(i))
      if (trim(requested_base) .eq. trim(listed_base)) then
        matches = matches + 1
        index = i
      endif
    enddo

    if (matches .ne. 1) then
      if (myid .eq. root) then
        write(6,'(A,1X,A)') &
          'RbXc dump is not uniquely listed in the manifest:', trim(path)
      endif
      ierr = 7622
      index = 0
    endif
  end subroutine SHG_find_rbxc_vector


  pure function basename(path) result(base)
    character(len=*), intent(in) :: path
    character(len=path_len) :: base

    integer :: i, first

    first = 1
    do i = len_trim(path), 1, -1
      if (path(i:i) .eq. '/') then
        first = i + 1
        exit
      endif
    enddo
    base = ''
    if (first .le. len_trim(path)) base = path(first:len_trim(path))
  end function basename


  pure logical function valid_axis(axis)
    character(len=1), intent(in) :: axis

    valid_axis = axis .eq. 'x' .or. axis .eq. 'y' .or. axis .eq. 'z'
  end function valid_axis

end module SHG_rbxc_manifest_io


program shg_y1_core
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
  use SHG_rbxc_manifest_io, only : SHG_rbxc_manifest, &
                                    SHG_find_rbxc_vector, &
                                    SHG_read_rbxc_manifest, &
                                    SHG_validate_rbxc_manifest

  implicit none

  integer, parameter :: path_len = 512

  type(o_system) :: sys
  type(SHG_rbxc_manifest) :: info
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
  character(len=path_len) :: input_manifest, input_vector
  character(len=path_len) :: output_vector, output_manifest
  logical :: eta_override, is_shg_target
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

  call read_command_line(input_manifest, input_vector, target_energy_ev, &
                         output_vector, eta_override, eta_ev, ierr)
  if (ierr .ne. 0) goto 900

  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900

  call SHG_read_rbxc_manifest(trim(input_manifest), info, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_validate_rbxc_manifest(sys, info, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_find_rbxc_vector(info, trim(input_vector), source_index, ierr)
  if (ierr .ne. 0) goto 900

  target_energy_ha = target_energy_ev*eV2Hartree
  if (eta_override) then
    eta_ha = eta_ev*eV2Hartree
  else
    eta_ha = info%eta_hartree
    eta_ev = eta_ha*Hartree2eV
  endif
  if (.not. ieee_is_finite(target_energy_ha) .or. &
      .not. ieee_is_finite(eta_ha) .or. eta_ha .le. tiny(1.0_DP)) then
    if (myid .eq. root) write(6,'(A)') 'Target energy must be finite and eta must be positive.'
    ierr = 7701
    goto 900
  endif

  if (len_trim(output_vector) .eq. 0) then
    call default_output_name(info%middle_axis, info%right_axis, source_index, &
                             output_vector)
  endif
  if (len_trim(output_vector) .gt. len(output_manifest)-5) then
    if (myid .eq. root) write(6,'(A)') 'Output path is too long for its companion manifest name.'
    ierr = 7705
    goto 900
  endif
  output_manifest = trim(output_vector)//'.info'
  if (trim(output_vector) .eq. trim(input_vector)) then
    if (myid .eq. root) write(6,'(A)') 'Refusing to overwrite the input RbXc vector.'
    ierr = 7702
    goto 900
  endif

  is_shg_target = abs(target_energy_ha-2.0_DP*info%omega_hartree(source_index)) &
                  .le. max(1.0e-10_DP, 1.0e-8_DP*max(1.0_DP, &
                  abs(target_energy_ha), abs(2.0_DP*info%omega_hartree(source_index))))

  if (myid .eq. root) then
    write(6,'(A)') 'Standalone core-SHG T1 second-resolvent solve'
    write(6,'(A)') '  [(Omega+i*eta)I-H_BSE^core] Y_bc = R_b X_c^-(omega)'
    write(6,'(A,1X,A)') '  RbXc manifest:', trim(input_manifest)
    write(6,'(A,1X,A)') '  RbXc vector:  ', trim(input_vector)
    write(6,'(A,F18.8,A)') '  source omega: ', info%omega_ev(source_index), ' eV'
    write(6,'(A,F18.8,A)') '  target Omega: ', target_energy_ev, ' eV'
    write(6,'(A,F18.8,A)') '  broadening:   ', eta_ev, ' eV'
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
  call write_output_manifest(sys, info, source_index, trim(input_manifest), &
       trim(input_vector), trim(output_vector), trim(output_manifest), &
       target_energy_ha, eta_ha, configured_eta_ha, preconditioner_ha, &
       convergence, gmres_depth, rhs_norm, y_norm, residual(1), &
       physical_residual, iterations(1), is_shg_target, ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A,1X,A)') 'Saved physical Y_bc vector:', trim(output_vector)
    write(6,'(A,1X,A)') 'Saved companion manifest:', trim(output_manifest)
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
      if (myid .eq. root) write(6,'(A,I0)') 'shg_y1_core failed with ierr = ', ierr
    else
      write(6,'(A,I0)') 'shg_y1_core failed during MPI initialization, ierr = ', ierr
    endif
    stop 1
  endif

contains

  subroutine read_command_line(manifest, rbxc_file, target_ev, output_file, &
                               have_eta, eta_value_ev, ierr)
    character(len=*), intent(out) :: manifest, rbxc_file, output_file
    real(DP), intent(out) :: target_ev, eta_value_ev
    logical, intent(out) :: have_eta
    integer, intent(inout) :: ierr

    character(len=128) :: number
    integer :: nargs, ios

    manifest = ''
    rbxc_file = ''
    output_file = ''
    target_ev = 0.0_DP
    eta_value_ev = 0.0_DP
    have_eta = .false.
    nargs = command_argument_count()

    if (nargs .lt. 3 .or. nargs .gt. 5) then
      if (myid .eq. root) call print_usage()
      ierr = 7710
      return
    endif

    call get_path_argument(1, manifest, ierr)
    if (ierr .ne. 0) return
    call get_path_argument(2, rbxc_file, ierr)
    if (ierr .ne. 0) return
    call get_command_argument(3, value=number)
    read(number,*,iostat=ios) target_ev
    if (ios .ne. 0 .or. .not. ieee_is_finite(target_ev)) then
      if (myid .eq. root) write(6,'(A)') 'TARGET_ENERGY_EV is not a finite number.'
      ierr = 7711
      return
    endif

    if (nargs .ge. 4) then
      call get_path_argument(4, output_file, ierr)
      if (ierr .ne. 0) return
    endif
    if (nargs .eq. 5) then
      call get_command_argument(5, value=number)
      read(number,*,iostat=ios) eta_value_ev
      if (ios .ne. 0 .or. .not. ieee_is_finite(eta_value_ev) .or. &
          eta_value_ev .le. 0.0_DP) then
        if (myid .eq. root) write(6,'(A)') 'ETA_EV must be a positive finite number.'
        ierr = 7712
        return
      endif
      have_eta = .true.
    endif
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
    write(6,'(A)') '  shg_y1_core.x RBXC_MANIFEST RBXC_DUMP TARGET_ENERGY_EV [OUTPUT_BIN] [ETA_EV]'
    write(6,'(A)') 'For T1_abc, set TARGET_ENERGY_EV = 2*omega_eV of the selected RbXc dump.'
    write(6,'(A)') 'When omitted, ETA_EV is read from RBXC_MANIFEST.'
  end subroutine print_usage


  subroutine default_output_name(b_axis, c_axis, index, path)
    character(len=1), intent(in) :: b_axis, c_axis
    integer, intent(in) :: index
    character(len=*), intent(out) :: path

    write(path,'(A,A1,A1,A,I6.6,A)') &
      'shg_cache_Y1_core_', b_axis, c_axis, '_', index, '.bin'
  end subroutine default_output_name


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


  subroutine write_output_manifest(sys, input, index, input_manifest_path, &
       input_vector_path, output_vector_path, manifest_path, target_ha, eta, &
       configured_eta, preconditioner, tolerance, depth, rbxc_norm, &
       output_norm, normalized_residual, physical_residual_value, niter, &
       target_is_two_omega, ierr)
    type(o_system), intent(in) :: sys
    type(SHG_rbxc_manifest), intent(in) :: input
    integer, intent(in) :: index, depth, niter
    character(len=*), intent(in) :: input_manifest_path, input_vector_path
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
        write(94,'(A)') 'status Y1_COMPLETE'
        write(94,'(A)') 'sector CORE_XAS'
        write(94,'(A,1X,A3)') 'component', input%component
        write(94,'(A,1X,A1)') 'left_axis', input%component(1:1)
        write(94,'(A,1X,A1)') 'middle_axis', input%middle_axis
        write(94,'(A,1X,A1)') 'right_axis', input%right_axis
        write(94,'(A)') 'equation Y1_bc=Gminus_core(Omega)*RbXc'
        write(94,'(A)') 'linear_system [(Omega+i*eta)I-H_BSE_core]Y1_bc=RbXc'
        write(94,'(A)') 't1_contraction T1_abc=<d_core^a|Y1_bc_at_Omega=2omega>'
        write(94,'(A)') 'hamiltonian_sign -1'
        write(94,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
        write(94,'(A,1X,A)') 'input_rbxc_operator', trim(input%operator_label)
        write(94,'(A,1X,A)') 'input_rbxc_binary_layout', trim(input%binary_layout)
        write(94,'(A,1X,A)') 'input_manifest', trim(input_manifest_path)
        write(94,'(A,1X,A)') 'input_vector', trim(input_vector_path)
        write(94,'(A,1X,A)') 'output_vector', trim(output_vector_path)
        write(94,'(A,1X,I0)') 'energy_count', 1
        write(94,'(A,1X,I0)') 'source_index', index
        write(94,'(A,1X,ES24.16)') 'source_omega_hartree', input%omega_hartree(index)
        write(94,'(A,1X,ES24.16)') 'source_omega_ev', input%omega_ev(index)
        write(94,'(A,1X,ES24.16)') 'target_energy_hartree', target_ha
        write(94,'(A,1X,ES24.16)') 'target_energy_ev', target_ha*Hartree2eV
        write(94,'(A,1X,L1)') 'target_is_two_omega', target_is_two_omega
        write(94,'(A,1X,ES24.16)') 'eta_hartree', eta
        write(94,'(A,1X,ES24.16)') 'eta_ev', eta*Hartree2eV
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
        write(94,'(A,3(1X,I0))') 'core_edge_ZNL', input%znl
        write(94,'(A,1X,A,1X,I0,1X,A,1X,I0)') 'core_site', &
             trim(input%element), input%site_index, trim(input%corelevel), input%photon
        write(94,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(94,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
             sys%cur_run%num_bands, sys%nkpts, sys%nalpha
        write(94,'(A,1X,A1,1X,A1,1X,I0,7(1X,ES24.16),1X,I0,1X,A)') &
             'Y1', input%middle_axis, input%right_axis, index, &
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

end program shg_y1_core
