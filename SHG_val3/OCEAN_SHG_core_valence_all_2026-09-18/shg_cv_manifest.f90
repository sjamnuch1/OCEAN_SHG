! Manifest reader and compatibility checks for the core <- valence SHG vertex.
module SHG_rbxc_cv_manifest_io
  use AI_kinds, only : DP
  implicit none
  private

  integer, parameter :: path_len = 512

  type, public :: SHG_rbxc_cv_manifest
    character(len=32)  :: status = ''
    character(len=32)  :: sector = ''
    character(len=32)  :: completion = ''
    character(len=128) :: operator_label = ''
    character(len=128) :: binary_layout = ''
    character(len=64)  :: units = ''
    character(len=64)  :: electron_frame = ''
    character(len=1)   :: middle_axis = '?'
    character(len=1)   :: right_axis = '?'
    character(len=32)  :: element = ''
    integer :: site_index = 0
    integer :: energy_count = 0
    integer :: kmesh(3) = 0
    integer :: znl(3) = 0
    integer :: brange(4) = 0
    integer :: ncon = 0
    integer :: nk = 0
    integer :: nalpha = 0
    real(DP) :: eta_valence_hartree = 0.0_DP
    real(DP) :: q_reduced(3) = 0.0_DP
    real(DP), allocatable :: omega_hartree(:)
    real(DP), allocatable :: omega_ev(:)
    real(DP), allocatable :: rbxc_norm(:)
    character(len=path_len), allocatable :: rbxc_file(:)
  end type SHG_rbxc_cv_manifest

  type, public :: SHG_basis_match
    character(len=32) :: status = ''
    character(len=32) :: method = ''
    character(len=64) :: valence_sha256 = ''
    character(len=64) :: core_sha256 = ''
    integer :: ncon = 0
    integer :: nk = 0
    integer :: kmesh(3) = 0
    integer :: conduction_brange(2) = 0
  end type SHG_basis_match

  public :: SHG_read_rbxc_cv_manifest, SHG_find_rbxc_cv_vector
  public :: SHG_read_basis_match, SHG_validate_rbxc_cv_manifest

contains

  subroutine SHG_read_rbxc_cv_manifest(path, info, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, &
                          MPI_DOUBLE_PRECISION, MPI_INTEGER
    character(len=*), intent(in) :: path
    type(SHG_rbxc_cv_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr
    integer :: local_ierr, mpi_ierr, n

    local_ierr = 0
    if (myid .eq. root) call read_cv_root(path, info, local_ierr)
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
    call MPI_BCAST(info%completion, len(info%completion), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%operator_label, len(info%operator_label), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%binary_layout, len(info%binary_layout), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%units, len(info%units), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%electron_frame, len(info%electron_frame), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%middle_axis, 1, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%right_axis, 1, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%element, len(info%element), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%site_index, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%energy_count, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%kmesh, 3, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%znl, 3, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%brange, 4, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%ncon, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%nk, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%nalpha, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%eta_valence_hartree, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%q_reduced, 3, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return

    n = info%energy_count
    if (myid .ne. root) call allocate_cv_arrays(info, n, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%omega_hartree, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%omega_ev, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%rbxc_norm, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%rbxc_file, path_len*n, MPI_CHARACTER, root, comm, ierr)
  end subroutine SHG_read_rbxc_cv_manifest


  subroutine read_cv_root(path, info, ierr)
    character(len=*), intent(in) :: path
    type(SHG_rbxc_cv_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr
    character(len=2048) :: line
    character(len=128) :: key, label1, label2, label3
    character(len=path_len) :: filename
    logical, allocatable :: seen(:)
    integer :: unit, ios, idx, count
    real(DP) :: omega_ha, omega_ev, vector_norm

    open(newunit=unit, file=trim(path), form='formatted', status='old', &
         action='read', iostat=ios)
    if (ios .ne. 0) then
      write(6,'(A,1X,A)') 'Cannot open core-valence RbXc manifest:', trim(path)
      ierr = 7801
      return
    endif
    count = 0
    do
      read(unit,'(A)',iostat=ios) line
      if (ios .lt. 0) exit
      if (ios .ne. 0) then
        ierr = 7802
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
        case('completion')
          read(line,*,iostat=ios) key, info%completion
        case('operator')
          read(line,*,iostat=ios) key, info%operator_label
        case('binary_layout')
          read(line,*,iostat=ios) key, info%binary_layout
        case('units')
          read(line,*,iostat=ios) key, info%units
        case('electron_frame')
          read(line,*,iostat=ios) key, info%electron_frame
        case('middle_axis')
          read(line,*,iostat=ios) key, info%middle_axis
        case('right_axis')
          read(line,*,iostat=ios) key, info%right_axis
        case('dimensions')
          read(line,*,iostat=ios) key, label1, label2, label3, &
                                   info%ncon, info%nk, info%nalpha
        case('input_brange')
          read(line,*,iostat=ios) key, info%brange
        case('kmesh')
          read(line,*,iostat=ios) key, info%kmesh
        case('ZNL')
          read(line,*,iostat=ios) key, info%znl
        case('site')
          read(line,*,iostat=ios) key, info%element, info%site_index
        case('eta_hartree_q_reduced')
          read(line,*,iostat=ios) key, info%eta_valence_hartree, info%q_reduced
        case('energy_count')
          read(line,*,iostat=ios) key, info%energy_count
          if (ios .eq. 0) then
            if (info%energy_count .lt. 1) then
              ios = 1
            else
              call allocate_cv_arrays(info, info%energy_count, ierr)
              if (ierr .ne. 0) exit
              allocate(seen(info%energy_count), stat=ierr)
              if (ierr .ne. 0) exit
              seen = .false.
            endif
          endif
        case('RbXc')
          if (.not. allocated(seen)) then
            ios = 1
          else
            read(line,*,iostat=ios) key, idx, omega_ha, omega_ev, vector_norm, filename
            if (ios .eq. 0) then
              if (idx .lt. 1 .or. idx .gt. info%energy_count .or. seen(idx)) then
                ios = 1
              else
                info%omega_hartree(idx) = omega_ha
                info%omega_ev(idx) = omega_ev
                info%rbxc_norm(idx) = vector_norm
                info%rbxc_file(idx) = filename
                seen(idx) = .true.
                count = count+1
              endif
            endif
          endif
        case default
          ios = 0
      end select
      if (ios .ne. 0) then
        write(6,'(A,1X,A)') 'Malformed core-valence manifest line:', trim(line)
        ierr = 7803
        exit
      endif
    enddo
    close(unit)

    info%middle_axis = lower_axis(info%middle_axis)
    info%right_axis = lower_axis(info%right_axis)

    if (ierr .eq. 0) then
      if (trim(info%status) .ne. 'RUNNING' .or. &
          trim(info%completion) .ne. 'RBXC_CV_COMPLETE' .or. &
          trim(info%sector) .ne. 'CORE_FROM_VALENCE') then
        write(6,'(A)') 'The input is not a completed CORE_FROM_VALENCE RbXc manifest.'
        ierr = 7804
      elseif (.not. allocated(seen)) then
        write(6,'(A)') 'The core-valence manifest has no energy table.'
        ierr = 7805
      elseif (count .ne. info%energy_count .or. .not. all(seen)) then
        write(6,'(A)') 'The core-valence manifest has missing RbXc entries.'
        ierr = 7805
      endif
    endif
    if (allocated(seen)) deallocate(seen)
  end subroutine read_cv_root


  subroutine allocate_cv_arrays(info, n, ierr)
    type(SHG_rbxc_cv_manifest), intent(inout) :: info
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    if (allocated(info%omega_hartree)) deallocate(info%omega_hartree)
    if (allocated(info%omega_ev)) deallocate(info%omega_ev)
    if (allocated(info%rbxc_norm)) deallocate(info%rbxc_norm)
    if (allocated(info%rbxc_file)) deallocate(info%rbxc_file)
    allocate(info%omega_hartree(n), info%omega_ev(n), info%rbxc_norm(n), &
             info%rbxc_file(n), stat=ierr)
    if (ierr .ne. 0) return
    info%omega_hartree = 0.0_DP
    info%omega_ev = 0.0_DP
    info%rbxc_norm = 0.0_DP
    info%rbxc_file = ''
  end subroutine allocate_cv_arrays


  subroutine SHG_read_basis_match(path, basis, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, MPI_INTEGER
    character(len=*), intent(in) :: path
    type(SHG_basis_match), intent(inout) :: basis
    integer, intent(inout) :: ierr
    integer :: local_ierr, mpi_ierr

    local_ierr = 0
    if (myid .eq. root) call read_basis_root(path, basis, local_ierr)
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, mpi_ierr)
    if (mpi_ierr .ne. 0) then
      ierr = mpi_ierr
      return
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
      return
    endif
    call MPI_BCAST(basis%status, len(basis%status), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%method, len(basis%method), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%valence_sha256, 64, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%core_sha256, 64, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%ncon, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%nk, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%kmesh, 3, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(basis%conduction_brange, 2, MPI_INTEGER, root, comm, ierr)
  end subroutine SHG_read_basis_match


  subroutine read_basis_root(path, basis, ierr)
    character(len=*), intent(in) :: path
    type(SHG_basis_match), intent(inout) :: basis
    integer, intent(inout) :: ierr
    character(len=1024) :: line
    character(len=128) :: key
    integer :: unit, ios

    open(newunit=unit, file=trim(path), form='formatted', status='old', &
         action='read', iostat=ios)
    if (ios .ne. 0) then
      write(6,'(A,1X,A)') 'Cannot open conduction-basis certificate:', trim(path)
      ierr = 7810
      return
    endif
    do
      read(unit,'(A)',iostat=ios) line
      if (ios .lt. 0) exit
      if (ios .ne. 0) then
        ierr = 7811
        exit
      endif
      if (len_trim(line) .eq. 0) cycle
      read(line,*,iostat=ios) key
      if (ios .ne. 0) cycle
      select case(trim(key))
        case('status'); read(line,*,iostat=ios) key, basis%status
        case('method'); read(line,*,iostat=ios) key, basis%method
        case('valence_sha256'); read(line,*,iostat=ios) key, basis%valence_sha256
        case('core_sha256'); read(line,*,iostat=ios) key, basis%core_sha256
        case('ncon'); read(line,*,iostat=ios) key, basis%ncon
        case('nk'); read(line,*,iostat=ios) key, basis%nk
        case('kmesh'); read(line,*,iostat=ios) key, basis%kmesh
        case('conduction_brange'); read(line,*,iostat=ios) key, basis%conduction_brange
        case default; ios = 0
      end select
      if (ios .ne. 0) then
        write(6,'(A,1X,A)') 'Malformed basis-certificate line:', trim(line)
        ierr = 7812
        exit
      endif
    enddo
    close(unit)
  end subroutine read_basis_root


  subroutine SHG_validate_rbxc_cv_manifest(sys, info, basis, component, ierr)
    use OCEAN_mpi, only : myid, root
    use OCEAN_system, only : o_system
    type(o_system), intent(in) :: sys
    type(SHG_rbxc_cv_manifest), intent(in) :: info
    type(SHG_basis_match), intent(in) :: basis
    character(len=3), intent(in) :: component
    integer, intent(inout) :: ierr
    integer :: j

    if (.not. associated(sys%cur_run)) then
      ierr = 7820
    elseif (sys%nruns .ne. 1) then
      ierr = 7821
    elseif (.not. sys%cur_run%have_core .or. sys%cur_run%calc_type .ne. 'XAS') then
      ierr = 7822
    elseif (trim(info%operator_label) .ne. 'R_CV=-conjg(xes)_hole') then
      ierr = 7823
    elseif (trim(info%binary_layout) .ne. &
            'raw_complex128_Fortran_ncon_nk_nalpha_no_padding') then
      ierr = 7824
    elseif (trim(info%units) .ne. 'bohr_squared_per_hartree') then
      ierr = 7825
    elseif (trim(info%electron_frame) .ne. 'unchanged_from_valence_k_plus_q') then
      ierr = 7826
    elseif (any(info%kmesh .ne. sys%kmesh) .or. info%nk .ne. sys%nkpts) then
      ierr = 7827
    elseif (maxval(abs(info%q_reduced-sys%qinunitsofbvectors)) .gt. &
            1.0e-10_DP) then
      ierr = 7839
    elseif (info%ncon .ne. sys%cur_run%num_bands .or. info%nalpha .ne. sys%nalpha) then
      ierr = 7828
    elseif (any(info%znl .ne. sys%cur_run%ZNL)) then
      ierr = 7829
    elseif (trim(info%element) .ne. trim(sys%cur_run%elname) .or. &
            info%site_index .ne. sys%cur_run%indx) then
      ierr = 7830
    elseif (any(info%brange(3:4) .ne. sys%brange(3:4))) then
      ierr = 7831
    elseif (component(2:2) .ne. info%middle_axis .or. &
            component(3:3) .ne. info%right_axis) then
      ierr = 7832
    elseif (trim(basis%status) .ne. 'BASIS_MATCH_COMPLETE' .or. &
            trim(basis%method) .ne. 'SHA256_EXACT_FILE_MATCH') then
      ierr = 7833
    elseif (basis%valence_sha256 .ne. basis%core_sha256 .or. &
            len_trim(basis%core_sha256) .ne. 64) then
      ierr = 7834
    elseif (basis%ncon .ne. info%ncon .or. basis%nk .ne. info%nk .or. &
            any(basis%kmesh .ne. info%kmesh) .or. &
            any(basis%conduction_brange .ne. info%brange(3:4))) then
      ierr = 7835
    endif
    do j = 1, 3
      if (.not. valid_axis(component(j:j))) ierr = 7836
    enddo
    if (ierr .eq. 0) then
      if (.not. allocated(info%rbxc_file) .or. &
          size(info%rbxc_file) .ne. info%energy_count) ierr = 7837
    endif
    if (ierr .ne. 0 .and. myid .eq. root) then
      write(6,'(A,I0)') 'CV RbXc/basis does not match the staged core XAS run, ierr = ', ierr
    endif
  end subroutine SHG_validate_rbxc_cv_manifest


  subroutine SHG_find_rbxc_cv_vector(info, path, index, ierr)
    use OCEAN_mpi, only : myid, root
    type(SHG_rbxc_cv_manifest), intent(in) :: info
    character(len=*), intent(in) :: path
    integer, intent(out) :: index
    integer, intent(inout) :: ierr
    character(len=path_len) :: wanted, listed
    integer :: i, matches

    index = 0
    do i = 1, info%energy_count
      if (trim(path) .eq. trim(info%rbxc_file(i))) then
        index = i
        return
      endif
    enddo
    wanted = basename(path)
    matches = 0
    do i = 1, info%energy_count
      listed = basename(info%rbxc_file(i))
      if (trim(wanted) .eq. trim(listed)) then
        index = i
        matches = matches+1
      endif
    enddo
    if (matches .ne. 1) then
      if (myid .eq. root) write(6,'(A,1X,A)') &
        'RbXc dump is not uniquely listed in the CV manifest:', trim(path)
      ierr = 7838
      index = 0
    endif
  end subroutine SHG_find_rbxc_cv_vector


  pure logical function valid_axis(axis)
    character(len=1), intent(in) :: axis
    valid_axis = axis .eq. 'x' .or. axis .eq. 'y' .or. axis .eq. 'z'
  end function valid_axis


  pure function lower_axis(axis) result(lower)
    character(len=1), intent(in) :: axis
    character(len=1) :: lower
    select case(axis)
      case('x','X'); lower = 'x'
      case('y','Y'); lower = 'y'
      case('z','Z'); lower = 'z'
      case default; lower = '?'
    end select
  end function lower_axis


  pure function basename(path) result(base)
    character(len=*), intent(in) :: path
    character(len=path_len) :: base
    integer :: i, first
    first = 1
    do i = len_trim(path), 1, -1
      if (path(i:i) .eq. '/') then
        first = i+1
        exit
      endif
    enddo
    base = ''
    if (first .le. len_trim(path)) base = path(first:len_trim(path))
  end function basename

end module SHG_rbxc_cv_manifest_io
