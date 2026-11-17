! Standalone final contraction for the core-valence SHG time ordering.
!
!   T1_gv_abc(omega) = <d_core^a | Y1_bc(2 omega;omega)>
!                 = sum_i conjg(d_a(i))*Y_bc(i)
!
! Both vectors must be raw complex(real64) streams written in
! OCEAN_psi_vtor order.  This program is serial and has no OCEAN, MPI,
! BLAS, or FFT dependency.
program shg_t1_gv_contract
  use, intrinsic :: iso_fortran_env, only : error_unit, int64, real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer, parameter :: path_len = 1024
  integer, parameter :: word_len = 128
  character(len=*), parameter :: expected_layout = &
       'raw_complex_DP_stream_in_OCEAN_psi_vtor_order'

  type :: manifest_info
    character(len=word_len) :: status = ''
    character(len=word_len) :: sector = ''
    character(len=word_len) :: binary_layout = ''
    character(len=3) :: component = '???'
    character(len=1) :: left_axis = '?'
    character(len=1) :: middle_axis = '?'
    character(len=1) :: right_axis = '?'
    character(len=word_len) :: element = ''
    character(len=word_len) :: corelevel = ''
    character(len=path_len) :: dipole_file = ''
    character(len=path_len) :: output_vector = ''
    integer :: kmesh(3) = 0
    integer :: znl(3) = 0
    integer :: ncon = 0
    integer :: nk = 0
    integer :: nalpha = 0
    integer :: site_index = 0
    integer :: photon = 0
    integer :: source_index = 0
    real(real64) :: dipole_norm = -1.0_real64
    real(real64) :: y1_norm = -1.0_real64
    real(real64) :: source_omega_hartree = 0.0_real64
    real(real64) :: source_omega_ev = 0.0_real64
    real(real64) :: target_energy_hartree = 0.0_real64
    real(real64) :: target_energy_ev = 0.0_real64
    real(real64) :: eta_hartree = 0.0_real64
    logical :: target_is_two_omega = .false.
    logical :: have_target_flag = .false.
    logical :: have_dimensions = .false.
    logical :: have_kmesh = .false.
    logical :: have_znl = .false.
    logical :: have_core_site = .false.
  end type manifest_info

  type(manifest_info) :: d_info, y_info
  character(len=path_len) :: d_manifest, d_file
  character(len=path_len) :: y_manifest, y_file, result_file
  complex(real64), allocatable :: d(:), y(:)
  complex(real64) :: t1
  real(real64) :: d_norm, y_norm, overlap_fraction
  integer(int64) :: nstored, nphysical
  logical :: padded_storage

  call read_command_line(d_manifest, d_file, y_manifest, y_file, result_file)
  call read_manifest(trim(d_manifest), d_info)
  call read_manifest(trim(y_manifest), y_info)
  call validate_manifests(d_info, y_info)

  nphysical = physical_size(d_info)
  call read_matching_streams(trim(d_file), trim(y_file), nphysical, &
                             d, y, nstored, padded_storage)

  call check_finite_vector('dipole', d)
  call check_finite_vector('Y1', y)
  call contract_vectors(d, y, t1, d_norm, y_norm)

  if (d_norm <= tiny(1.0_real64)) call fatal('The dipole vector has zero norm.')
  if (y_norm <= tiny(1.0_real64)) call fatal('The Y1 vector has zero norm.')
  overlap_fraction = abs(t1)/(d_norm*y_norm)

  call check_manifest_norm('dipole', d_norm, d_info%dipole_norm)
  call check_manifest_norm('Y1', y_norm, y_info%y1_norm)
  call check_listed_file('dipole', trim(d_file), trim(d_info%dipole_file))
  call check_listed_file('Y1', trim(y_file), trim(y_info%output_vector))

  call write_result(trim(result_file), d_info, y_info, trim(d_manifest), &
                    trim(d_file), trim(y_manifest), trim(y_file), nphysical, &
                    nstored, padded_storage, d_norm, y_norm, t1, &
                    overlap_fraction)

  write(*,'(A)') 'Standalone core-valence SHG T1 contraction complete'
  write(*,'(A,1X,A3)') '  component:', y_info%component
  write(*,'(A,1X,ES24.16)') '  Re(T1):', real(t1,real64)
  write(*,'(A,1X,ES24.16)') '  Im(T1):', aimag(t1)
  write(*,'(A,1X,A)') '  output:', trim(result_file)

contains

  subroutine read_command_line(dm, df, ym, yf, output)
    character(len=*), intent(out) :: dm, df, ym, yf, output
    integer :: nargs

    dm = ''
    df = ''
    ym = ''
    yf = ''
    output = ''
    nargs = command_argument_count()
    if (nargs < 4 .or. nargs > 5) then
      write(error_unit,'(A)') 'Usage:'
      write(error_unit,'(A)') &
        '  shg_t1_gv_contract.x DIPOLE_MANIFEST DIPOLE_BIN Y1_MANIFEST Y1_BIN [OUTPUT_DAT]'
      error stop 2
    endif

    call get_path_argument(1, dm)
    call get_path_argument(2, df)
    call get_path_argument(3, ym)
    call get_path_argument(4, yf)
    if (nargs == 5) call get_path_argument(5, output)
    if (len_trim(output) == 0) output = 'shg_T1_gv.dat'

    if (trim(output) == trim(dm) .or. trim(output) == trim(df) .or. &
        trim(output) == trim(ym) .or. trim(output) == trim(yf)) then
      call fatal('Refusing to overwrite an input file.')
    endif
  end subroutine read_command_line


  subroutine get_path_argument(position, value)
    integer, intent(in) :: position
    character(len=*), intent(out) :: value
    integer :: nchar

    value = ''
    call get_command_argument(position, length=nchar)
    if (nchar < 1) call fatal('A required path argument is empty.')
    if (nchar > len(value)) call fatal('A path argument exceeds the internal path limit.')
    call get_command_argument(position, value=value)
  end subroutine get_path_argument


  subroutine read_manifest(path, info)
    character(len=*), intent(in) :: path
    type(manifest_info), intent(out) :: info

    character(len=4096) :: line
    character(len=word_len) :: key, label1, label2, label3
    integer :: unit, ios

    open(newunit=unit, file=path, form='formatted', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) call fatal('Cannot open manifest: '//trim(path))

    do
      read(unit,'(A)',iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) then
        close(unit)
        call fatal('Error while reading manifest: '//trim(path))
      endif
      if (len_trim(line) == 0) cycle
      key = ''
      read(line,*,iostat=ios) key
      if (ios /= 0) cycle

      select case(trim(key))
        case('status')
          read(line,*,iostat=ios) key, info%status
        case('sector')
          read(line,*,iostat=ios) key, info%sector
        case('component')
          read(line,*,iostat=ios) key, info%component
        case('left_axis')
          read(line,*,iostat=ios) key, info%left_axis
        case('middle_axis')
          read(line,*,iostat=ios) key, info%middle_axis
        case('right_axis')
          read(line,*,iostat=ios) key, info%right_axis
        case('binary_layout')
          read(line,*,iostat=ios) key, info%binary_layout
        case('kmesh')
          read(line,*,iostat=ios) key, info%kmesh
          if (ios == 0) info%have_kmesh = .true.
        case('core_edge_ZNL')
          read(line,*,iostat=ios) key, info%znl
          if (ios == 0) info%have_znl = .true.
        case('core_site')
          read(line,*,iostat=ios) key, info%element, info%site_index, &
                                   info%corelevel, info%photon
          if (ios == 0) info%have_core_site = .true.
        case('dimensions')
          read(line,*,iostat=ios) key, label1, label2, label3, &
                                   info%ncon, info%nk, info%nalpha
          if (ios == 0) info%have_dimensions = .true.
        case('dipole')
          read(line,*,iostat=ios) key, info%dipole_norm, info%dipole_file
        case('output_vector')
          read(line,*,iostat=ios) key, info%output_vector
        case('source_index')
          read(line,*,iostat=ios) key, info%source_index
        case('source_omega_hartree')
          read(line,*,iostat=ios) key, info%source_omega_hartree
        case('source_omega_ev')
          read(line,*,iostat=ios) key, info%source_omega_ev
        case('target_energy_hartree')
          read(line,*,iostat=ios) key, info%target_energy_hartree
        case('target_energy_ev')
          read(line,*,iostat=ios) key, info%target_energy_ev
        case('target_is_two_omega')
          read(line,*,iostat=ios) key, info%target_is_two_omega
          if (ios == 0) info%have_target_flag = .true.
        case('eta_hartree')
          read(line,*,iostat=ios) key, info%eta_hartree
        case('y1_norm')
          read(line,*,iostat=ios) key, info%y1_norm
        case default
          ios = 0
      end select

      if (ios /= 0) then
        close(unit)
        call fatal('Malformed manifest line: '//trim(line))
      endif
    enddo
    close(unit)

    info%left_axis = lower_axis(info%left_axis)
    info%middle_axis = lower_axis(info%middle_axis)
    info%right_axis = lower_axis(info%right_axis)
    info%component(1:1) = lower_axis(info%component(1:1))
    info%component(2:2) = lower_axis(info%component(2:2))
    info%component(3:3) = lower_axis(info%component(3:3))
  end subroutine read_manifest


  subroutine validate_manifests(dmeta, ymeta)
    type(manifest_info), intent(in) :: dmeta, ymeta

    if (trim(dmeta%status) /= 'XC_MINUS_COMPLETE') &
      call fatal('The dipole manifest is not an XC_MINUS_COMPLETE manifest.')
    if (trim(ymeta%status) /= 'Y1_GV_COMPLETE') &
      call fatal('The Y manifest is not a Y1_GV_COMPLETE manifest.')
    if (trim(dmeta%sector) /= 'CORE_XAS') &
      call fatal('The dipole manifest must describe the CORE_XAS sector.')
    if (trim(ymeta%sector) /= 'CORE_FROM_VALENCE') &
      call fatal('The Y manifest must describe the CORE_FROM_VALENCE sector.')
    if (trim(dmeta%binary_layout) /= expected_layout .or. &
        trim(ymeta%binary_layout) /= expected_layout) &
      call fatal('Unsupported binary layout in one or both manifests.')

    if (.not. dmeta%have_dimensions .or. .not. ymeta%have_dimensions) &
      call fatal('Both manifests must contain dimensions ncon nk nalpha.')
    if (dmeta%ncon <= 0 .or. dmeta%nk <= 0 .or. dmeta%nalpha <= 0) &
      call fatal('The physical vector dimensions must be positive.')
    if (dmeta%ncon /= ymeta%ncon .or. dmeta%nk /= ymeta%nk .or. &
        dmeta%nalpha /= ymeta%nalpha) &
      call fatal('Dipole and Y dimensions do not match.')

    if (.not. dmeta%have_kmesh .or. .not. ymeta%have_kmesh) &
      call fatal('Both manifests must contain kmesh.')
    if (any(dmeta%kmesh /= ymeta%kmesh)) &
      call fatal('Dipole and Y k meshes do not match.')
    if (product(int(dmeta%kmesh,int64)) /= int(dmeta%nk,int64)) &
      call fatal('The dipole manifest has nk inconsistent with kmesh.')

    if (.not. dmeta%have_znl .or. .not. ymeta%have_znl) &
      call fatal('Both manifests must identify the core edge.')
    if (any(dmeta%znl /= ymeta%znl)) &
      call fatal('Dipole and Y core edges do not match.')
    if (.not. dmeta%have_core_site .or. .not. ymeta%have_core_site) &
      call fatal('Both manifests must identify the core site.')
    if (trim(dmeta%element) /= trim(ymeta%element) .or. &
        dmeta%site_index /= ymeta%site_index .or. &
        trim(dmeta%corelevel) /= trim(ymeta%corelevel)) &
      call fatal('Dipole and Y core sites do not match.')

    if (.not. valid_axis(dmeta%right_axis)) &
      call fatal('The dipole manifest has an invalid right_axis.')
    if (.not. valid_axis(ymeta%left_axis) .or. &
        .not. valid_axis(ymeta%middle_axis) .or. &
        .not. valid_axis(ymeta%right_axis)) &
      call fatal('The Y manifest has an invalid tensor axis.')
    if (dmeta%right_axis /= ymeta%left_axis) &
      call fatal('d_a polarization does not match the Y manifest left_axis a.')
    if (ymeta%component(1:1) /= ymeta%left_axis .or. &
        ymeta%component(2:2) /= ymeta%middle_axis .or. &
        ymeta%component(3:3) /= ymeta%right_axis) &
      call fatal('The Y component and its a,b,c axis labels are inconsistent.')

    if (dmeta%dipole_norm <= 0.0_real64) &
      call fatal('The dipole manifest does not contain a positive dipole norm.')
    if (ymeta%y1_norm <= 0.0_real64) &
      call fatal('The Y manifest does not contain a positive y1_norm.')

    ! The photon indices need not match: d_a can use a different polarization
    ! file from the d_c source that originally generated Y_bc.
    if (ymeta%have_target_flag .and. .not. ymeta%target_is_two_omega) then
      call warn('Y was not marked as an Omega=2*omega solve; the dot product is '// &
                'still computed, but it is not the standard T1 SHG point.')
    endif
  end subroutine validate_manifests


  function physical_size(info) result(n)
    type(manifest_info), intent(in) :: info
    integer(int64) :: n

    n = int(info%ncon,int64)*int(info%nk,int64)*int(info%nalpha,int64)
    if (n <= 0_int64) call fatal('Invalid or overflowing physical vector size.')
  end function physical_size


  subroutine read_matching_streams(dpath, ypath, nminimum, dvec, yvec, &
                                   n, padded)
    character(len=*), intent(in) :: dpath, ypath
    integer(int64), intent(in) :: nminimum
    complex(real64), allocatable, intent(out) :: dvec(:), yvec(:)
    integer(int64), intent(out) :: n
    logical, intent(out) :: padded

    integer(int64) :: dbytes, ybytes, element_bytes
    integer :: ndefault

    element_bytes = int(storage_size(cmplx(0.0_real64,0.0_real64,real64))/8, &
                        int64)
    if (element_bytes /= 16_int64) &
      call fatal('This compiler does not store complex(real64) in 16 bytes.')

    call stream_size(dpath, dbytes)
    call stream_size(ypath, ybytes)
    if (dbytes /= ybytes) &
      call fatal('Dipole and Y stream byte sizes differ.')
    if (mod(dbytes,element_bytes) /= 0_int64) &
      call fatal('Input stream size is not a multiple of complex(real64).')

    n = dbytes/element_bytes
    if (n < nminimum) &
      call fatal('Stored vector is smaller than ncon*nk*nalpha.')
    if (n > int(huge(0),int64)) &
      call fatal('Stored vector is too large for this serial executable.')
    padded = n /= nminimum
    ndefault = int(n)

    allocate(dvec(ndefault), yvec(ndefault))
    call read_stream(dpath, dvec)
    call read_stream(ypath, yvec)
  end subroutine read_matching_streams


  subroutine stream_size(path, nbytes)
    character(len=*), intent(in) :: path
    integer(int64), intent(out) :: nbytes
    integer :: ios
    logical :: exists

    inquire(file=path, exist=exists, size=nbytes, iostat=ios)
    if (ios /= 0 .or. .not. exists) &
      call fatal('Cannot inspect input stream: '//trim(path))
    if (nbytes <= 0_int64) call fatal('Input stream is empty: '//trim(path))
  end subroutine stream_size


  subroutine read_stream(path, vector)
    character(len=*), intent(in) :: path
    complex(real64), intent(out) :: vector(:)
    integer :: unit, ios

    open(newunit=unit, file=path, form='unformatted', access='stream', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('Cannot open input stream: '//trim(path))
    read(unit,iostat=ios) vector
    close(unit)
    if (ios /= 0) call fatal('Cannot read complete input stream: '//trim(path))
  end subroutine read_stream


  subroutine check_finite_vector(label, vector)
    character(len=*), intent(in) :: label
    complex(real64), intent(in) :: vector(:)
    integer :: i

    do i = 1, size(vector)
      if (.not. ieee_is_finite(real(vector(i),real64)) .or. &
          .not. ieee_is_finite(aimag(vector(i)))) then
        call fatal(trim(label)//' vector contains NaN or infinity.')
      endif
    enddo
  end subroutine check_finite_vector


  subroutine contract_vectors(dvec, yvec, value, norm_d, norm_y)
    complex(real64), intent(in) :: dvec(:), yvec(:)
    complex(real64), intent(out) :: value
    real(real64), intent(out) :: norm_d, norm_y

    real(real64) :: sum_re, sum_im, corr_re, corr_im
    real(real64) :: d2, y2, corr_d2, corr_y2
    real(real64) :: dr, di, yr, yi
    integer :: i

    sum_re = 0.0_real64
    sum_im = 0.0_real64
    corr_re = 0.0_real64
    corr_im = 0.0_real64
    d2 = 0.0_real64
    y2 = 0.0_real64
    corr_d2 = 0.0_real64
    corr_y2 = 0.0_real64

    do i = 1, size(dvec)
      dr = real(dvec(i),real64)
      di = aimag(dvec(i))
      yr = real(yvec(i),real64)
      yi = aimag(yvec(i))

      ! conjg(d)*y = (dr*yr+di*yi) + i*(dr*yi-di*yr)
      call kahan_add(dr*yr+di*yi, sum_re, corr_re)
      call kahan_add(dr*yi-di*yr, sum_im, corr_im)
      call kahan_add(dr*dr+di*di, d2, corr_d2)
      call kahan_add(yr*yr+yi*yi, y2, corr_y2)
    enddo

    value = cmplx(sum_re,sum_im,real64)
    norm_d = sqrt(max(0.0_real64,d2))
    norm_y = sqrt(max(0.0_real64,y2))
  end subroutine contract_vectors


  subroutine kahan_add(x, sum_value, correction)
    real(real64), intent(in) :: x
    real(real64), intent(inout) :: sum_value, correction
    real(real64) :: adjusted, next_sum

    adjusted = x-correction
    next_sum = sum_value+adjusted
    correction = (next_sum-sum_value)-adjusted
    sum_value = next_sum
  end subroutine kahan_add


  subroutine check_manifest_norm(label, calculated, recorded)
    character(len=*), intent(in) :: label
    real(real64), intent(in) :: calculated, recorded
    real(real64) :: relative_error
    character(len=256) :: message

    relative_error = abs(calculated-recorded)/ &
                     max(calculated,recorded,tiny(1.0_real64))
    if (relative_error > 1.0e-8_real64) then
      write(message,'(A,A,A,2(1X,ES16.8))') 'Recomputed ', trim(label), &
        ' norm differs from manifest:', calculated, recorded
      call warn(trim(message))
    endif
  end subroutine check_manifest_norm


  subroutine check_listed_file(label, supplied, listed)
    character(len=*), intent(in) :: label, supplied, listed

    if (len_trim(listed) == 0) then
      call warn(trim(label)//' manifest does not list its vector filename.')
    elseif (trim(basename(supplied)) /= trim(basename(listed))) then
      call warn('Supplied '//trim(label)//' filename differs from the manifest entry.')
    endif
  end subroutine check_listed_file


  function basename(path) result(base)
    character(len=*), intent(in) :: path
    character(len=path_len) :: base
    integer :: i, first

    first = 1
    do i = len_trim(path), 1, -1
      if (path(i:i) == '/') then
        first = i+1
        exit
      endif
    enddo
    base = ''
    if (first <= len_trim(path)) base = path(first:len_trim(path))
  end function basename


  subroutine write_result(path, dmeta, ymeta, dm_path, dpath, ym_path, &
                          ypath, nminimum, n, padded, norm_d, norm_y, &
                          value, fraction)
    character(len=*), intent(in) :: path, dm_path, dpath, ym_path, ypath
    type(manifest_info), intent(in) :: dmeta, ymeta
    integer(int64), intent(in) :: nminimum, n
    logical, intent(in) :: padded
    real(real64), intent(in) :: norm_d, norm_y, fraction
    complex(real64), intent(in) :: value
    integer :: unit, ios

    open(newunit=unit, file=path, form='formatted', status='new', &
         action='write', iostat=ios)
    if (ios /= 0) call fatal('Cannot create result file: '//trim(path))

    write(unit,'(A)') 'status T1_GV_COMPLETE'
    write(unit,'(A)') 'sector CORE_FROM_VALENCE'
    write(unit,'(A)') 'pathway GV'
    write(unit,'(A,1X,A3)') 'component', ymeta%component
    write(unit,'(A,1X,A1)') 'left_axis', ymeta%left_axis
    write(unit,'(A,1X,A1)') 'middle_axis', ymeta%middle_axis
    write(unit,'(A,1X,A1)') 'right_axis', ymeta%right_axis
    write(unit,'(A)') 'equation T1_gv_abc=sum_i(conjg(d_core_a(i))*Y_gv_bc(i))'
    write(unit,'(A)') 'normalization direct_OCEAN_psi_dot_no_extra_1_over_Nk'
    write(unit,'(A)') 'prefactor none_raw_time_ordering_term_not_chi2'
    write(unit,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
    write(unit,'(A,1X,A)') 'dipole_manifest', trim(dm_path)
    write(unit,'(A,1X,A)') 'dipole_vector', trim(dpath)
    write(unit,'(A,1X,A)') 'y1_manifest', trim(ym_path)
    write(unit,'(A,1X,A)') 'y1_vector', trim(ypath)
    write(unit,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
                               dmeta%ncon, dmeta%nk, dmeta%nalpha
    write(unit,'(A,3(1X,I0))') 'kmesh', dmeta%kmesh
    write(unit,'(A,3(1X,I0))') 'core_edge_ZNL', dmeta%znl
    write(unit,'(A,1X,A,1X,I0,1X,A)') 'core_site', trim(dmeta%element), &
                                      dmeta%site_index, trim(dmeta%corelevel)
    write(unit,'(A,1X,I0)') 'physical_complex_count', nminimum
    write(unit,'(A,1X,I0)') 'stored_complex_count', n
    write(unit,'(A,1X,L1)') 'padded_storage', padded
    write(unit,'(A,1X,I0)') 'source_index', ymeta%source_index
    write(unit,'(A,1X,ES24.16)') 'source_omega_hartree', &
                                 ymeta%source_omega_hartree
    write(unit,'(A,1X,ES24.16)') 'source_omega_ev', ymeta%source_omega_ev
    write(unit,'(A,1X,ES24.16)') 'target_energy_hartree', &
                                 ymeta%target_energy_hartree
    write(unit,'(A,1X,ES24.16)') 'target_energy_ev', ymeta%target_energy_ev
    if (ymeta%have_target_flag) &
      write(unit,'(A,1X,L1)') 'target_is_two_omega', ymeta%target_is_two_omega
    write(unit,'(A,1X,ES24.16)') 'eta_hartree', ymeta%eta_hartree
    write(unit,'(A,1X,ES24.16)') 'dipole_norm_recomputed', norm_d
    write(unit,'(A,1X,ES24.16)') 'y1_norm_recomputed', norm_y
    write(unit,'(A,1X,ES24.16)') 'T1_real', real(value,real64)
    write(unit,'(A,1X,ES24.16)') 'T1_imag', aimag(value)
    write(unit,'(A,1X,ES24.16)') 'T1_abs', abs(value)
    write(unit,'(A,1X,ES24.16)') 'T1_phase_radian', atan2(aimag(value), &
                                                        real(value,real64))
    write(unit,'(A,1X,ES24.16)') 'normalized_overlap_abs', fraction
    close(unit,iostat=ios)
    if (ios /= 0) call fatal('Error closing result file: '//trim(path))
  end subroutine write_result


  pure function lower_axis(axis) result(lower)
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
  end function lower_axis


  pure logical function valid_axis(axis)
    character(len=1), intent(in) :: axis
    valid_axis = axis == 'x' .or. axis == 'y' .or. axis == 'z'
  end function valid_axis


  subroutine warn(message)
    character(len=*), intent(in) :: message
    write(error_unit,'(A)') 'WARNING: '//trim(message)
  end subroutine warn


  subroutine fatal(message)
    character(len=*), intent(in) :: message
    write(error_unit,'(A)') 'FATAL: '//trim(message)
    error stop 1
  end subroutine fatal

end program shg_t1_gv_contract
