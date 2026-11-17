! Copyright (C) 2026 OCEAN collaboration
!
! Standalone second-layer source construction for core SHG:
!   |q_b,c(omega)> = R_core^b |X_core,c^-(omega)>
!
! Input vectors are the physical first-layer caches written by shg_xc_minus.x.
! The program links only unmodified OCEAN objects; all SHG-specific work is in
! SHG_core_vertex.
program shg_rbxc_core
  use AI_kinds, only : DP
  use OCEAN_energies, only : OCEAN_energies_allow_full, OCEAN_energies_init, &
                             OCEAN_energies_load
  use OCEAN_mpi, only : myid, root, comm, OCEAN_mpi_init, OCEAN_mpi_finalize, &
                        MPI_BCAST, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, MPI_INTEGER
  use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_init, OCEAN_psi_kill, &
                        OCEAN_psi_new, OCEAN_psi_rtov, OCEAN_psi_size_full, &
                        OCEAN_psi_vtor, OCEAN_psi_zero_full
  use OCEAN_system, only : o_system, OCEAN_sys_init
  use OCEAN_timekeeper, only : OCEAN_tk_init, OCEAN_tk_printtimes, &
                               OCEAN_tk_start, OCEAN_tk_stop, tk_total
  use SHG_core_vertex, only : SHG_xc_manifest, SHG_read_xc_manifest, &
                              SHG_validate_xc_manifest, SHG_load_conduction_states, &
                              SHG_build_forward_links, SHG_apply_core_R_direction, &
                              SHG_axis_index, SHG_cartesian_weight

  implicit none

  type(o_system) :: sys
  type(SHG_xc_manifest) :: info
  type(OCEAN_vector) :: xc, rbxc
  complex(DP), allocatable :: conduction_states(:,:,:,:), uplus(:,:,:,:)
  real(DP), allocatable :: output_norm(:)
  real(DP) :: raw_norm_min, raw_norm_max, weight(3)
  integer :: active_direction(3), nactive, iactive, direction, cart_axis
  integer :: nx_local, x_start, i, ierr, cleanup_ierr, mpi_ierr
  integer :: arg_length
  character(len=256) :: input_manifest, output_manifest
  character(len=256) :: output_file, work_file
  character(len=1) :: b_axis, c_axis
  logical :: mpi_started, timing_started, have_vectors, complex_bse

  ierr = 0
  cleanup_ierr = 0
  mpi_ierr = 0
  mpi_started = .false.
  timing_started = .false.
  have_vectors = .false.
  complex_bse = .false.

  input_manifest = 'shg_xc_minus_core_z.info'
  call get_command_argument(1, length=arg_length)
  if (arg_length .gt. 0) call get_command_argument(1, value=input_manifest)

  call OCEAN_mpi_init(ierr)
  if (ierr .ne. 0) goto 900
  mpi_started = .true.

  if (myid .eq. root) then
    write(6,'(A)') 'Standalone core-SHG pair-position vertex'
    write(6,'(A)') '  q_b,c(omega) = R_core^b X_core,c^-(omega)'
    write(6,'(A,1X,A)') '  input manifest:', trim(input_manifest)
  endif

  call OCEAN_sys_init(sys, ierr)
  if (ierr .ne. 0) goto 900

  call SHG_read_xc_manifest(trim(input_manifest), info, ierr)
  if (ierr .ne. 0) goto 900
  call SHG_validate_xc_manifest(sys, info, ierr)
  if (ierr .ne. 0) goto 900

  b_axis = info%component(2:2)
  c_axis = info%component(3:3)
  cart_axis = SHG_axis_index(b_axis)
  if (cart_axis .eq. 0) then
    ierr = 7501
    goto 900
  endif

  call OCEAN_tk_init()
  call OCEAN_tk_start(tk_total)
  timing_started = .true.

  ! Initialize only the original OCEAN vector layout and core admission mask.
  ! The BSE Hamiltonian and GMRES are not needed for this vertex-only stage.
  call OCEAN_psi_init(sys, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_energies_init(sys, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_energies_load(sys, complex_bse, ierr)
  if (ierr .ne. 0) goto 900

  call OCEAN_psi_new(xc, ierr)
  if (ierr .ne. 0) goto 900
  call OCEAN_psi_new(rbxc, ierr)
  if (ierr .ne. 0) goto 900
  have_vectors = .true.

  call SHG_load_conduction_states(sys, 'con.u2.dat', conduction_states, &
       nx_local, x_start, raw_norm_min, raw_norm_max, ierr)
  if (ierr .ne. 0) goto 900
  if (myid .eq. root) then
    write(6,'(A,2(1X,ES16.8))') 'Raw conduction-state norm range:', &
                                 raw_norm_min, raw_norm_max
  endif

  nactive = 0
  do direction = 1, 3
    weight(direction) = SHG_cartesian_weight(sys, cart_axis, direction)
    if (abs(weight(direction)) .gt. 100.0_DP*epsilon(1.0_DP)) then
      nactive = nactive + 1
      active_direction(nactive) = direction
    endif
  enddo
  if (nactive .lt. 1) then
    if (myid .eq. root) write(6,'(A)') 'No finite reduced-direction weight for the requested Cartesian vertex.'
    ierr = 7502
    goto 900
  endif

  allocate(output_norm(info%energy_count), stat=ierr)
  if (ierr .ne. 0) goto 900
  output_norm = 0.0_DP

  do iactive = 1, nactive
    direction = active_direction(iactive)
    if (myid .eq. root) then
      write(6,'(A,I0,A,ES16.8)') 'Building +delta-k links for reduced direction ', &
                                  direction, ', weight = ', weight(direction)
    endif
    call SHG_build_forward_links(sys, conduction_states, nx_local, x_start, direction, uplus, ierr)
    if (ierr .ne. 0) goto 900

    do i = 1, info%energy_count
      call read_vector_stream(xc, trim(info%vector_file(i)), ierr)
      if (ierr .ne. 0) goto 900

      call cache_names(b_axis, c_axis, i, output_file, work_file)
      if (iactive .eq. 1) then
        call OCEAN_psi_zero_full(rbxc, ierr)
      else
        call read_vector_stream(rbxc, trim(work_file), ierr)
      endif
      if (ierr .ne. 0) goto 900

      call SHG_apply_core_R_direction(sys, direction, cart_axis, uplus, xc, rbxc, ierr)
      if (ierr .ne. 0) goto 900

      if (iactive .eq. nactive) then
        ! The input was already admitted.  Reapply W_C after R_b because the
        ! neighbor transport can generate components outside the core basis.
        call OCEAN_energies_allow_full(sys, rbxc, ierr)
        if (ierr .ne. 0) goto 900
        call write_vector_stream(rbxc, trim(output_file), ierr)
        if (ierr .ne. 0) goto 900
        call full_vector_norm(sys, rbxc, output_norm(i), ierr)
        if (ierr .ne. 0) goto 900
        if (nactive .gt. 1) call delete_stream_file(trim(work_file), ierr)
      else
        call write_vector_stream(rbxc, trim(work_file), ierr)
      endif
      if (ierr .ne. 0) goto 900
    enddo

    deallocate(uplus)
  enddo

  write(output_manifest,'(A,A1,A1,A)') 'shg_rbxc_core_', b_axis, c_axis, '.info'
  call write_output_manifest(sys, info, trim(input_manifest), b_axis, c_axis, &
       weight, raw_norm_min, raw_norm_max, output_norm, trim(output_manifest), ierr)
  if (ierr .ne. 0) goto 900

  if (myid .eq. root) then
    write(6,'(A,I0,A,1X,A)') 'Saved ', info%energy_count, &
         ' physical R_b X_core,c^-(omega) vectors; manifest:', trim(output_manifest)
  endif

900 continue
  if (allocated(uplus)) deallocate(uplus)
  if (allocated(conduction_states)) deallocate(conduction_states)
  if (allocated(output_norm)) deallocate(output_norm)
  if (have_vectors) then
    call OCEAN_psi_kill(xc, cleanup_ierr)
    call OCEAN_psi_kill(rbxc, cleanup_ierr)
  endif

  if (timing_started) then
    call OCEAN_tk_stop(tk_total)
    call OCEAN_tk_printtimes(myid, comm)
  endif

  if (ierr .eq. 0 .and. cleanup_ierr .ne. 0) ierr = cleanup_ierr
  if (mpi_started) call OCEAN_mpi_finalize(mpi_ierr)
  if (ierr .eq. 0 .and. mpi_ierr .ne. 0) ierr = mpi_ierr

  if (ierr .ne. 0) then
    if (mpi_started) then
      if (myid .eq. root) write(6,'(A,I0)') 'shg_rbxc_core failed with ierr = ', ierr
    else
      write(6,'(A,I0)') 'shg_rbxc_core failed during MPI initialization, ierr = ', ierr
    endif
    stop 1
  endif

contains

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
        write(6,'(A,1X,A)') 'Missing or incompatible exciton cache:', trim(path)
        local_ierr = 7510
      else
        open(unit=92, file=trim(path), form='unformatted', access='stream', &
             status='old', iostat=ios)
        if (ios .ne. 0) then
          local_ierr = 7511
        else
          read(92,iostat=ios) buffer
          close(92)
          if (ios .ne. 0) local_ierr = 7512
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
    type(OCEAN_vector), intent(in) :: vector
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: buffer(:)
    integer :: n, ios, local_ierr, bcast_ierr

    local_ierr = 0
    if (myid .eq. root) then
      n = OCEAN_psi_size_full(vector)
      allocate(buffer(n), stat=ios)
      if (ios .ne. 0) then
        local_ierr = 7520
      else
        call OCEAN_psi_vtor(vector, buffer)
        open(unit=93, file=trim(path), form='unformatted', access='stream', &
             status='replace', iostat=ios)
        if (ios .ne. 0) then
          local_ierr = 7521
        else
          write(93,iostat=ios) buffer
          close(93)
          if (ios .ne. 0) local_ierr = 7522
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


  subroutine full_vector_norm(sys, vector, norm, ierr)
    type(o_system), intent(in) :: sys
    type(OCEAN_vector), intent(in) :: vector
    real(DP), intent(out) :: norm
    integer, intent(inout) :: ierr

    integer :: bcast_ierr

    norm = 0.0_DP
    if (myid .eq. root) then
      norm = sqrt(sum(vector%r(1:sys%cur_run%num_bands,1:sys%nkpts,1:sys%nalpha)**2 + &
                      vector%i(1:sys%cur_run%num_bands,1:sys%nkpts,1:sys%nalpha)**2))
    endif
    call MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) ierr = bcast_ierr
  end subroutine full_vector_norm


  subroutine cache_names(b_axis, c_axis, index, output_file, work_file)
    character(len=1), intent(in) :: b_axis, c_axis
    integer, intent(in) :: index
    character(len=*), intent(out) :: output_file, work_file

    write(output_file,'(A,A1,A,A1,A,I6.6,A)') &
      'shg_cache_R_', b_axis, '_Gm_core_w_', c_axis, '_', index, '.bin'
    write(work_file,'(A,A1,A,A1,A,I6.6,A)') &
      '.shg_work_R_', b_axis, '_Gm_core_w_', c_axis, '_', index, '.bin'
  end subroutine cache_names


  subroutine delete_stream_file(path, ierr)
    character(len=*), intent(in) :: path
    integer, intent(inout) :: ierr

    integer :: ios, local_ierr, bcast_ierr
    logical :: exists

    local_ierr = 0
    if (myid .eq. root) then
      inquire(file=trim(path), exist=exists)
      if (exists) then
        open(unit=94, file=trim(path), form='unformatted', access='stream', &
             status='old', iostat=ios)
        if (ios .eq. 0) then
          close(94,status='delete',iostat=ios)
        endif
        if (ios .ne. 0) local_ierr = 7530
      endif
    endif
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, bcast_ierr)
    if (bcast_ierr .ne. 0) then
      ierr = bcast_ierr
    elseif (local_ierr .ne. 0) then
      ierr = local_ierr
    endif
  end subroutine delete_stream_file


  subroutine write_output_manifest(sys, input, input_path, b_axis, c_axis, &
       weights, raw_norm_min, raw_norm_max, out_norm, path, ierr)
    use OCEAN_constants, only : PI_DP

    type(o_system), intent(in) :: sys
    type(SHG_xc_manifest), intent(in) :: input
    character(len=*), intent(in) :: input_path, path
    character(len=1), intent(in) :: b_axis, c_axis
    real(DP), intent(in) :: weights(3), raw_norm_min, raw_norm_max, out_norm(:)
    integer, intent(inout) :: ierr

    character(len=256) :: vector_file, work
    real(DP) :: theta
    integer :: i, j, ios, local_ierr, bcast_ierr

    local_ierr = 0
    if (myid .eq. root) then
      open(unit=95, file=trim(path), form='formatted', status='replace', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7540
      else
        write(95,'(A)') 'status RBXC_COMPLETE'
        write(95,'(A)') 'sector CORE_XAS'
        write(95,'(A,1X,A3)') 'component', input%component
        write(95,'(A,1X,A1)') 'middle_axis', b_axis
        write(95,'(A,1X,A1)') 'right_axis', c_axis
        write(95,'(A)') 'equation RbXc=R_core^b*Xcore_c_minus'
        write(95,'(A)') 'operator covariant_central_difference_electron_minus_core_hole'
        write(95,'(A)') 'binary_layout raw_complex_DP_stream_in_OCEAN_psi_vtor_order'
        write(95,'(A,1X,A)') 'input_manifest', trim(input_path)
        write(95,'(A,1X,ES24.16)') 'eta_hartree', input%eta_hartree
        write(95,'(A,3(1X,I0))') 'core_edge_ZNL', input%znl
        write(95,'(A,1X,A,1X,I0,1X,A,1X,I0)') 'core_site', &
             trim(input%element), input%site_index, trim(input%corelevel), input%photon
        write(95,'(A)') 'conduction_wavefunctions con.u2.dat'
        write(95,'(A)') 'conduction_mesh base'
        write(95,'(A)') 'conduction_mesh_offset_file k0.ipt'
        write(95,'(A,3(1X,ES24.16))') 'upstream_photon_q_reduced',sys%qinunitsofbvectors
        write(95,'(A)') 'conduction_photon_q_shift_reduced 0 0 0'
        write(95,'(A)') 'core_operator_q_reduced 0 0 0'
        write(95,'(A,2(1X,I0))') 'conduction_active_brange', &
          sys%brange(3),sys%brange(3)+sys%cur_run%num_bands-1
        write(95,'(A)') 'conduction_normalization unit_discrete_norm_before_links'
        write(95,'(A,2(1X,ES24.16))') 'raw_conduction_norm_range', raw_norm_min, raw_norm_max
        write(95,'(A,3(1X,I0))') 'kmesh', sys%kmesh
        write(95,'(A,3(1X,ES24.16))') 'core_tau_fractional', sys%cur_run%tau
        write(95,'(A,3(1X,ES24.16))') 'cartesian_weights_bohr', weights
        do j = 1, 3
          theta = 2.0_DP*PI_DP*sys%cur_run%tau(j)/real(sys%kmesh(j),DP)
          write(95,'(A,1X,I0,3(1X,ES24.16))') 'hole_link_direction', j, theta, cos(theta), sin(theta)
        enddo
        write(95,'(A,3(1X,I0))') 'dimensions ncon nk nalpha', &
             sys%cur_run%num_bands, sys%nkpts, sys%nalpha
        write(95,'(A,1X,I0)') 'energy_count', input%energy_count
        do i = 1, input%energy_count
          call cache_names(b_axis, c_axis, i, vector_file, work)
          write(95,'(A,1X,A1,1X,A1,1X,I0,4(1X,ES24.16),2(1X,A))') &
               'RbXc', b_axis, c_axis, i, input%omega_hartree(i), input%omega_ev(i), &
               input%input_norm(i), out_norm(i), trim(input%vector_file(i)), trim(vector_file)
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
  end subroutine write_output_manifest

end program shg_rbxc_core
