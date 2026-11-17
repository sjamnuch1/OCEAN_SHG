! Copyright (C) 2026 OCEAN collaboration
!
! Standalone core-sector pair-position vertex for SHG.
!
! This module is deliberately kept outside OCEAN2.  It consumes public OCEAN
! types and the unmodified PREP conduction wavefunctions staged as con.u2.dat.
module SHG_core_vertex
  use AI_kinds, only : DP

  implicit none
  private

  integer, parameter :: path_len = 256

  type, public :: SHG_xc_manifest
    character(len=32) :: status = ''
    character(len=32) :: sector = ''
    character(len=3) :: component = '???'
    character(len=1) :: right_axis = '?'
    character(len=2) :: element = '__'
    character(len=2) :: corelevel = '__'
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
    real(DP), allocatable :: residual(:)
    real(DP), allocatable :: physical_residual(:)
    real(DP), allocatable :: input_norm(:)
    integer, allocatable :: iterations(:)
    character(len=path_len), allocatable :: vector_file(:)
  end type SHG_xc_manifest

  public :: SHG_read_xc_manifest, SHG_validate_xc_manifest
  public :: SHG_load_conduction_states, SHG_build_forward_links
  public :: SHG_apply_core_R_direction, SHG_axis_index
  public :: SHG_cartesian_weight, SHG_neighbor_k

contains

  subroutine SHG_read_xc_manifest(path, info, ierr)
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_CHARACTER, &
                          MPI_DOUBLE_PRECISION, MPI_INTEGER

    character(len=*), intent(in) :: path
    type(SHG_xc_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr

    integer :: local_ierr, mpi_ierr, n

    local_ierr = 0
    if (myid .eq. root) call read_xc_manifest_root(path, info, local_ierr)

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
    call MPI_BCAST(info%right_axis, 1, MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%element, len(info%element), MPI_CHARACTER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%corelevel, len(info%corelevel), MPI_CHARACTER, root, comm, ierr)
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
    call MPI_BCAST(info%residual, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%physical_residual, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%input_norm, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%iterations, n, MPI_INTEGER, root, comm, ierr)
    if (ierr .ne. 0) return
    call MPI_BCAST(info%vector_file, path_len*n, MPI_CHARACTER, root, comm, ierr)
  end subroutine SHG_read_xc_manifest


  subroutine read_xc_manifest_root(path, info, ierr)
    character(len=*), intent(in) :: path
    type(SHG_xc_manifest), intent(inout) :: info
    integer, intent(inout) :: ierr

    character(len=1024) :: line
    character(len=64) :: key, word2, word3, word4
    character(len=path_len) :: temp_file
    character(len=1) :: axis
    integer :: ios, idx, count
    integer :: temp_iterations
    real(DP) :: temp_omega_ha, temp_omega_ev, temp_residual
    real(DP) :: temp_physical_residual, temp_input_norm

    open(unit=91, file=trim(path), form='formatted', status='old', iostat=ios)
    if (ios .ne. 0) then
      write(6,'(A,1X,A)') 'Failed to open first-layer manifest:', trim(path)
      ierr = 7401
      return
    endif

    count = 0
    do
      read(91,'(A)',iostat=ios) line
      if (ios .lt. 0) exit
      if (ios .ne. 0) then
        ierr = 7402
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
        case('right_axis')
          read(line,*,iostat=ios) key, info%right_axis
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
        case('Gminus_omega')
          if (.not. allocated(info%omega_hartree)) then
            ios = 1
          else
            read(line,*,iostat=ios) key, axis, idx, &
                 temp_omega_ha, temp_omega_ev, temp_residual, &
                 temp_physical_residual, temp_input_norm, temp_iterations, temp_file
            if (ios .eq. 0) then
              if (idx .lt. 1 .or. idx .gt. info%energy_count) then
                ios = 1
              elseif (axis .ne. info%right_axis) then
                ios = 1
              else
                info%omega_hartree(idx) = temp_omega_ha
                info%omega_ev(idx) = temp_omega_ev
                info%residual(idx) = temp_residual
                info%physical_residual(idx) = temp_physical_residual
                info%input_norm(idx) = temp_input_norm
                info%iterations(idx) = temp_iterations
                info%vector_file(idx) = temp_file
                count = count + 1
              endif
            endif
          endif
        case default
          ! The first-layer manifest contains additional diagnostic fields.
          ! They are deliberately ignored here.
          ios = 0
      end select

      if (ios .ne. 0) then
        write(6,'(A,1X,A)') 'Malformed line in first-layer manifest:', trim(line)
        ierr = 7403
        exit
      endif
    enddo
    close(91)

    if (ierr .ne. 0) return
    if (trim(info%status) .ne. 'XC_MINUS_COMPLETE' .or. &
        trim(info%sector) .ne. 'CORE_XAS') then
      write(6,'(A)') 'The input manifest is not a completed core Xc-minus calculation.'
      ierr = 7404
      return
    endif
    if (count .ne. info%energy_count) then
      write(6,'(A,2(1X,I0))') 'Gminus entry count mismatch:', count, info%energy_count
      ierr = 7405
    endif
  end subroutine read_xc_manifest_root


  subroutine allocate_manifest_arrays(info, n, ierr)
    type(SHG_xc_manifest), intent(inout) :: info
    integer, intent(in) :: n
    integer, intent(inout) :: ierr

    if (allocated(info%omega_hartree)) deallocate(info%omega_hartree)
    if (allocated(info%omega_ev)) deallocate(info%omega_ev)
    if (allocated(info%residual)) deallocate(info%residual)
    if (allocated(info%physical_residual)) deallocate(info%physical_residual)
    if (allocated(info%input_norm)) deallocate(info%input_norm)
    if (allocated(info%iterations)) deallocate(info%iterations)
    if (allocated(info%vector_file)) deallocate(info%vector_file)

    allocate(info%omega_hartree(n), info%omega_ev(n), info%residual(n), &
             info%physical_residual(n), info%input_norm(n), &
             info%iterations(n), info%vector_file(n), stat=ierr)
    if (ierr .ne. 0) return

    info%omega_hartree = 0.0_DP
    info%omega_ev = 0.0_DP
    info%residual = 0.0_DP
    info%physical_residual = 0.0_DP
    info%input_norm = 0.0_DP
    info%iterations = 0
    info%vector_file = ''
  end subroutine allocate_manifest_arrays


  subroutine SHG_validate_xc_manifest(sys, info, ierr)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use OCEAN_mpi, only : myid, root
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    type(SHG_xc_manifest), intent(in) :: info
    integer, intent(inout) :: ierr

    integer :: j

    if (.not. associated(sys%cur_run)) then
      ierr = 7410
    elseif (sys%nruns .ne. 1) then
      ierr = 7411
    elseif (.not. sys%cur_run%have_core .or. sys%cur_run%calc_type .ne. 'XAS') then
      ierr = 7412
    elseif (any(info%kmesh .ne. sys%kmesh)) then
      ierr = 7413
    elseif (info%nk .ne. sys%nkpts .or. info%ncon .ne. sys%cur_run%num_bands .or. &
            info%nalpha .ne. sys%nalpha) then
      ierr = 7414
    elseif (any(info%znl .ne. sys%cur_run%ZNL)) then
      ierr = 7415
    elseif (trim(info%element) .ne. trim(sys%cur_run%elname) .or. &
            info%site_index .ne. sys%cur_run%indx .or. &
            trim(info%corelevel) .ne. trim(sys%cur_run%corelevel) .or. &
            info%photon .ne. sys%cur_run%photon) then
      ierr = 7416
    elseif (info%right_axis .ne. info%component(3:3)) then
      ierr = 7417
    elseif (.not. all(ieee_is_finite(sys%qinunitsofbvectors))) then
      ierr = 7419
    elseif (sys%nspn .ne. 1 .and. sys%nspn .ne. 2) then
      ierr = 7423
    endif

    if (ierr .eq. 0) then
      if (.not. allocated(info%vector_file)) then
        ierr = 7424
      else
        do j = 1, size(info%vector_file)
          if (len_trim(info%vector_file(j)) .eq. 0) ierr = 7424
        enddo
      endif
    endif

    do j = 1, 3
      if (SHG_axis_index(info%component(j:j)) .eq. 0) ierr = 7418
    enddo

    if (ierr .ne. 0 .and. myid .eq. root) then
      write(6,'(A,I0)') 'First-layer manifest does not match the staged XAS run, ierr = ', ierr
    endif
    ! OCEAN QE/PREP keeps conduction states on the base mesh even when the
    ! occupied states use k-q. sys%kshift describes the upstream photon-q
    ! metadata; it is not a shift selector for the con.u2.dat stream.
    ! Preserve that metadata and all basis checks above. Both core vertices
    ! continue to use the base-mesh conduction basis of the incoming core X.
    if (ierr .eq. 0 .and. myid .eq. root .and. sys%kshift) then
      write(6,'(A,3(1X,ES16.8))') 'Upstream photon q (reduced):',sys%qinunitsofbvectors
      write(6,'(A)') 'Using base-mesh con.u2.dat for the q=0 core position operator.'
      write(6,'(A)') 'The incoming core X and core projectors must use this same conduction basis.'
    endif
  end subroutine SHG_validate_xc_manifest


  subroutine SHG_load_conduction_states(sys, path, states, nx_local, x_start, &
                                         norm_min, norm_max, ierr)
    use OCEAN_mpi
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    character(len=*), intent(in) :: path
    complex(DP), allocatable, intent(out) :: states(:,:,:,:)
    integer, intent(out) :: nx_local, x_start
    real(DP), intent(out) :: norm_min, norm_max
    integer, intent(inout) :: ierr

    real(DP), allocatable :: norms(:,:,:)
    integer :: bands_in_file, nx_left, rank, ib, ik, ispn, fh
    integer :: mpi_ierr, local_ierr, nbytes_complex
    integer :: status_array(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: displacement, offset
    integer(8) :: expected_bytes, actual_bytes

    bands_in_file = sys%brange(4) - sys%brange(3) + 1
    if (sys%num_bands .gt. bands_in_file) then
      if (myid .eq. root) write(6,'(A)') 'con.u2.dat contains fewer bands than the core basis.'
      ierr = 7420
      return
    endif

    nbytes_complex = storage_size(cmplx(0.0_DP,0.0_DP,DP))/8
    expected_bytes = int(sys%nxpts,8) * int(bands_in_file,8) * &
                     int(sys%nkpts,8) * int(sys%nspn,8) * int(nbytes_complex,8)
    actual_bytes = -1_8
    local_ierr = 0
    if (myid .eq. root) then
      inquire(file=trim(path), size=actual_bytes, iostat=local_ierr)
      if (local_ierr .eq. 0 .and. actual_bytes .ne. expected_bytes) local_ierr = 7421
    endif
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, mpi_ierr)
    if (mpi_ierr .ne. 0) then
      ierr = mpi_ierr
      return
    elseif (local_ierr .ne. 0) then
      if (myid .eq. root) then
        write(6,'(A,2(1X,I0))') 'Unexpected con.u2.dat byte count:', actual_bytes, expected_bytes
      endif
      ierr = local_ierr
      return
    endif

    nx_left = sys%nxpts
    x_start = 1
    nx_local = 0
    do rank = 0, myid
      x_start = x_start + nx_local
      nx_local = nx_left/(nproc-rank)
      nx_left = nx_left - nx_local
    enddo

    allocate(states(max(1,nx_local),sys%num_bands,sys%nkpts,sys%nspn), &
             norms(sys%num_bands,sys%nkpts,sys%nspn), stat=ierr)
    if (ierr .ne. 0) return
    states = cmplx(0.0_DP,0.0_DP,DP)
    norms = 0.0_DP

    call MPI_FILE_OPEN(comm, trim(path), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if (ierr .ne. 0) return
    displacement = int(0,MPI_OFFSET_KIND)
    call MPI_FILE_SET_VIEW(fh, displacement, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, &
                           'native', MPI_INFO_NULL, ierr)
    if (ierr .ne. 0) then
      call MPI_FILE_CLOSE(fh, mpi_ierr)
      return
    endif

    do ispn = 1, sys%nspn
      do ik = 1, sys%nkpts
        do ib = 1, sys%num_bands
          offset = int((((ispn-1)*sys%nkpts + (ik-1))*bands_in_file + (ib-1)), &
                       MPI_OFFSET_KIND) * int(sys%nxpts,MPI_OFFSET_KIND) &
                 + int(x_start-1,MPI_OFFSET_KIND)
          call MPI_FILE_READ_AT_ALL(fh, offset, states(1,ib,ik,ispn), nx_local, &
                                    MPI_DOUBLE_COMPLEX, status_array, ierr)
          if (ierr .ne. 0) then
            call MPI_FILE_CLOSE(fh, mpi_ierr)
            return
          endif
        enddo
      enddo
    enddo
    call MPI_FILE_CLOSE(fh, ierr)
    if (ierr .ne. 0) return

    if (nx_local .gt. 0) then
      do ispn = 1, sys%nspn
        do ik = 1, sys%nkpts
          do ib = 1, sys%num_bands
            norms(ib,ik,ispn) = sum(real(conjg(states(1:nx_local,ib,ik,ispn)) * &
                                               states(1:nx_local,ib,ik,ispn),DP))
          enddo
        enddo
      enddo
    endif
    call MPI_ALLREDUCE(MPI_IN_PLACE, norms, size(norms), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, comm, ierr)
    if (ierr .ne. 0) return

    norm_min = minval(norms)
    norm_max = maxval(norms)
    if (norm_min .le. tiny(1.0_DP)) then
      if (myid .eq. root) write(6,'(A)') 'A staged conduction state has zero norm.'
      ierr = 7422
      return
    endif

    if (nx_local .gt. 0) then
      do ispn = 1, sys%nspn
        do ik = 1, sys%nkpts
          do ib = 1, sys%num_bands
            states(1:nx_local,ib,ik,ispn) = states(1:nx_local,ib,ik,ispn) &
                                           / sqrt(norms(ib,ik,ispn))
          enddo
        enddo
      enddo
    endif
    deallocate(norms)
  end subroutine SHG_load_conduction_states


  subroutine SHG_build_forward_links(sys, states, nx_local, x_start, direction, uplus, ierr)
    use OCEAN_constants, only : PI_DP
    use OCEAN_mpi
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    complex(DP), intent(in) :: states(:,:,:,:)
    integer, intent(in) :: nx_local, x_start, direction
    complex(DP), allocatable, intent(out) :: uplus(:,:,:,:)
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: local_link(:,:), neighbor_states(:,:)
    complex(DP) :: boundary_phase
    real(DP) :: fraction
    integer :: ik, ikp, ispn, ncon, remainder, ix_local, ix_global
    integer :: kindex(3), xindex(3)

    ncon = sys%cur_run%num_bands
    allocate(uplus(ncon,ncon,sys%nkpts,sys%nspn), local_link(ncon,ncon), &
             neighbor_states(max(1,nx_local),ncon), stat=ierr)
    if (ierr .ne. 0) return
    uplus = cmplx(0.0_DP,0.0_DP,DP)

    do ispn = 1, sys%nspn
      do ik = 1, sys%nkpts
        ikp = SHG_neighbor_k(ik, direction, 1, sys%kmesh)
        local_link = cmplx(0.0_DP,0.0_DP,DP)
        if (nx_local .gt. 0) then
          remainder = ik - 1
          kindex(1) = remainder/(sys%kmesh(2)*sys%kmesh(3)) + 1
          remainder = mod(remainder,sys%kmesh(2)*sys%kmesh(3))
          kindex(2) = remainder/sys%kmesh(3) + 1
          kindex(3) = mod(remainder,sys%kmesh(3)) + 1

          if (kindex(direction) .eq. sys%kmesh(direction)) then
            ! Across a Brillouin-zone face, the full-Bloch transport carries
            ! exp(-i G_j.r).  con.u2.dat stores cell-periodic u, so restore
            ! this one reciprocal-lattice boundary phase explicitly.
            do ix_local = 1, nx_local
              ix_global = x_start + ix_local - 1
              remainder = ix_global - 1
              xindex(1) = mod(remainder,sys%xmesh(1))
              xindex(2) = mod(remainder/sys%xmesh(1),sys%xmesh(2))
              xindex(3) = remainder/(sys%xmesh(1)*sys%xmesh(2))
              fraction = real(xindex(direction),DP)/real(sys%xmesh(direction),DP)
              boundary_phase = cmplx(cos(-2.0_DP*PI_DP*fraction), &
                                     sin(-2.0_DP*PI_DP*fraction),DP)
              neighbor_states(ix_local,:) = boundary_phase*states(ix_local,:,ikp,ispn)
            enddo
            local_link = matmul(conjg(transpose(states(1:nx_local,:,ik,ispn))), &
                                neighbor_states(1:nx_local,:))
          else
            local_link = matmul(conjg(transpose(states(1:nx_local,:,ik,ispn))), &
                                states(1:nx_local,:,ikp,ispn))
          endif
        endif
        call MPI_ALLREDUCE(local_link, uplus(:,:,ik,ispn), ncon*ncon, &
                           MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
        if (ierr .ne. 0) return
      enddo
    enddo

    deallocate(local_link,neighbor_states)
  end subroutine SHG_build_forward_links


  subroutine SHG_apply_core_R_direction(sys, direction, cart_axis, uplus, x, y, ierr)
    use OCEAN_constants, only : PI_DP
    use OCEAN_psi, only : OCEAN_vector
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    integer, intent(in) :: direction, cart_axis
    complex(DP), intent(in) :: uplus(:,:,:,:)
    type(OCEAN_vector), intent(in) :: x
    type(OCEAN_vector), intent(inout) :: y
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: xplus(:), xminus(:), term_plus(:), term_minus(:)
    complex(DP) :: phase_plus, phase_minus, iw
    real(DP) :: theta, weight
    integer :: ialpha, ik, ikp, ikm, ispn, ncon

    ncon = sys%cur_run%num_bands
    if (sys%nspn .ne. 1 .and. sys%nspn .ne. 2) then
      ierr = 7430
      return
    endif

    weight = SHG_cartesian_weight(sys, cart_axis, direction)
    if (abs(weight) .le. 100.0_DP*epsilon(1.0_DP)) return

    theta = 2.0_DP*PI_DP*sys%cur_run%tau(direction) / real(sys%kmesh(direction),DP)
    phase_plus = cmplx(cos(theta),sin(theta),DP)
    phase_minus = conjg(phase_plus)
    iw = cmplx(0.0_DP,weight,DP)

    allocate(xplus(ncon), xminus(ncon), term_plus(ncon), term_minus(ncon), stat=ierr)
    if (ierr .ne. 0) return

    do ialpha = 1, sys%nalpha
      if (sys%nspn .eq. 1) then
        ispn = 1
      else
        ! This is the same core-alpha to electron-spin map used by original
        ! OCEAN_multiplet and OCEAN_energies.
        ispn = 2 - mod(ialpha,2)
      endif

      do ik = 1, sys%nkpts
        ikp = SHG_neighbor_k(ik, direction, 1, sys%kmesh)
        ikm = SHG_neighbor_k(ik, direction, -1, sys%kmesh)

        xplus = cmplx(x%r(1:ncon,ikp,ialpha), x%i(1:ncon,ikp,ialpha), DP)
        xminus = cmplx(x%r(1:ncon,ikm,ialpha), x%i(1:ncon,ikm,ialpha), DP)

        term_plus = matmul(uplus(:,:,ik,ispn), xplus)
        ! U^-(k) = [U^+(k-delta_k)]^dagger.  Constructing the backward
        ! transport this way enforces the adjoint relation exactly.
        term_minus = matmul(conjg(transpose(uplus(:,:,ikm,ispn))), xminus)

        term_plus = iw*(phase_plus*term_plus - phase_minus*term_minus)
        y%r(1:ncon,ik,ialpha) = y%r(1:ncon,ik,ialpha) + real(term_plus,DP)
        y%i(1:ncon,ik,ialpha) = y%i(1:ncon,ik,ialpha) + aimag(term_plus)
      enddo
    enddo

    deallocate(xplus, xminus, term_plus, term_minus)
  end subroutine SHG_apply_core_R_direction


  integer function SHG_neighbor_k(ik, direction, step, kmesh) result(neighbor)
    integer, intent(in) :: ik, direction, step, kmesh(3)

    integer :: ij(3), remainder

    remainder = ik - 1
    ij(1) = remainder/(kmesh(2)*kmesh(3)) + 1
    remainder = mod(remainder,kmesh(2)*kmesh(3))
    ij(2) = remainder/kmesh(3) + 1
    ij(3) = mod(remainder,kmesh(3)) + 1

    ij(direction) = modulo(ij(direction)-1+step,kmesh(direction)) + 1
    neighbor = ((ij(1)-1)*kmesh(2) + (ij(2)-1))*kmesh(3) + ij(3)
  end function SHG_neighbor_k


  real(DP) function SHG_cartesian_weight(sys, cart_axis, direction) result(weight)
    use OCEAN_constants, only : PI_DP
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    integer, intent(in) :: cart_axis, direction

    ! w_j^mu = (N_j/2)(B^{-1} e_mu)_j.  Since a_j dot b_l =
    ! 2*pi*delta_jl, (B^{-1}e_mu)_j = a_j^mu/(2*pi).
    weight = real(sys%kmesh(direction),DP)*sys%avec(cart_axis,direction) &
             /(4.0_DP*PI_DP)
  end function SHG_cartesian_weight


  integer function SHG_axis_index(axis) result(index)
    character(len=1), intent(in) :: axis

    select case(axis)
      case('x','X')
        index = 1
      case('y','Y')
        index = 2
      case('z','Z')
        index = 3
      case default
        index = 0
    end select
  end function SHG_axis_index

end module SHG_core_vertex
