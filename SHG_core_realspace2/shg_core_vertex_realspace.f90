! Copyright (C) 2026 OCEAN collaboration
!
! Standalone real-space representation of the core-exciton pair-position
! vertex.  This is the Fourier-dual implementation of
!
!   R_core^mu X = i sum_j w_j^mu ( T_j^+ - T_j^- ) X .
!
! The primitive-cell functions in con.u2.dat are contracted with X, phased
! into full Bloch functions, transformed from the product k mesh to the
! Born-von Karman cell mesh, multiplied by the periodic sine coordinate, and
! projected back into the retained conduction subspace.
module SHG_core_vertex_realspace
  use AI_kinds, only : DP
  use FFT_wrapper, only : fft_obj, FFT_wrapper_init, FFT_wrapper_single, &
                          FFT_wrapper_delete, OCEAN_FORWARD, OCEAN_BACKWARD

  implicit none
  private

  type, public :: SHG_realspace_workspace
    private
    type(fft_obj) :: plan
    complex(DP), allocatable :: scratch(:)
    complex(DP), allocatable :: bloch_phase(:,:)
    real(DP), allocatable :: rho_fd(:,:)
    real(DP) :: k0(3) = 0.0_DP
    integer :: fft_dims(3) = 0
    integer :: nx_local = 0
    integer :: x_start = 0
    integer :: cart_axis = 0
    logical :: plan_ready = .false.
    logical :: initialized = .false.
  end type SHG_realspace_workspace

  public :: SHG_realspace_init, SHG_apply_core_R_realspace
  public :: SHG_realspace_clean, SHG_realspace_get_k0

contains

  pure function SHG_realspace_get_k0(work) result(k0)
    type(SHG_realspace_workspace), intent(in) :: work
    real(DP) :: k0(3)

    ! Return a copy for provenance without exposing the private workspace.
    k0 = work%k0
  end function SHG_realspace_get_k0

  subroutine SHG_realspace_init(sys, nx_local, x_start, cart_axis, work, ierr)
    use OCEAN_constants, only : PI_DP
    use OCEAN_mpi, only : myid, root, comm, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_INTEGER
    use OCEAN_system, only : o_system
    use SHG_core_vertex, only : SHG_cartesian_weight

    type(o_system), intent(in) :: sys
    integer, intent(in) :: nx_local, x_start, cart_axis
    type(SHG_realspace_workspace), intent(inout) :: work
    integer, intent(inout) :: ierr

    real(DP) :: xi(3), qred(3), angle, rho, weight(3)
    integer :: ix_local, ix_global, ik, ir, ios, local_ierr, mpi_ierr
    integer :: kindex(3), rindex(3)

    if (ierr .ne. 0) return
    if (work%initialized .or. work%plan_ready) call SHG_realspace_clean(work)

    if (cart_axis .lt. 1 .or. cart_axis .gt. 3) then
      ierr = 7601
      return
    elseif (nx_local .lt. 0 .or. x_start .lt. 1) then
      ierr = 7602
      return
    elseif (product(sys%kmesh) .ne. sys%nkpts .or. product(sys%xmesh) .ne. sys%nxpts) then
      ierr = 7603
      return
    endif

    do ir = 1, 3
      weight(ir) = SHG_cartesian_weight(sys, cart_axis, ir)
      if (abs(weight(ir)) .gt. 100.0_DP*epsilon(1.0_DP) .and. sys%kmesh(ir) .lt. 3) then
        if (myid .eq. root) then
          write(6,'(A,I0,A,I0)') 'Real-space central vertex needs kmesh >= 3 in direction ', &
                                  ir, '; found ', sys%kmesh(ir)
        endif
        ierr = 7604
        return
      endif
    enddo

    local_ierr = 0
    work%k0 = 0.0_DP
    if (myid .eq. root) then
      open(unit=89, file='k0.ipt', form='formatted', status='old', iostat=ios)
      if (ios .ne. 0) then
        local_ierr = 7605
      else
        read(89,*,iostat=ios) work%k0
        close(89)
        if (ios .ne. 0) local_ierr = 7606
      endif
    endif
    call MPI_BCAST(local_ierr, 1, MPI_INTEGER, root, comm, mpi_ierr)
    if (mpi_ierr .ne. 0) then
      ierr = mpi_ierr
      return
    elseif (local_ierr .ne. 0) then
      if (myid .eq. root) write(6,'(A)') 'Failed to read the three reduced k offsets from k0.ipt.'
      ierr = local_ierr
      return
    endif
    call MPI_BCAST(work%k0, 3, MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr .ne. 0) return

    work%nx_local = nx_local
    work%x_start = x_start
    work%cart_axis = cart_axis
    ! The OCEAN vector index has k3 fastest, then k2, then k1.
    work%fft_dims = (/ sys%kmesh(3), sys%kmesh(2), sys%kmesh(1) /)

    allocate(work%scratch(sys%nkpts), &
             work%bloch_phase(max(1,nx_local),sys%nkpts), &
             work%rho_fd(max(1,nx_local),sys%nkpts), stat=ierr)
    if (ierr .ne. 0) return
    work%scratch = cmplx(0.0_DP,0.0_DP,DP)
    work%bloch_phase = cmplx(1.0_DP,0.0_DP,DP)
    work%rho_fd = 0.0_DP

    ! OCEAN_FORWARD is the negative-exponent transform.  The k -> R synthesis
    ! below therefore uses OCEAN_BACKWARD, whose default 1/Nk normalization is
    ! exactly undone by the subsequent unnormalized OCEAN_FORWARD transform.
    call FFT_wrapper_init(work%fft_dims, work%plan, work%scratch)
    work%plan_ready = .true.

    do ix_local = 1, nx_local
      ix_global = x_start + ix_local - 1
      call file_grid_fraction(ix_global, sys%xmesh, xi)

      do ik = 1, sys%nkpts
        call mesh_index(ik, sys%kmesh, kindex)
        ! con.u2.dat is on OCEAN's base conduction mesh. The upstream photon
        ! q shifts the occupied mesh by -q, not this conduction mesh. Retain
        ! the ordinary k0 offset and do not add photon q to the Bloch phase.
        qred = (work%k0 + real(kindex,DP))/real(sys%kmesh,DP)
        angle = 2.0_DP*PI_DP*dot_product(qred,xi)
        work%bloch_phase(ix_local,ik) = cmplx(cos(angle),sin(angle),DP)
      enddo

      do ir = 1, sys%nkpts
        call mesh_index(ir, sys%kmesh, rindex)
        rho = 0.0_DP
        do ik = 1, 3
          if (abs(weight(ik)) .le. 100.0_DP*epsilon(1.0_DP)) cycle
          angle = 2.0_DP*PI_DP * &
                  (real(rindex(ik),DP) + xi(ik) - sys%cur_run%tau(ik)) &
                  / real(sys%kmesh(ik),DP)
          rho = rho + 2.0_DP*weight(ik)*sin(angle)
        enddo
        work%rho_fd(ix_local,ir) = rho
      enddo
    enddo

    work%initialized = .true.
  end subroutine SHG_realspace_init


  subroutine SHG_apply_core_R_realspace(sys, states, nx_local, work, x, y, ierr)
    use OCEAN_mpi, only : comm, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_zero_full
    use OCEAN_system, only : o_system

    type(o_system), intent(in) :: sys
    complex(DP), intent(in) :: states(:,:,:,:)
    integer, intent(in) :: nx_local
    type(SHG_realspace_workspace), intent(inout) :: work
    type(OCEAN_vector), intent(in) :: x
    type(OCEAN_vector), intent(inout) :: y
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: local_result(:,:,:), x_band(:)
    complex(DP) :: phase, field
    integer :: ialpha, ispn, ix_local, ik, ic, ncon

    if (ierr .ne. 0) return
    if (.not. work%initialized .or. .not. work%plan_ready) then
      ierr = 7610
      return
    elseif (nx_local .ne. work%nx_local) then
      ierr = 7611
      return
    elseif (size(states,1) .lt. max(1,nx_local) .or. &
            size(states,2) .lt. sys%cur_run%num_bands .or. &
            size(states,3) .lt. sys%nkpts .or. size(states,4) .lt. sys%nspn) then
      ierr = 7612
      return
    endif

    ncon = sys%cur_run%num_bands
    allocate(local_result(ncon,sys%nkpts,sys%nalpha), x_band(ncon), stat=ierr)
    if (ierr .ne. 0) return
    local_result = cmplx(0.0_DP,0.0_DP,DP)

    do ialpha = 1, sys%nalpha
      if (sys%nspn .eq. 1) then
        ispn = 1
      else
        ! OCEAN alpha order has the conduction spin as the innermost index.
        ispn = mod(ialpha-1,sys%nspn) + 1
      endif

      do ix_local = 1, nx_local
        ! Band synthesis at fixed primitive-cell point:
        ! Phi(k,x) = exp(i k.x) sum_c u(c,k,x) X(c,k,alpha).
        do ik = 1, sys%nkpts
          x_band = cmplx(x%r(1:ncon,ik,ialpha),x%i(1:ncon,ik,ialpha),DP)
          field = sum(states(ix_local,1:ncon,ik,ispn)*x_band)
          work%scratch(ik) = work%bloch_phase(ix_local,ik)*field
        enddo

        call FFT_wrapper_single(work%scratch, OCEAN_BACKWARD, work%plan)
        work%scratch = work%rho_fd(ix_local,:)*work%scratch
        call FFT_wrapper_single(work%scratch, OCEAN_FORWARD, work%plan)

        ! Projection at fixed output k.  The discrete normalization of the
        ! states was established once by SHG_load_conduction_states, so no
        ! additional cell-volume or 1/Nx factor is inserted here.
        do ik = 1, sys%nkpts
          phase = work%bloch_phase(ix_local,ik)
          do ic = 1, ncon
            local_result(ic,ik,ialpha) = local_result(ic,ik,ialpha) + &
                 conjg(phase*states(ix_local,ic,ik,ispn))*work%scratch(ik)
          enddo
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, local_result, size(local_result), &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
    if (ierr .ne. 0) then
      deallocate(local_result,x_band)
      return
    endif

    call OCEAN_psi_zero_full(y, ierr)
    if (ierr .eq. 0) then
      y%r(1:ncon,1:sys%nkpts,1:sys%nalpha) = real(local_result,DP)
      y%i(1:ncon,1:sys%nkpts,1:sys%nalpha) = aimag(local_result)
    endif
    deallocate(local_result,x_band)
  end subroutine SHG_apply_core_R_realspace


  subroutine SHG_realspace_clean(work)
    type(SHG_realspace_workspace), intent(inout) :: work

    if (work%plan_ready) call FFT_wrapper_delete(work%plan)
    if (allocated(work%scratch)) deallocate(work%scratch)
    if (allocated(work%bloch_phase)) deallocate(work%bloch_phase)
    if (allocated(work%rho_fd)) deallocate(work%rho_fd)
    work%plan_ready = .false.
    work%initialized = .false.
    work%fft_dims = 0
    work%nx_local = 0
    work%x_start = 0
    work%cart_axis = 0
    work%k0 = 0.0_DP
  end subroutine SHG_realspace_clean


  pure subroutine mesh_index(linear_index, mesh, index3)
    integer, intent(in) :: linear_index, mesh(3)
    integer, intent(out) :: index3(3)

    integer :: remainder

    remainder = linear_index - 1
    index3(1) = remainder/(mesh(2)*mesh(3))
    remainder = mod(remainder,mesh(2)*mesh(3))
    index3(2) = remainder/mesh(3)
    index3(3) = mod(remainder,mesh(3))
  end subroutine mesh_index


  pure subroutine file_grid_fraction(linear_index, xmesh, fraction)
    integer, intent(in) :: linear_index, xmesh(3)
    real(DP), intent(out) :: fraction(3)

    integer :: zero_index, grid_index(3)

    ! New con.u2.dat stores the primitive grid with x fastest, then y, then z.
    zero_index = linear_index - 1
    grid_index(1) = mod(zero_index,xmesh(1))
    grid_index(2) = mod(zero_index/xmesh(1),xmesh(2))
    grid_index(3) = zero_index/(xmesh(1)*xmesh(2))
    fraction = real(grid_index,DP)/real(xmesh,DP)
  end subroutine file_grid_fraction

end module SHG_core_vertex_realspace
