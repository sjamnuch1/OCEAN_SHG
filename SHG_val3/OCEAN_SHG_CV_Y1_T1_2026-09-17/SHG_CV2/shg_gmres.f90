! Copyright (C) 2017 OCEAN collaboration
! Standalone SHG adaptation Copyright (C) 2026 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!> @brief Private standalone copy of OCEAN's recycled GMRES machinery.
!!
!! This module carries the arbitrary-RHS, signed-resolvent, and converged-vector
!! extensions required by the SHG test.  Keeping it under SHG means that the
!! OCEAN_gmres module in an upstream OCEAN checkout remains byte-for-byte
!! unchanged.
module SHG_gmres
  use OCEAN_psi, only : ocean_vector
  use AI_kinds

  implicit none
  private
  save
    
  complex(DP), allocatable :: c_matrix( : )

  type(ocean_vector), allocatable, target :: au_matrix(:)
  type(ocean_vector), allocatable, target :: u_matrix(:)

  real(DP) :: gmres_resolution !> imaginary part added to inverse, stored in Ha
  real(DP) :: gmres_preconditioner !> width of preconditioner 
  real(DP) :: gmres_convergence 
  integer :: gmres_depth ! max size of gmres Arnoldi? space before restart

  integer :: gmres_nsteps
  real(DP), allocatable :: gmres_energy_list(:)

  logical :: echamp
  logical :: do_precondition = .true.
  logical :: allow_reuse_x = .true.
  logical :: force_reuse_x = .true.
  logical :: allow_recycle = .true.

  !SJ BEGIN expose the configured linear-system parameters to the SHG preflight
  public :: SHG_gmres_do, SHG_gmres_setup, SHG_gmres_clean, SHG_gmres_do_recycle, &
            SHG_gmres_get_parameters
  !SJ END expose the configured linear-system parameters to the SHG preflight

  contains

  subroutine SHG_gmres_clean( ierr )
    use OCEAN_psi, only : OCEAN_psi_kill
    !
    integer, intent( inout ) :: ierr
    !
    integer :: i, n
    !
    if( allocated( c_matrix ) ) deallocate( c_matrix ) 
    if( allocated( gmres_energy_list ) ) deallocate( gmres_energy_list ) 

    if( allocated( au_matrix ) ) then
      n = size( au_matrix, 1 )
      do i = 1, n
        call OCEAN_psi_kill( au_matrix( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
      deallocate( au_matrix )
    endif

    if( allocated( u_matrix ) ) then
      n = size( u_matrix, 1 )
      do i = 1, n
        call OCEAN_psi_kill( u_matrix( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
      deallocate( u_matrix )
    endif

  end subroutine SHG_gmres_clean

  subroutine SHG_gmres_setup( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi!, only : myid, root, comm, MPI_BCAST, MPI_INTEGER, MPI_DOUBLE_PRECISION
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_constants, only : eV2Hartree
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    real(dp) :: e_start, e_stop, e_step
    character(len=4) :: inv_style
    integer :: ierr_, i
    !
    if( myid .eq. root ) then
!      open(unit=99,file='gmres.in',form='formatted',status='old',iostat=ierr,ERR=100)
      open(unit=99,file='bse.in',form='formatted',status='old',iostat=ierr,ERR=100)
      read(99,*)
      read(99,*)
      read(99,*)

      read(99,*) gmres_depth, gmres_resolution, gmres_preconditioner, gmres_convergence
      !
      ! if gres is negative fill it with core-hole lifetime broadening
      if( gmres_resolution .lt. 1.0d-20 ) then
        call returnLifetime( sys%cur_run%ZNL(1), sys%cur_run%ZNL(2), sys%cur_run%ZNL(3), gmres_resolution )
        if( gmres_resolution .le. 0.0_dp ) gmres_resolution = 0.1_dp
      endif
      gmres_resolution = gmres_resolution * eV2Hartree
      gmres_preconditioner = gmres_preconditioner * eV2Hartree
      !
      read(99,*) inv_style
      select case( inv_style )

        case('list')
          read(99,*) gmres_nsteps
          allocate( gmres_energy_list( gmres_nsteps ) )
          do i = 1, gmres_nsteps
            read(99,*) gmres_energy_list( i )
          enddo

        case('loop')
          read(99,*) e_start, e_stop, e_step
          gmres_nsteps = 1 + ceiling( ( e_stop - e_start - e_step*0.0001_DP ) / e_step )
          allocate( gmres_energy_list( gmres_nsteps ) )
          do i = 1, gmres_nsteps
            gmres_energy_list( i ) = e_start + ( i - 1 ) * e_step
          enddo

        case default
          write(6,*) 'Error in gmres.in! Unsupported energy point style', inv_style
          ierr = -1

      end select
      close( 99 )

      write(6,*) 'GMRES energy steps: ', gmres_nsteps

      inquire( file="echamp.inp", exist=echamp )
      if( echamp ) then
        open( unit=99, file="echamp.inp", form="formatted", status="old" )
        read( 99, * ) echamp
        close( 99 )
      endif
    endif
100 continue
    if( myid .eq. 0 .and. ierr .ne. 0 ) then
      write(6,*) 'Error reading gmres.in/bse.in!'
    endif
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    call MPI_BCAST( gmres_depth, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_resolution, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_convergence, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_preconditioner, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_nsteps, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( myid .ne. root ) allocate( gmres_energy_list( gmres_nsteps ) )
    call MPI_BCAST( gmres_energy_list, gmres_nsteps, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( echamp, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'
    call initialize_gmres_storage( sys, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'

  end subroutine SHG_gmres_setup

  !SJ BEGIN read-only GMRES parameter interface for SHG
  !> @brief Return the exact energy shifts and controls configured by setup.
  !>
  !> @details This is intentionally read-only. It lets an SHG driver assemble
  !! and validate [(omega+i*eta)-H] before entering the iterative solver,
  !! without duplicating the bse.in parser or changing the GMRES algorithm.
  subroutine SHG_gmres_get_parameters( energy_ha, resolution_ha, &
                                         preconditioner_ha, convergence, depth, ierr )
    use OCEAN_constants, only : eV2Hartree

    implicit none

    real(DP), allocatable, intent(out) :: energy_ha(:)
    real(DP), intent(out) :: resolution_ha
    real(DP), intent(out) :: preconditioner_ha
    real(DP), intent(out) :: convergence
    integer, intent(out) :: depth
    integer, intent(inout) :: ierr

    if( .not. allocated(gmres_energy_list) .or. gmres_nsteps .lt. 1 ) then
      ierr = 7111
      return
    endif

    allocate( energy_ha(gmres_nsteps), stat=ierr )
    if( ierr .ne. 0 ) return

    energy_ha(:) = gmres_energy_list(:) * eV2Hartree
    resolution_ha = gmres_resolution
    preconditioner_ha = gmres_preconditioner
    convergence = gmres_convergence
    depth = gmres_depth

  end subroutine SHG_gmres_get_parameters
  !SJ END read-only GMRES parameter interface for SHG

  subroutine initialize_gmres_storage( sys, ierr )
    use OCEAN_psi, only : OCEAN_psi_init, OCEAN_psi_new, OCEAN_psi_size_min
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: i

    allocate( au_matrix( gmres_depth ) )
    allocate( u_matrix( gmres_depth ) )

    call OCEAN_psi_init( sys, ierr )
    if( ierr .ne. 0 ) return

    do i = 1, gmres_depth
      call OCEAN_psi_new( au_matrix( i ), ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( u_matrix( i ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    if( myid .eq. root ) then
      i = OCEAN_psi_size_min(au_matrix(1))
      write(6,'(A,F8.2,A)') 'GMRES requires additional: ', gmres_depth*dble(i)*(32.0_dp/1048576.0_dp), ' MB'
    endif

    if( myid .eq. root ) then
      allocate( c_matrix( gmres_depth*(gmres_depth+1)/2 ) )
    else
      allocate( c_matrix( 1 ) )
    endif

  end subroutine initialize_gmres_storage

  subroutine make_recycle( current_iter, ener, prev_ener, psi_tmp, sys, ierr )
    use OCEAN_psi, only : OCEAN_psi_axpy, OCEAN_psi_min_set_prec, OCEAN_psi_3element_mult
    use OCEAN_energies, only : OCEAN_energies_sfact_copy, OCEAN_energies_allow
    use OCEAN_system, only : o_system

    integer, intent( in ) :: current_iter
    real(DP), intent( in ) :: ener, prev_ener
    type(OCEAN_vector), intent(inout) :: psi_tmp
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: i
    real(DP) :: delta

    delta = ener - prev_ener

    if( sys%bwflg ) then
      do i = 1, current_iter
        call OCEAN_energies_sfact_copy( sys, au_matrix(i), psi_tmp, ierr )
        if( ierr .ne. 0 ) return
        call OCEAN_psi_axpy( delta, u_matrix( i ), psi_tmp, ierr )
        if( ierr .ne. 0 ) return
        call OCEAN_energies_sfact_copy( sys, psi_tmp, au_matrix(i), ierr )
        if( ierr .ne. 0 ) return
      enddo
    else
      do i = 1, current_iter
        call OCEAN_psi_axpy( delta, u_matrix( i ), au_matrix( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
    endif

  end subroutine make_recycle

  subroutine update_gmres( current_iter, psi_g, psi_x, psi_ax, ierr)
    use OCEAN_mpi
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_dot, OCEAN_psi_axpy, OCEAN_psi_nrm
    use OCEAN_timekeeper, only : tk_inv, OCEAN_tk_start, OCEAN_tk_stop
    !
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent(inout ) :: psi_g, psi_x, psi_ax
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: c_temp( : ), coeff( : )
    real(DP), allocatable, dimension(:) :: re_coeff_vec, im_coeff_vec, re_rvec, im_rvec 
    integer, allocatable :: ipiv( : )
    real(DP) :: tmp_r, tmp_i
    integer :: local_gmres_size, info, iter, c_size, iter_start
    !
    call OCEAN_tk_start( tk_inv )
    !
    allocate( re_coeff_vec( current_iter ), im_coeff_vec( current_iter ), & 
              re_rvec( current_iter ), im_rvec( current_iter ) )    

    !
    !
    ! For iter = current_iter we can call psi_nrm instead of psi_dot
    !  this will save a small amount of time, but need to correctly 
    !  set the imaginary value and request to 0/null 
    re_rvec( : ) = 0.0_DP
    im_rvec( : ) = 0.0_DP
    re_coeff_vec( : ) = 0.0_DP
    im_coeff_vec( : ) = 0.0_DP
    !
    do iter = 1, current_iter - 1
      call OCEAN_psi_dot( au_matrix( current_iter ), au_matrix( iter ), &
                          re_coeff_vec( iter ), ierr, &
                          ival=im_coeff_vec( iter ), defer=.true. )

      call OCEAN_psi_dot( au_matrix( iter ), psi_g, &
                          re_rvec( iter ), ierr , &
                          ival=im_rvec( iter ), defer=.true. )
    enddo

    call OCEAN_psi_nrm( re_coeff_vec( current_iter ), au_matrix( current_iter ), ierr, defer = .true. ) 
    call OCEAN_psi_dot( au_matrix( current_iter ), psi_g, &
                        re_rvec( current_iter ), ierr , &
                        ival=im_rvec( current_iter ), defer=.true. )

    if( myid .eq. root ) then
      call MPI_REDUCE( MPI_IN_PLACE, re_coeff_vec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( MPI_IN_PLACE, im_coeff_vec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( MPI_IN_PLACE, re_rvec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( MPI_IN_PLACE, im_rvec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
    else
      call MPI_REDUCE( re_coeff_vec, re_coeff_vec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( im_coeff_vec, im_coeff_vec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( re_rvec, re_rvec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
      call MPI_REDUCE( im_rvec, im_rvec, current_iter, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr )
    endif
    !
    allocate( coeff( current_iter ) )
    if( myid .eq. root ) then
      iter_start = (current_iter * ( current_iter - 1 ) ) / 2
      do iter = 1, current_iter
        c_matrix( iter + iter_start ) = cmplx( re_coeff_vec( iter ), -im_coeff_vec( iter ), DP )
      enddo
      
      c_size = current_iter * (current_iter + 1 ) / 2 
      allocate( c_temp( c_size ), ipiv( current_iter ) )
      c_temp( : ) = c_matrix( 1 : c_size )

      call zhptrf( 'U', current_iter, c_temp, ipiv, info )

      ! Want to test after waitall so that all procs are on the same page if we abort
      if( info .ne. 0 ) then
        ierr = info
      else
        do iter = 1, current_iter
          coeff( iter ) = -cmplx( re_rvec( iter ), im_rvec( iter ), DP )
        enddo


        call zhptrs( 'U', current_iter, 1, c_temp, ipiv, coeff, current_iter, info )
      endif
      deallocate( c_temp, ipiv )

    endif

    call MPI_BCAST( info, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    ierr = info
    if( ierr .ne. 0 ) return

    call MPI_BCAST( coeff, current_iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return

    ! update ax and x
    do iter = 1, current_iter
      tmp_r = real( coeff( iter ), DP )
      tmp_i = aimag( coeff( iter ) )
      call OCEAN_psi_axpy( tmp_r, u_matrix( iter ), psi_x, ierr, tmp_i )
      call OCEAN_psi_axpy( tmp_r, au_matrix( iter ), psi_ax, ierr, tmp_i )
    enddo

    call OCEAN_tk_stop( tk_inv )
    !
  end subroutine update_gmres

  subroutine eigensystem_gmres( iter, keep, ierr )
!    use ocean_mpi, only : root, myid, MPI_REQUEST_NULL,  MPI_STATUSES_IGNORE
    use ocean_mpi
    use OCEAN_psi
    integer, intent( in ) :: iter
    integer, intent( inout ) :: keep
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: hess( :,: ), w(:), z(:,:), work(:), vl(:,:), vr(:,:)
    real(DP), allocatable :: re_hess(:,:), im_hess(:,:), re_coeff_vec(:), im_coeff_vec(:), rwork(:)
    integer, allocatable :: re_request(:,:), im_request(:,:), keep_list(:)
    real(DP) :: itmp, rtmp
    integer :: i, j, lwork, k, l, iter_start, ldvl
    type(ocean_vector), allocatable :: new_u(:)

    
    allocate( re_hess( iter, iter ), im_hess( iter, iter ),  &
              re_request( iter, iter ), im_request( iter, iter ), STAT=ierr )
    if( ierr .ne. 0 ) return

    re_hess = 0.0_DP
    im_hess = 0.0_DP
    re_request(:,:) = MPI_REQUEST_NULL
    im_request(:,:) = MPI_REQUEST_NULL
    
    do i = 1, iter
      do j = 1, min( i+1,iter)
!TODO: This hasn't been re-thought in terms of minimizing comms
        call OCEAN_psi_dot( u_matrix( j ), au_matrix( i ), &
                            re_hess( j,i ), ierr, &
                            ival=im_hess( j,i ), rrequest=re_request(j,i), &
                            irequest=im_request( j,i ) )
      enddo
    enddo

    do i = 1, iter
      call MPI_WAITALL( iter, re_request(:,i), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( iter, im_request(:,i), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
    enddo
    
    allocate( w(iter), z(iter,iter) )
    if( myid .eq. root ) then
      lwork = iter * 12
      allocate( hess( iter, iter ), work( lwork ) )
      hess(:,:) = cmplx( re_hess(:,:), im_hess(:,:), DP )

      call zhseqr('S', 'I', iter, 1, iter, hess, iter, w, z, iter, work, lwork, ierr )
      if( ierr .ne. 0 ) then
        keep = 0 
      endif
      ldvl = iter
      allocate( vl(iter,iter), vr( iter, iter ), rwork( iter ) )
      call ztrevc( 'R', 'A', .true., iter, z, iter, vl, ldvl, vr, iter, iter, l, work, rwork, ierr)
      ierr = 0
      z(:,:) = vr(:,:)
      deallocate( hess, work, vl, vr )
    endif
    call MPI_BCAST( w, iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( z, iter*iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( keep, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    keep = min( keep, iter )
    allocate( keep_list( keep ) )
    write(6,*) myid, keep, iter
    keep_list(:) = 1
    do i = 2, keep
      k = i
      do j = 1, i-1
        if( abs(w(k)) .lt. abs(w(keep_list(j))) ) then
          l = keep_list(j)
          keep_list( j ) = k
          k = l
        endif
      enddo
      keep_list( i ) = k
    enddo

    do i = keep+1, iter
      k = i
      do j = 1, keep
        if( abs(w(k)) .lt. abs(w(keep_list(j))) ) then
          keep_list( j ) = k
          k = j
        endif
      enddo
    enddo


    do j = 1, keep
      i = keep_list(j)
      if( myid .eq. root ) write(6,*) 'EIG-S:', i, real(abs(w(i)))
       write(6,*) 'EIG-S:', myid, i, real(abs(w(i)))
    enddo

    call MPI_BARRIER( comm, ierr )

    allocate( new_u(keep) )
    do j = 1, keep
      call OCEAN_psi_new( new_u( j ), ierr )
      call OCEAN_psi_zero_min( new_u(j), ierr )
      if( ierr .ne. 0 ) return
      k = keep_list(j)
      do i = 1, iter
        if( k.lt.1 .or. k.gt.iter ) write(6,*) myid, j, k, keep_list(:)
        rtmp = real( z(i,k), DP )
        itmp = aimag( z(i,k) )
        call OCEAN_psi_axpy( rtmp, u_matrix(i), new_u(j), ierr, itmp )
        if( ierr .ne. 0 ) return
      enddo
    enddo

    ! update u and au 
    do j = 1, keep
      call OCEAN_psi_copy_min( u_matrix(j), new_u(j), ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_zero_min( au_matrix(j), ierr )
      if( ierr .ne. 0 ) return
      rtmp = real( w(keep_list(j) ), DP )
      itmp = aimag( w(keep_list(j) ) )
      call OCEAN_psi_axpy( rtmp, new_u(j), au_matrix(j), ierr, itmp )
      if( ierr .ne. 0 ) return
      
    enddo

    deallocate(keep_list, new_u)
    deallocate( re_hess, im_hess, re_request, im_request, w, z )

    keep = 1
    return
    allocate( re_coeff_vec(keep), im_coeff_vec(keep), re_request(keep,1), im_request(keep,1) )
    re_request(:,:) = MPI_REQUEST_NULL
    im_request(:,:) = MPI_REQUEST_NULL
    do j = 1, keep
      re_coeff_vec(:) = 0.0_DP
      im_coeff_vec(:) = 0.0_DP
      do i = 1, j-1
        call OCEAN_psi_dot( au_matrix( j ), au_matrix( i ), &
                            re_coeff_vec( i ), ierr, &
                            ival=im_coeff_vec( i ), & 
                            rrequest=re_request( i,1 ), irequest=im_request( i,1 ) )
      enddo
      call OCEAN_psi_nrm( re_coeff_vec( j ), au_matrix( j ), ierr, &
                        re_request( j,1 ) )
      call MPI_WAITALL( j, re_request(1:j,1), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( j, im_request(1:j,1), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      iter_start = ( j * ( j - 1 ) ) / 2
      do i = 1, j
        c_matrix( i + iter_start ) = cmplx( re_coeff_vec( i ), -im_coeff_vec( i ), DP )
      enddo
    enddo

    deallocate( re_coeff_vec, im_coeff_vec, re_request, im_request )

  end subroutine eigensystem_gmres

  !SJ BEGIN return converged GMRES vectors to nonlinear-response callers
  subroutine SHG_gmres_do_recycle( sys, hay_vec, ierr, solution_vectors, &
                                     solution_residual_norm, solution_iterations, &
                                     energy_override_ha, resolution_override_ha, &
                                     hamiltonian_sign_override, apply_sfact_override, &
                                     solver_only_override )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_action, only : OCEAN_xact, OCEAN_action_h1
    use OCEAN_constants, only : eV2Hartree, Hartree2eV
    use OCEAN_mpi, only : myid, root, comm, MPI_STATUSES_IGNORE, MPI_REQUEST_NULL, MPI_INTEGER
    use OCEAN_energies, only : OCEAN_energies_allow_full, OCEAN_energies_sfact_copy, OCEAN_energies_allow
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( inout ) :: ierr
    ! Optional outputs leave the legacy spectrum-only call unchanged.  When
    ! requested, solution_vectors(j) contains the converged response to the
    ! normalized solver RHS at gmres_energy_list(j).  The SHG layer restores
    ! the physical optical-source norm after this routine returns.
    type( ocean_vector ), allocatable, intent( out ), optional :: solution_vectors(:)
    real(DP), allocatable, intent( out ), optional :: solution_residual_norm(:)
    integer, allocatable, intent( out ), optional :: solution_iterations(:)
    ! Arbitrary-source SHG solves use these controls to reuse the same GMRES
    ! implementation for G_-(Omega) and G_+(Omega), without changing the
    ! configured linear-spectrum defaults.
    real(DP), intent(in), optional :: energy_override_ha(:)
    real(DP), intent(in), optional :: resolution_override_ha
    integer, intent(in), optional :: hamiltonian_sign_override
    logical, intent(in), optional :: apply_sfact_override
    logical, intent(in), optional :: solver_only_override
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g, psi_pcdiv, psi_tmp
    type( ocean_vector), pointer :: psi_pg, psi_apg
    !
    real(DP) :: ener, rval, ival, gval, rel_error, fact, prev_energy
    real(DP) :: active_resolution
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: zero = 0.0_dp
    integer :: iter, step_iter, complete_iter, requests( 2 ), ierr_, prev_iter, global_iter, keep, i
    logical :: loud = .true.
    character( len=25 ) :: abs_filename
    integer, parameter :: abs_fh = 76
    type( ocean_vector ) :: hay_psi_x, psi_rhs
    integer :: solution_i, active_nsteps, active_hamiltonian_sign
    integer :: nhflag(6)
    logical :: active_apply_sfact, solver_only
    real(DP), allocatable :: active_energy_ha(:)
  !SJ END return converged GMRES vectors to nonlinear-response callers

    !SJ BEGIN configure signed, arbitrary-source SHG resolvent solves
    active_nsteps = gmres_nsteps
    allocate( active_energy_ha(active_nsteps), stat=ierr )
    if( ierr .ne. 0 ) return
    active_energy_ha(:) = gmres_energy_list(:) * eV2Hartree
    if( present(energy_override_ha) ) then
      deallocate( active_energy_ha )
      active_nsteps = size(energy_override_ha)
      if( active_nsteps .lt. 1 ) then
        ierr = 7131
        return
      endif
      allocate( active_energy_ha(active_nsteps), stat=ierr )
      if( ierr .ne. 0 ) return
      active_energy_ha(:) = energy_override_ha(:)
    endif

    active_resolution = gmres_resolution
    if( present(resolution_override_ha) ) active_resolution = resolution_override_ha
    if( active_resolution .le. tiny(1.0_DP) ) then
      ierr = 7132
      return
    endif

    active_hamiltonian_sign = -1
    if( present(hamiltonian_sign_override) ) then
      active_hamiltonian_sign = hamiltonian_sign_override
    endif
    if( abs(active_hamiltonian_sign) .ne. 1 ) then
      ierr = 7133
      return
    endif

    active_apply_sfact = .true.
    if( present(apply_sfact_override) ) active_apply_sfact = apply_sfact_override
    solver_only = .false.
    if( present(solver_only_override) ) solver_only = solver_only_override
    !SJ END configure signed, arbitrary-source SHG resolvent solves

    nhflag(:) = sys%nhflag

    requests( : ) = MPI_REQUEST_NULL

    if( sys%cur_run%backf ) then
      write(6,*) 'BACKF not implemented for GMRES (try Haydock)'
      ierr = 9963
      return
    endif

    if( sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    !SJ BEGIN suppress linear-spectrum files for internal SHG solves
    if( myid .eq. root .and. .not. solver_only ) then

      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref

      call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
      open( unit=abs_fh,file=abs_filename,form='formatted',status='unknown' )
      rewind( abs_fh )

    endif
    !SJ END suppress linear-spectrum files for internal SHG solves
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif


    ! Create all the psi vectors we need
    call OCEAN_psi_new( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( psi_g, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( psi_rhs, ierr )
    if( ierr .ne. 0 ) return
    !SJ BEGIN allow a prepared nonlinear-response vector to be the exact RHS
    if( active_apply_sfact ) then
      call OCEAN_energies_sfact_copy( sys, hay_vec, psi_rhs, ierr )
    else
      call OCEAN_psi_copy_min( psi_rhs, hay_vec, ierr )
    endif
    if( ierr .ne. 0 ) return
    !SJ END allow a prepared nonlinear-response vector to be the exact RHS
    !!!!

    !SJ BEGIN allocate optional converged-vector outputs once per energy grid
    if( present(solution_vectors) ) then
      allocate( solution_vectors(active_nsteps), stat=ierr )
      if( ierr .ne. 0 ) return
      do solution_i = 1, active_nsteps
        call OCEAN_psi_new( solution_vectors(solution_i), ierr )
        if( ierr .ne. 0 ) return
      enddo
    endif
    if( present(solution_residual_norm) ) then
      allocate( solution_residual_norm(active_nsteps), stat=ierr )
      if( ierr .ne. 0 ) return
      solution_residual_norm(:) = huge(1.0_DP)
    endif
    if( present(solution_iterations) ) then
      allocate( solution_iterations(active_nsteps), stat=ierr )
      if( ierr .ne. 0 ) return
      solution_iterations(:) = 0
    endif
    !SJ END allocate optional converged-vector outputs once per energy grid
    
    if( sys%cur_run%bwflg ) then
      call OCEAN_psi_new( psi_tmp, ierr )
      if( ierr .ne. 0 ) return
    endif


    ! Create the psi_min_arrays we will need
    ! AX etc

    ! Keep around H * (1,0) for making pcdiv for each loop
    call OCEAN_psi_one_full( psi_x, ierr )
    if( ierr .ne. 0 ) return

    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow_full( sys, psi_x, ierr, sfact=.true. )
      if( ierr .ne. 0 ) return
    endif

    call OCEAN_xact( sys, sys%interactionScale, psi_x, hpsi1, ierr,.false. , nhflag)
    if( ierr .ne. 0 ) return
    !SJ BEGIN match the diagonal preconditioner to (Omega+i*eta)+/-H
    call OCEAN_psi_scal( -real(active_hamiltonian_sign,DP), hpsi1, ierr )
    if( ierr .ne. 0 ) return
    !SJ END match the diagonal preconditioner to (Omega+i*eta)+/-H
    

    global_iter = 0
    prev_iter = 0
    prev_energy = -huge(prev_energy)

    !SJ BEGIN iterate over the active signed SHG frequency grid
    do step_iter = 1, active_nsteps
      complete_iter = 0

!      ener = e_list( iter ) * eV2Hartree 
      ener = active_energy_ha( step_iter )
      if( myid .eq. root ) write(6,'(F20.8)') ener * Hartree2eV


      if( .not. allow_recycle ) then
        prev_iter = 0
      endif
      ! are we new or recycling?
      !   need to tune when we throw out because the preconditioner is no longer valid
      if( do_precondition .and. prev_iter .eq. 0 ) then 
        prev_energy = ener
        call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
          return
        endif
        if( sys%bwflg ) then
          call OCEAN_energies_allow( sys, psi_pcdiv, ierr, sfact=.true. )
          if( ierr .ne. 0 ) return
        endif
      endif


      if( prev_iter .eq. 0 ) then
!        if( myid .eq. root ) write( 6, * ) 'no recycle: '
        !SJ BEGIN pass the active SHG shift/sign into the initial residual
        call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, psi_rhs, psi_tmp, ierr, nhflag, &
                                 ener, active_resolution, active_hamiltonian_sign )
        !SJ END pass the active SHG shift/sign into the initial residual
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write( 6,*) 'set_initial_vector failed'
          return
        endif
        if( loud ) then
          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return

          call OCEAN_psi_dot( psi_rhs, psi_x, rval, ierr, ival=ival )
          if( ierr .ne. 0 ) return

          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, prev_iter, gmres_depth, &
                  gval, gmres_convergence, ener, (1.0_dp - rval), -ival
          endif
        endif
      else
        ! recycle au_matrix
        if( myid .eq. root ) write( 6, * ) 'recycle: ', prev_iter, &
              (active_energy_ha( step_iter ) - active_energy_ha( step_iter -1)), &
              (ener - prev_energy)
        call make_recycle( prev_iter, active_energy_ha( step_iter ), &
                           active_energy_ha( step_iter -1), psi_tmp, sys, ierr )
        if( ierr .ne. 0 ) return

        !SJ BEGIN pass the active SHG shift/sign into the recycled residual
        call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, psi_rhs, psi_tmp, ierr, nhflag, &
                                 ener, active_resolution, active_hamiltonian_sign )
        !SJ END pass the active SHG shift/sign into the recycled residual
        if( ierr .ne. 0 ) return
        do iter = 1, prev_iter
          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_axmz( psi_ax, psi_g, psi_rhs, ierr )
        enddo
          !

        call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
        if( ierr .ne. 0 ) return
        if( loud ) then
          call OCEAN_psi_dot( psi_rhs, psi_x, rval, ierr, ival=ival )
          if( ierr .ne. 0 ) return
          call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                  gval, gmres_convergence, ener, (1.0_dp - rval), -ival
          endif
        else
          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,3(1x,1e15.8),i8)' ) complete_iter, iter, gmres_depth, &
                  gval, gmres_convergence, ener, global_iter
          endif
        endif
!          if( myid.eq.root) write(*,*) 'conv check'
        if( gval .lt. gmres_convergence ) then
          !SJ BEGIN use the true residual criterion for arbitrary SHG sources
          if( solver_only ) goto 200
          !SJ END use the true residual criterion for arbitrary SHG sources
          call OCEAN_psi_dot( psi_rhs, psi_x, rval, ierr, ival=ival )
          if( ierr .ne. 0 ) return
          call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          rel_error = -gval / ival

          if( rel_error .lt. 0.001 .or. (gval .lt. gmres_convergence/10.0_dp ) ) then
            goto 200 ! if convergered goto 200
          endif
        endif
        if( abs( prev_energy - ener ) > gmres_preconditioner/4.0_DP ) then
          if( myid .eq. root ) write(6,*) 'PRECONDITIONER OUT OF RANGE'
          prev_iter = 0
          prev_energy = ener
          call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
            return
          endif
          if( sys%bwflg ) then
            call OCEAN_energies_allow( sys, psi_pcdiv, ierr, sfact=.true. )
            if( ierr .ne. 0 ) return
          endif
        endif
        !
      endif

      do ! outerloop for restarted GMRES
      
        do iter = prev_iter+1, gmres_depth
          complete_iter = complete_iter + 1
          global_iter = global_iter + 1
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          if( do_precondition) then
            call OCEAN_psi_3element_mult( psi_pg, psi_g, psi_pcdiv, ierr )
          else
            call OCEAN_psi_copy_min( psi_pg, psi_g, ierr )
          endif
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'do_precondition failed'
            return
          endif

          call OCEAN_psi_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return

          ! apg = H . pg
          !TODO
          ! if bwflg then apg = H . S . pg
          if( sys%bwflg ) then
            !TODO
            call OCEAN_energies_sfact_copy( sys, psi_pg, psi_tmp, ierr )
            if( ierr .ne. 0 ) return
            call OCEAN_xact( sys, sys%interactionScale, psi_tmp, psi_apg, ierr )
          else
            call OCEAN_xact( sys, sys%interactionScale, psi_pg, psi_apg, ierr )
          endif
          if( ierr .ne. 0 ) return

          ! apg = (e+iG) * pg - apg
          !SJ BEGIN apply the selected G-minus or G-plus shifted operator
          call SHG_gmres_signed_shift( psi_pg, psi_apg, ener, active_resolution, &
                                         active_hamiltonian_sign, ierr )
          !SJ END apply the selected G-minus or G-plus shifted operator
          if( ierr .ne. 0 ) return

          !TODO:
          !if bwflg then apg = S. apg          
          if( sys%bwflg ) then
            call OCEAN_energies_allow( sys, psi_apg, ierr, sfact=.true. )
            if( ierr .ne. 0 ) return
          endif

          call OCEAN_psi_free_fbe( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_apg, ierr )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )

          if( iter .eq. gmres_depth ) then
            call OCEAN_psi_min2full( psi_x, ierr )
            if( ierr .ne. 0 ) return
            ! newi2loop section here

            ! if bwflg then ax = H . S . x
            if( sys%cur_run%bwflg ) then
              call OCEAN_energies_sfact_copy( sys, psi_x, psi_tmp, ierr )
              if( ierr .ne. 0 ) return
              call OCEAN_xact( sys, sys%interactionScale, psi_tmp, psi_ax, ierr )
              if( ierr .ne. 0 ) return
            else
              call OCEAN_xact( sys, sys%interactionScale, psi_x, psi_ax, ierr )
              if( ierr .ne. 0 ) return
            endif

            ! psi_ax = (ener + iG ) * psi_x - psi_ax
            !SJ BEGIN apply the selected G-minus or G-plus restart operator
            call SHG_gmres_signed_shift( psi_x, psi_ax, ener, active_resolution, &
                                           active_hamiltonian_sign, ierr )
            !SJ END apply the selected G-minus or G-plus restart operator
            if( ierr .ne. 0 ) return
            !
            ! if bwflg then ax = S . ax          
            if( sys%bwflg ) then
              call OCEAN_energies_allow( sys, psi_ax, ierr, sfact=.true. )
              if( ierr .ne. 0 ) return
            endif
            !
            if( ener .ne. prev_energy ) then
              call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
              if( ierr .ne. 0 ) then
                if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
                return
              endif
              if( sys%bwflg ) then
                call OCEAN_energies_allow( sys, psi_pcdiv, ierr, sfact=.true. )
                if( ierr .ne. 0 ) return
              endif
              prev_energy = ener
            endif

            prev_iter = 0
          endif

          ! get new g
          ! g = ax - b   !!! y = x - z
          call OCEAN_psi_axmz( psi_ax, psi_g, psi_rhs, ierr )
          !

          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return
          if( loud ) then
            call OCEAN_psi_dot( psi_rhs, psi_x, rval, ierr, ival=ival )
            if( ierr .ne. 0 ) return

            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, (1.0_dp - rval), -ival
            endif
          else
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,3(1x,1e15.8),i8)' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, global_iter
            endif
          endif

          if( gval .lt. gmres_convergence ) then 

            !SJ BEGIN use the true residual criterion for arbitrary SHG sources
            if( solver_only ) then
              prev_iter = iter
              if( iter .eq. gmres_depth ) prev_iter = 0
              goto 200
            endif
            !SJ END use the true residual criterion for arbitrary SHG sources

            call OCEAN_psi_dot( psi_rhs, psi_x, rval, ierr, ival=ival )
            if( ierr .ne. 0 ) return
            rel_error = -gval / ival

            if( rel_error .lt. 0.001 .or. (gval .lt. gmres_convergence/10.0_dp ) ) then
              prev_iter = iter
              if( iter .eq. gmres_depth ) prev_iter = 0 
              goto 200 ! if convergered goto 200
            endif
          endif
        enddo


      enddo

!     Exit to here if we have converged
200   continue
!      prev_iter = iter

      !SJ BEGIN preserve each converged response before psi_x is reused
      if( present(solution_vectors) ) then
        call OCEAN_psi_copy_min( solution_vectors(step_iter), psi_x, ierr )
        if( ierr .ne. 0 ) return
      endif
      if( present(solution_residual_norm) ) then
        ! OCEAN_psi_nrm returns the squared Euclidean norm.
        solution_residual_norm(step_iter) = sqrt( max(0.0_DP, gval) )
      endif
      if( present(solution_iterations) ) solution_iterations(step_iter) = complete_iter
      !SJ END preserve each converged response before psi_x is reused

      !SJ BEGIN keep internal SHG resolvent solves free of linear-spectrum output
      if(myid.eq. root .and. .not. solver_only) then
        rel_error = -gval / ival
        write( abs_fh, '(1p,4(1e15.8,1x),1i6,1x,i10)' ) ener * Hartree2eV, &
                    (1.0_dp - rval*fact), -ival*fact, rel_error, complete_iter, global_iter
        flush(abs_fh)
      endif
      !SJ END keep internal SHG resolvent solves free of linear-spectrum output

      ! TODO: make this a separate flag
      !SJ BEGIN suppress field-amplitude products for internal SHG solves
      if ( echamp .and. .not. solver_only ) then

        rval = 0.0_DP
        call OCEAN_psi_3element_mult( hay_psi_x, psi_rhs, psi_x, ierr, rval, .true. )
        if( ierr .ne. 0 ) return
      endif

      if( echamp .and. .not. solver_only ) then

        call write_out_echamp(sys, step_iter, psi_x, ierr)
        if( ierr .ne. 0 ) return
      endif
      !SJ END suppress field-amplitude products for internal SHG solves
    enddo
    !SJ END iterate over the active signed SHG frequency grid

    call OCEAN_psi_kill( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_x, ierr )
       if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_g, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return


    if( sys%cur_run%bwflg ) then
      call OCEAN_psi_kill( psi_tmp, ierr )
      if( ierr .ne. 0 ) return
    endif

    !SJ BEGIN close only the spectrum file opened by this invocation
    if( myid .eq. root .and. .not. solver_only ) close( abs_fh )
    if( allocated(active_energy_ha) ) deallocate(active_energy_ha)
    !SJ END close only the spectrum file opened by this invocation


  end subroutine SHG_gmres_do_recycle

  subroutine SHG_gmres_do( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_action, only : OCEAN_xact, OCEAN_action_h1
    use OCEAN_constants, only : eV2Hartree, Hartree2eV
    use OCEAN_mpi, only : myid, root, comm, MPI_STATUSES_IGNORE, MPI_REQUEST_NULL, MPI_INTEGER
    use OCEAN_energies, only : OCEAN_energies_allow
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g, psi_pcdiv
    !AK-NOTE:
    !! the strategy to output intermediate eigenvalues is:
    !! inefficient route: save as a new ocean_vector after the right combination
    !of hamiltonian elements have been applied and allow the minimization to
    !proceed as programmed.
    !! efficient route: use previously minimized eigenvalues as initial guesses
    !for the subsequent gmres of similar eigenvalues -- should provide a shorter
    !gmres time. 
    !! The elements of the BSE hamiltonian applied in SHG_gmres_do are: 
    type( ocean_vector), pointer :: psi_pg, psi_apg
    !
    real(DP) :: ener, rval, ival, gval, rel_error, fact
    real(dp), parameter :: one = 1.0_dp
    integer :: iter, step_iter, complete_iter, requests( 2 ), ierr_
    logical :: loud = .false.
    character( len=25 ) :: abs_filename
    integer, parameter :: abs_fh = 76
    type( ocean_vector ) :: hay_psi_x, psi_tmp
    
    integer :: nhflag(6)
    nhflag(:) = sys%nhflag

    requests( : ) = MPI_REQUEST_NULL

    if( sys%cur_run%bwflg .or. sys%cur_run%backf ) then
      write(6,*) 'BWFLG not implemented for OCEAN_do_gmres (see recycle)'
      ierr = 9962
      return
    endif
    
    if( sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    if( myid .eq. root ) then

      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref

      call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
      open( unit=abs_fh,file=abs_filename,form='formatted',status='unknown' )
      rewind( abs_fh )

    endif
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    
    ! Create all the psi vectors we need
    call OCEAN_psi_new( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_g, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    !!!!
    call OCEAN_psi_new( psi_tmp, ierr )
    if( ierr .ne. 0 ) return



    ! Create the psi_min_arrays we will need
    ! AX etc

    ! Keep around H * (1,0) for making pcdiv for each loop
    call OCEAN_psi_one_full( psi_x, ierr )
    if( ierr .ne. 0 ) return
    
    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, psi_x, ierr )
      if( ierr .ne. 0 ) return
    endif

    call OCEAN_xact( sys, sys%interactionScale, psi_x, hpsi1, ierr,.false. , nhflag)
    if( ierr .ne. 0 ) return




    do step_iter = 1, gmres_nsteps
      complete_iter = 0
      iter = 0

!      ener = e_list( iter ) * eV2Hartree
      ener = gmres_energy_list( step_iter ) * eV2Hartree
      if( myid .eq. root ) write(6,*) ener * Hartree2eV

      if( do_precondition ) then
        call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
          return
        endif
      endif
!      write(6,*) 'prec:', psi_pcdiv%min_r(1,1), psi_pcdiv%min_i(1,1)

!      zener = cmplx( ener, gres, DP )

      ! do we start with x = 0 or with x = previous or something else?
      call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, hay_vec, psi_tmp, ierr, nhflag)
      if( ierr .ne. 0 ) then
        if( myid .eq. root ) write( 6,*) 'set_initial_vector failed'
        return
      endif

!      write(6,*) hay_vec%r(1,1,1), hay_vec%i(1,1,1)
!      write(6,*) psi_g%min_r(1,1), psi_g%min_i(1,1)

      if( loud ) then
        call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
        if( ierr .ne. 0 ) return

        call OCEAN_psi_dot( hay_vec, psi_x, rval, ierr, ival=ival )
        if( ierr .ne. 0 ) return
!        call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
!        if( ierr .ne. 0 ) return
        if( myid .eq. 0 ) then
          write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                gval, gmres_convergence, ener, (1.0_dp - rval), -ival
        endif
      endif


      ! could have some out maximum here, like size_of_psi / gmres_depth 
      do ! outerloop for restarted GMRES
    ! if(myid .eq. root)  write(6, *) 'write out unmodified hay_vec with itter+200'
    !call write_out_echamp(sys, step_iter+200, hay_vec, ierr)
    !if( ierr .ne. 0 ) return
    
        do iter = 1, gmres_depth ! inner loop for restarted GMRES
          complete_iter = complete_iter + 1
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          if( do_precondition) then
            call OCEAN_psi_3element_mult( psi_pg, psi_g, psi_pcdiv, ierr )
          else
            call OCEAN_psi_copy_min( psi_pg, psi_g, ierr )
          endif
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'do_precondition failed'
            return
          endif

          call OCEAN_psi_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return
!          call OCEAN_psi_start_min2full( psi_pg, ierr )
!          if( ierr .ne. 0 ) return
!          call OCEAN_psi_finish_min2full( psi_pg, ierr )
!          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, sys%interactionScale, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return
        
          ! apg = (e+iG) * pg - apg
          call OCEAN_psi_axmy( psi_pg, psi_apg, ierr, ener, gmres_resolution )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_apg, ierr )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )


          if( iter .eq. gmres_depth ) then
            call OCEAN_psi_min2full( psi_x, ierr )
            if( ierr .ne. 0 ) return
            ! newi2loop section here
            call OCEAN_xact( sys, sys%interactionScale, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            
            ! psi_ax = (ener + iG ) * psi_x - psi_ax
            call OCEAN_psi_axmy( psi_x, psi_ax, ierr, ener, gmres_resolution )
            if( ierr .ne. 0 ) return
            !
          endif

          ! get new g
          ! g = ax - b   !!! y = x - z
          call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
          !

          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return
          if( loud ) then
            call OCEAN_psi_dot( hay_vec, psi_x, rval, ierr, ival=ival )
            if( ierr .ne. 0 ) return
!            call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
!            if( ierr .ne. 0 ) return
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, (1.0_dp - rval), -ival
            endif
          else
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,3(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener
            endif
          endif
!          if( myid.eq.root) write(*,*) 'conv check'
          if( gval .lt. gmres_convergence ) goto 200 ! if convergered goto 200
!          if(myid.eq.root) write(*,*) 'no conv: ', step_iter
        enddo


      enddo

!     Exit to here if we have converged
200   continue

      call OCEAN_psi_dot( hay_vec, psi_x, rval, ierr, ival=ival )
      if( ierr .ne. 0 ) return
!      call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
!      if( ierr .ne. 0 ) return
      if(myid.eq. root) then
        rel_error = -gval / ival
        write( abs_fh, '(1p,4(1e15.8,1x),1i6)' ) ener * Hartree2eV, & 
                    (1.0_dp - rval*fact), -ival*fact, rel_error, complete_iter
        flush(abs_fh)
      endif

      ! TODO: make this a separate flag
      if ( echamp ) then
        ! Move create above
!        call OCEAN_psi_new( hay_psi_x, ierr, psi_x)
!        if( ierr .ne. 0 ) return
    
        rval = 0.0_DP
        call OCEAN_psi_3element_mult( hay_psi_x, hay_vec, psi_x, ierr, rval, .true. )
        if( ierr .ne. 0 ) return
!        call OCEAN_psi_dot_write(hay_vec,psi_x,hay_psi_x,requests(1),rval,ierr,requests(2),ival)
!        call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
!        if( ierr .ne. 0 ) return

!JTV
!        call write_out_energy(sys, step_iter, hay_psi_x, nhflag, ierr)
!        if( ierr .ne. 0 ) return
      endif

 
      if( echamp ) then
        
        call write_out_echamp(sys, step_iter, psi_x, ierr)
        if( ierr .ne. 0 ) return
      endif
    enddo    
!enddo 
   ! Create all the psi vectors we need
  
    call OCEAN_psi_kill( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_x, ierr )
       if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_g, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_tmp, ierr )
    if( ierr .ne. 0 ) return
    

    if( myid .eq. root ) close( abs_fh )
  
  end subroutine SHG_gmres_do

!SJ BEGIN signed shifted-operator and initial-vector support for SHG
  !> @brief Apply [(Omega+i*eta)+sign*H] to a vector after Hx is available.
  subroutine SHG_gmres_signed_shift( psi_x, psi_ax, energy_ha, resolution_ha, &
                                       hamiltonian_sign, ierr )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_scal, OCEAN_psi_axmy

    type(ocean_vector), intent(inout) :: psi_x, psi_ax
    real(DP), intent(in) :: energy_ha, resolution_ha
    integer, intent(in) :: hamiltonian_sign
    integer, intent(inout) :: ierr

    ! OCEAN_psi_axmy performs ax=(energy+i*resolution)*x-ax.  Scaling
    ! Hx by -sign first therefore gives the requested plus/minus branch.
    call OCEAN_psi_scal( -real(hamiltonian_sign,DP), psi_ax, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_axmy( psi_x, psi_ax, ierr, energy_ha, resolution_ha )
  end subroutine SHG_gmres_signed_shift

!> @brief Initializes X (and hence AX)
  subroutine set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec, psi_tmp, ierr, nhflag, &
                                 energy_override_ha, resolution_override_ha, hamiltonian_sign_override )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_zero_min, OCEAN_psi_axmz, OCEAN_psi_axmy, &
                          OCEAN_psi_min2full
    use OCEAN_action, only : OCEAN_xact
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : eV2Hartree
    use OCEAN_mpi, only : myid, root
    use OCEAN_energies, only : OCEAN_energies_sfact_copy, OCEAN_energies_allow
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( in ) :: iter
    type( ocean_vector ), intent( inout ) :: psi_x, psi_g, psi_ax, psi_tmp
    integer, intent( inout ) :: ierr
    !
    real(DP) :: ener, resolution
    character(len=4) :: initial = 'zero'
    integer, intent ( in ) :: nhflag(6)
    real(DP), intent(in), optional :: energy_override_ha, resolution_override_ha
    integer, intent(in), optional :: hamiltonian_sign_override
    integer :: hamiltonian_sign

    ener = gmres_energy_list(iter) * eV2Hartree
    if( present(energy_override_ha) ) ener = energy_override_ha
    resolution = gmres_resolution
    if( present(resolution_override_ha) ) resolution = resolution_override_ha
    hamiltonian_sign = -1
    if( present(hamiltonian_sign_override) ) hamiltonian_sign = hamiltonian_sign_override

    initial = 'zero'
    if( allow_reuse_x .and. .not. present(energy_override_ha) ) then
      ! option to reuse X
      if( iter .gt. 1 ) then
        if( ( gmres_energy_list( iter ) - gmres_energy_list( iter -1 ) ) * eV2Hartree & 
                .lt. 3.0_DP * gmres_resolution ) then
          initial = 'keep'
        endif
        if( force_reuse_x ) initial = 'keep'
      endif
    endif


    select case ( initial )

      case ( 'keep' )
        ! Keep X, generate psi_ax
        ! if( myid .eq. root ) write(6,*) '   Re-use previous X'
        call OCEAN_psi_min2full( psi_x, ierr )
        if( ierr .ne. 0 ) return
    
        ! if bwflg then ax = H . S . x
        if( sys%cur_run%bwflg ) then
          call OCEAN_energies_sfact_copy( sys, psi_x, psi_tmp, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_xact( sys, sys%interactionScale, psi_tmp, psi_ax, ierr )
          if( ierr .ne. 0 ) return
        else
          call OCEAN_xact( sys, sys%interactionScale, psi_x, psi_ax, ierr )
          if( ierr .ne. 0 ) return
        endif

        call SHG_gmres_signed_shift( psi_x, psi_ax, ener, resolution, hamiltonian_sign, ierr )
        if( ierr .ne. 0 ) return

        ! if bwflg then ax = S . ax          
        if( sys%cur_run%bwflg ) then
          call OCEAN_energies_allow( sys, psi_ax, ierr, sfact=.true. )
          if( ierr .ne. 0 ) return
        endif

      case default 
        ! or 'zero'
        call OCEAN_psi_zero_min( psi_x, ierr )
        call OCEAN_psi_zero_min( psi_ax, ierr )

    end select
    !
    call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
    !
  end subroutine set_initial_vector
!SJ END signed shifted-operator and initial-vector support for SHG

  subroutine write_out_energy( sys, current_iter, psi_x,hflag, ierr )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_min2full, OCEAN_psi_vtor, OCEAN_psi_size_full
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_filenames, only : OCEAN_filenames_energy
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent( inout ) :: psi_x
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: out_vec( : )
    integer :: hflag(6)
    character(len=35) :: filename
    !
    call OCEAN_filenames_energy( sys, filename, current_iter,hflag, ierr )
    if( ierr .ne. 0 ) return
    !
   
    call OCEAN_psi_min2full( psi_x, ierr )
   !  write(6, *) 'FLAGAK 12 ', ierr
        if( ierr .ne. 0 ) return
    !
    if( myid .eq. root ) then
      allocate( out_vec( OCEAN_psi_size_full( psi_x ) ) )
      call OCEAN_psi_vtor( psi_x, out_vec ) 
      !
      open( unit=99, file=filename, form='formatted', status='unknown' )
      rewind( 99 )
      write(99, *) out_vec
      close( 99 )
      if( ierr .ne. 0 ) return
      deallocate( out_vec )
    endif
  end subroutine write_out_energy

  subroutine write_out_echamp( sys, current_iter, psi_x, ierr )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_min2full, OCEAN_psi_vtor, OCEAN_psi_size_full
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_filenames, only : OCEAN_filenames_ehamp
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent( inout ) :: psi_x
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: out_vec( : )
    character(len=28) :: filename
    !
    call OCEAN_filenames_ehamp( sys, filename, current_iter, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_min2full( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    if( myid .eq. root ) then
      allocate( out_vec( OCEAN_psi_size_full( psi_x ) ) )
      call OCEAN_psi_vtor( psi_x, out_vec ) 
      !
      open( unit=99, file=filename, form='unformatted', status='unknown' )
      rewind( 99 )
      write(99) out_vec
      close( 99 )
      deallocate( out_vec )
    endif
  end subroutine write_out_echamp

end module SHG_gmres

