! MPI core <- valence hole vertex.
!
! Parallelization is over contiguous conduction-band blocks.  Each rank reads
! only its local n slab for one k point at a time and writes the corresponding
! slab of the compact core-sector output.  The cv.in format and output binary
! layout are identical to the serial shg_rbxc_cv program.
program shg_rbxc_cv_mpi
  use iso_fortran_env, only: dp=>real64, int32, int64
  use ieee_arithmetic, only: ieee_is_finite
  use mpi
  implicit none

  integer :: nc,nv,nk,nb,l,ns,nvs,vskip,ne,nm,na
  integer :: u,ios,p,v,k,m,h,e,b,a,j,nproj,ntot,nspin
  integer :: world_rank,world_size,active_rank,active_size,active_comm,color
  integer :: mpierr,abort_ierr,base,extra,nlocal,nfirst,nlast
  integer :: xfh,yfh,xtype,ytype,status_mpi(MPI_STATUS_SIZE)
  integer :: xsizes(4),xsubs(4),xstarts(4)
  integer :: ysizes(3),ysubs(3),ystarts(3)
  integer(int32) :: header(3)
  integer(int64) :: expected,complex_count
  integer(MPI_OFFSET_KIND) :: zero_offset,output_bytes
  character(4096) :: control,cksfile,melfile,line
  character(4096), allocatable :: infile(:),outfile(:)
  character(1) :: bax,cax
  character(32) :: element
  complex(dp), allocatable :: proj(:,:,:,:),mel(:,:),d(:,:,:,:)
  complex(dp), allocatable :: x(:,:,:),y(:,:)
  real(dp), allocatable :: omega(:),omegaev(:)
  real(dp) :: re,im,eta,q(3),local_norm2,global_norm2,ynorm
  integer :: bands(4),mesh(3),znl(3),site
  logical :: bad_local,bad_global

  call MPI_Init(mpierr)
  if(mpierr/=MPI_SUCCESS) error stop 'MPI_Init failed'
  call MPI_Comm_rank(MPI_COMM_WORLD,world_rank,mpierr)
  call mpi_check(mpierr,'MPI_Comm_rank')
  call MPI_Comm_size(MPI_COMM_WORLD,world_size,mpierr)
  call mpi_check(mpierr,'MPI_Comm_size')

  if(command_argument_count()/=1) then
    if(world_rank==0) write(*,'(A)') &
      'Usage: shg_rbxc_cv.x cv.in (run in prepared directory)'
    call abort_world('invalid command line',7001)
  endif
  call get_command_argument(1,control)

  ! Root parses the complete control file once.  The fixed-length path arrays
  ! are then broadcast so every active rank can participate in MPI-IO.
  if(world_rank==0) then
    open(newunit=u,file=trim(control),status='old',iostat=ios)
    call io_check(ios,'open control')
    read(u,*,iostat=ios) nc,nv,nk,nb,l,ns,nvs,vskip,ne
    call io_check(ios,'control dimensions')
    read(u,*,iostat=ios) bax,cax
    call io_check(ios,'axes')
    read(u,*,iostat=ios) bands,mesh,znl,site,element
    call io_check(ios,'metadata')
    read(u,*,iostat=ios) eta,q
    call io_check(ios,'energy/q metadata')
    read(u,'(A)',iostat=ios) cksfile
    call io_check(ios,'projector path')
    read(u,'(A)',iostat=ios) melfile
    call io_check(ios,'mels path')
    allocate(omega(ne),omegaev(ne),infile(ne),outfile(ne),stat=ios)
    if(ios/=0) call abort_world('allocate control arrays',7002)
    do j=1,ne
      read(u,*,iostat=ios) omega(j),omegaev(j)
      call io_check(ios,'energy entry')
      read(u,'(A)',iostat=ios) infile(j)
      call io_check(ios,'input file entry')
      read(u,'(A)',iostat=ios) outfile(j)
      call io_check(ios,'output file entry')
    enddo
    close(u,iostat=ios)
    call io_check(ios,'close control')
  endif

  call MPI_Bcast(nc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nc')
  call MPI_Bcast(nv,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nv')
  call MPI_Bcast(nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nk')
  call MPI_Bcast(nb,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nbeta')
  call MPI_Bcast(l,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast l')
  call MPI_Bcast(ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nspin')
  call MPI_Bcast(nvs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nval stored')
  call MPI_Bcast(vskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast valence offset')
  call MPI_Bcast(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast energy count')

  if(world_rank/=0) then
    allocate(omega(ne),omegaev(ne),infile(ne),outfile(ne),stat=ios)
    if(ios/=0) call abort_world('allocate broadcast control arrays',7003)
  endif
  call MPI_Bcast(bax,1,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast middle axis')
  call MPI_Bcast(cax,1,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast right axis')
  call MPI_Bcast(bands,4,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast bands')
  call MPI_Bcast(mesh,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast mesh')
  call MPI_Bcast(znl,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast ZNL')
  call MPI_Bcast(site,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast site')
  call MPI_Bcast(element,len(element),MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast element')
  call MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast eta')
  call MPI_Bcast(q,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast q')
  call MPI_Bcast(cksfile,len(cksfile),MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast projector path')
  call MPI_Bcast(melfile,len(melfile),MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast mels path')
  call MPI_Bcast(omega,ne,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast omega')
  call MPI_Bcast(omegaev,ne,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast omega eV')
  call MPI_Bcast(infile,len(infile(1))*ne,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast input paths')
  call MPI_Bcast(outfile,len(outfile(1))*ne,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast output paths')

  if(min(nc,nv,nk,nvs,ne)<1.or.l<0.or.l>8) &
    call abort_world('invalid dimensions',7010)
  if(ns/=1.and.ns/=2) call abort_world('only collinear nspin=1,2 supported',7011)
  if(nb/=1.and.nb/=4) call abort_world('only nbeta=1 or 4 supported',7012)
  if(nb==1.and.ns/=1) call abort_world('scalar valence requires nspin=1',7013)
  if(vskip<0.or.vskip+nv>nvs) &
    call abort_world('valence slice outside projector file',7014)
  if(product(mesh)/=nk.or.znl(3)/=l) &
    call abort_world('metadata dimensions disagree',7015)
  if(.not.ieee_is_finite(eta).or..not.all(ieee_is_finite(q))) &
    call abort_world('nonfinite eta/q metadata',7016)

  nm=2*l+1
  na=4*nm

  ! Root reads the site projectors and radial/angular matrix elements, builds
  ! the much smaller d(v,k,m,spin), and broadcasts d.  The large conduction
  ! dimension never appears in this replicated object.
  if(world_rank==0) then
    open(newunit=p,file=trim(cksfile),status='old',access='stream', &
         form='unformatted',iostat=ios)
    call io_check(ios,'open parcksv')
    read(p,iostat=ios) header
    call io_check(ios,'parcksv header')
    nproj=int(header(1)); ntot=int(header(2)); nspin=int(header(3))
    if(nproj<1.or.ntot/=nvs*nk.or.nspin/=ns) &
      call abort_world('projector dimensions disagree',7020)
    expected=12_int64+16_int64*nproj*nvs*nk*ns
    call require_bytes_root(cksfile,expected)
  endif
  call MPI_Bcast(nproj,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast nproj')

  complex_count=int(nv,int64)*nk*nm*ns
  if(complex_count>int(huge(0),int64)) &
    call abort_world('dipole broadcast exceeds MPI count range',7021)
  allocate(d(nv,nk,nm,ns),stat=ios)
  if(ios/=0) call abort_world('allocate local dipoles',7022)

  if(world_rank==0) then
    allocate(proj(nproj,nvs,nk,ns),mel(nproj,nm),stat=ios)
    if(ios/=0) call abort_world('allocate projector work arrays',7023)
    read(p,iostat=ios) proj
    call io_check(ios,'projector payload')
    close(p,iostat=ios)
    call io_check(ios,'close parcksv')
    if(.not.all(ieee_is_finite(real(proj))).or. &
       .not.all(ieee_is_finite(aimag(proj)))) &
      call abort_world('nonfinite projectors',7024)

    open(newunit=p,file=trim(melfile),status='old',iostat=ios)
    call io_check(ios,'open mels')
    do m=1,nm
      do j=1,nproj
        read(p,*,iostat=ios) re,im
        call io_check(ios,'mels entry')
        if(.not.ieee_is_finite(re).or..not.ieee_is_finite(im)) &
          call abort_world('nonfinite mels',7025)
        mel(j,m)=cmplx(re,im,dp)
      enddo
    enddo
    do
      read(p,'(A)',iostat=ios) line
      if(ios<0) exit
      call io_check(ios,'mels trailing data')
      if(len_trim(line)>0) call abort_world('extra entries in mels',7026)
    enddo
    close(p,iostat=ios)
    call io_check(ios,'close mels')

    ! OCEAN cksv_read convention:
    !   xes = sum_p mel(p,m)*proj(p,v,k,spin)
    !   d   = <v|r_b|g> = conjg(xes)
    do h=1,ns
      do m=1,nm
        do k=1,nk
          do v=1,nv
            d(v,k,m,h)=conjg(sum(mel(:,m)*proj(:,vskip+v,k,h)))
          enddo
        enddo
      enddo
    enddo
    deallocate(proj,mel)
  endif
  call MPI_Bcast(d,int(complex_count),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
  call mpi_check(mpierr,'broadcast local dipoles')

  ! At most nc ranks are useful for a conduction-band decomposition.  Extra
  ! ranks remain in MPI_COMM_WORLD and wait in MPI_Finalize.
  active_size=min(world_size,nc)
  if(world_rank<active_size) then
    color=1
  else
    color=MPI_UNDEFINED
  endif
  call MPI_Comm_split(MPI_COMM_WORLD,color,world_rank,active_comm,mpierr)
  call mpi_check(mpierr,'MPI_Comm_split')

  if(color/=MPI_UNDEFINED) then
    call MPI_Comm_rank(active_comm,active_rank,mpierr)
    call mpi_check(mpierr,'active rank')
    call MPI_Comm_size(active_comm,active_size,mpierr)
    call mpi_check(mpierr,'active size')

    base=nc/active_size
    extra=mod(nc,active_size)
    if(active_rank<extra) then
      nlocal=base+1
      nfirst=active_rank*(base+1)+1
    else
      nlocal=base
      nfirst=extra*(base+1)+(active_rank-extra)*base+1
    endif
    nlast=nfirst+nlocal-1

    allocate(x(nlocal,nv,nb),y(nlocal,na),stat=ios)
    if(ios/=0) call abort_world('allocate distributed X/RX buffers',7030)

    if(active_rank==0) then
      ! A new invocation always supersedes a completed or interrupted RbXc
      ! stage.  Replacing the manifest immediately records the current run as
      ! RUNNING; successful completion is marked at the end of the program.
      open(newunit=p,file='shg_rbxc_cv.info',status='replace',iostat=ios)
      call io_check(ios,'replace manifest')
      write(p,'(A)') 'status RUNNING'
      write(p,'(A)') 'sector CORE_FROM_VALENCE'
      write(p,'(A)') 'operator R_CV=-conjg(xes)_hole'
      write(p,'(A)') 'binary_layout raw_complex128_Fortran_ncon_nk_nalpha_no_padding'
      write(p,'(A)') 'units bohr_squared_per_hartree'
      write(p,'(A)') 'electron_frame unchanged_from_valence_k_plus_q'
      write(p,'(A)') 'spin_convention nbeta1_per_spin_diagonal_replication_nbeta4_explicit'
      write(p,'(A)') 'parallel_decomposition contiguous_conduction_band_blocks_one_k_at_a_time'
      write(p,'(A,1X,I0,1X,A,1X,I0)') 'mpi_ranks requested',world_size,'active',active_size
      do j=0,active_size-1
        call block_bounds(nc,active_size,j,nfirst,nlast)
        write(p,'(A,1X,I0,2(1X,I0))') 'mpi_conduction_block',j,nfirst,nlast
      enddo
      call block_bounds(nc,active_size,active_rank,nfirst,nlast)
      nlocal=nlast-nfirst+1
      write(p,'(A,1X,A)') 'middle_axis',bax
      write(p,'(A,1X,A)') 'right_axis',cax
      write(p,'(A,3(1X,I0))') 'dimensions ncon nk nalpha',nc,nk,na
      write(p,'(A,4(1X,I0))') 'input_brange',bands
      write(p,'(A,3(1X,I0))') 'kmesh',mesh
      write(p,'(A,3(1X,I0))') 'ZNL',znl
      write(p,'(A,1X,A,1X,I0)') 'site',trim(element),site
      write(p,'(A,4(1X,ES24.16))') 'eta_hartree_q_reduced',eta,q
      write(p,'(A,1X,I0)') 'energy_count',ne
      write(p,'(A)') 'columns index omega_Ha omega_eV norm file'
      flush(p)
    endif

    ! Restore every rank's own block after root generated the mapping table.
    call block_bounds(nc,active_size,active_rank,nfirst,nlast)
    nlocal=nlast-nfirst+1

    xsizes=(/nc,nv,nk,nb/)
    xsubs=(/nlocal,nv,1,nb/)
    ysizes=(/nc,nk,na/)
    ysubs=(/nlocal,1,na/)
    zero_offset=0_MPI_OFFSET_KIND
    output_bytes=16_MPI_OFFSET_KIND*int(nc,MPI_OFFSET_KIND) &
                               *int(nk,MPI_OFFSET_KIND)*int(na,MPI_OFFSET_KIND)

    do j=1,ne
      if(active_rank==0) then
        call require_bytes_root(infile(j),16_int64*nc*nv*nk*nb)
      endif
      call MPI_Barrier(active_comm,mpierr)
      call mpi_check(mpierr,'pre-I/O barrier')

      call MPI_File_open(active_comm,trim(infile(j)),MPI_MODE_RDONLY, &
                         MPI_INFO_NULL,xfh,mpierr)
      call mpi_check(mpierr,'open distributed X')
      call MPI_File_open(active_comm,trim(outfile(j)), &
                         ior(MPI_MODE_WRONLY,MPI_MODE_CREATE), &
                         MPI_INFO_NULL,yfh,mpierr)
      call mpi_check(mpierr,'create or replace distributed RX')
      ! Truncate or extend a pre-existing file to the exact current layout.
      ! Every element is subsequently rewritten by the collective k-point
      ! loop, so no data from an earlier or interrupted run can survive.
      call MPI_File_set_size(yfh,output_bytes,mpierr)
      call mpi_check(mpierr,'set RX file size')

      local_norm2=0.0_dp
      do k=1,nk
        xstarts=(/nfirst-1,0,k-1,0/)
        call MPI_Type_create_subarray(4,xsizes,xsubs,xstarts,MPI_ORDER_FORTRAN, &
                                      MPI_DOUBLE_COMPLEX,xtype,mpierr)
        call mpi_check(mpierr,'create X subarray type')
        call MPI_Type_commit(xtype,mpierr)
        call mpi_check(mpierr,'commit X subarray type')
        call MPI_File_set_view(xfh,zero_offset,MPI_DOUBLE_COMPLEX,xtype, &
                               'native',MPI_INFO_NULL,mpierr)
        call mpi_check(mpierr,'set X file view')
        call MPI_File_read_all(xfh,x,nlocal*nv*nb,MPI_DOUBLE_COMPLEX, &
                               status_mpi,mpierr)
        call mpi_check(mpierr,'read local X band slab')
        call MPI_Type_free(xtype,mpierr)
        call mpi_check(mpierr,'free X subarray type')

        bad_local=.not.all(ieee_is_finite(real(x))).or. &
                  .not.all(ieee_is_finite(aimag(x)))
        call MPI_Allreduce(bad_local,bad_global,1,MPI_LOGICAL,MPI_LOR, &
                           active_comm,mpierr)
        call mpi_check(mpierr,'check finite X')
        if(bad_global) call abort_world('nonfinite X',7032)

        y=cmplx(0.0_dp,0.0_dp,dp)
        ! OCEAN alpha order: hole spin, core m, electron spin (fastest).
        do h=1,2
          do m=1,nm
            do e=1,2
              a=((h-1)*nm+m-1)*2+e
              if(nb==1) then
                if(h/=e) cycle
                b=1
              else
                b=(h-1)*2+e
              endif
              y(:,a)=-matmul(x(:,:,b),d(:,k,m,min(h,ns)))
            enddo
          enddo
        enddo

        bad_local=.not.all(ieee_is_finite(real(y))).or. &
                  .not.all(ieee_is_finite(aimag(y)))
        call MPI_Allreduce(bad_local,bad_global,1,MPI_LOGICAL,MPI_LOR, &
                           active_comm,mpierr)
        call mpi_check(mpierr,'check finite RX')
        if(bad_global) call abort_world('nonfinite RX',7033)

        local_norm2=local_norm2+sum(real(conjg(y)*y,dp))

        ystarts=(/nfirst-1,k-1,0/)
        call MPI_Type_create_subarray(3,ysizes,ysubs,ystarts,MPI_ORDER_FORTRAN, &
                                      MPI_DOUBLE_COMPLEX,ytype,mpierr)
        call mpi_check(mpierr,'create RX subarray type')
        call MPI_Type_commit(ytype,mpierr)
        call mpi_check(mpierr,'commit RX subarray type')
        call MPI_File_set_view(yfh,zero_offset,MPI_DOUBLE_COMPLEX,ytype, &
                               'native',MPI_INFO_NULL,mpierr)
        call mpi_check(mpierr,'set RX file view')
        call MPI_File_write_all(yfh,y,nlocal*na,MPI_DOUBLE_COMPLEX, &
                                status_mpi,mpierr)
        call mpi_check(mpierr,'write local RX band slab')
        call MPI_Type_free(ytype,mpierr)
        call mpi_check(mpierr,'free RX subarray type')
      enddo

      call MPI_File_close(xfh,mpierr)
      call mpi_check(mpierr,'close distributed X')
      call MPI_File_close(yfh,mpierr)
      call mpi_check(mpierr,'close distributed RX')
      call MPI_Reduce(local_norm2,global_norm2,1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,active_comm,mpierr)
      call mpi_check(mpierr,'reduce RX norm')

      if(active_rank==0) then
        call require_bytes_root(outfile(j),16_int64*nc*nk*na)
        ynorm=sqrt(max(0.0_dp,global_norm2))
        write(p,'(A,1X,I0,3(1X,ES24.16),1X,A)') &
          'RbXc',j,omega(j),omegaev(j),ynorm,trim(outfile(j))
        flush(p)
      endif
    enddo

    if(active_rank==0) then
      write(p,'(A)') 'completion RBXC_CV_COMPLETE'
      close(p,iostat=ios)
      call io_check(ios,'close manifest')
    endif
    deallocate(x,y)
    call MPI_Comm_free(active_comm,mpierr)
    call mpi_check(mpierr,'free active communicator')
  endif

  deallocate(d,omega,omegaev,infile,outfile)
  call MPI_Finalize(mpierr)
  if(mpierr/=MPI_SUCCESS) error stop 'MPI_Finalize failed'

contains

  subroutine block_bounds(total,nblocks,irank,first,last)
    integer,intent(in) :: total,nblocks,irank
    integer,intent(out) :: first,last
    integer :: nbase,nextra,nmine
    nbase=total/nblocks
    nextra=mod(total,nblocks)
    if(irank<nextra) then
      nmine=nbase+1
      first=irank*(nbase+1)+1
    else
      nmine=nbase
      first=nextra*(nbase+1)+(irank-nextra)*nbase+1
    endif
    last=first+nmine-1
  end subroutine block_bounds

  subroutine io_check(code,where)
    integer,intent(in) :: code
    character(*),intent(in) :: where
    if(code/=0) call abort_world(trim(where)//' iostat failure',7100+abs(code))
  end subroutine io_check

  subroutine mpi_check(code,where)
    integer,intent(in) :: code
    character(*),intent(in) :: where
    if(code/=MPI_SUCCESS) then
      write(*,'(A,I0,2A,I0)') 'MPI rank ',world_rank,' failed at ',trim(where),code
      call MPI_Abort(MPI_COMM_WORLD,code,abort_ierr)
      error stop 1
    endif
  end subroutine mpi_check

  subroutine abort_world(message,code)
    character(*),intent(in) :: message
    integer,intent(in) :: code
    if(world_rank==0) write(*,'(A,I0,2A)') 'Fatal error ',code,': ',trim(message)
    call MPI_Abort(MPI_COMM_WORLD,code,abort_ierr)
    error stop 1
  end subroutine abort_world

  subroutine require_bytes_root(path,want)
    character(*),intent(in) :: path
    integer(int64),intent(in) :: want
    integer(int64) :: actual
    integer :: stat
    inquire(file=trim(path),size=actual,iostat=stat)
    call io_check(stat,'file size '//trim(path))
    if(actual/=want) then
      write(*,'(A,2(1X,I0))') trim(path)//' bytes actual/expected',actual,want
      call abort_world('file size mismatch',7110)
    endif
  end subroutine require_bytes_root

end program shg_rbxc_cv_mpi
