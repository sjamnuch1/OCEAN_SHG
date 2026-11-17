! Serial core <- valence hole vertex. No OCEAN modules or MPI required.
program shg_rbxc_cv
  use iso_fortran_env, only: dp=>real64, int32, int64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: nc,nv,nk,nb,l,ns,nvs,vskip,ne,nm,na
  integer :: u,ios,p,v,k,m,h,e,b,a,j,nproj,ntot,nspin
  integer(int32) :: header(3)
  integer(int64) :: expected
  character(4096) :: control,cksfile,melfile,infile,outfile,line
  character(1) :: bax,cax
  complex(dp), allocatable :: proj(:,:,:,:),mel(:,:),d(:,:,:,:),x(:,:,:,:),y(:,:,:)
  real(dp) :: re,im,omega,omegaev,eta,q(3),ynorm
  integer :: bands(4),mesh(3),znl(3),site
  character(32) :: element

  if(command_argument_count()/=1) error stop 'Usage: shg_rbxc_cv.x cv.in (run in prepared directory)'
  call get_command_argument(1,control)
  open(newunit=u,file=trim(control),status='old',iostat=ios)
  call check(ios,'open control')
  read(u,*,iostat=ios) nc,nv,nk,nb,l,ns,nvs,vskip,ne
  call check(ios,'control dimensions')
  if(min(nc,nv,nk,nvs,ne)<1.or.l<0.or.l>8) error stop 'Invalid dimensions'
  if(ns/=1.and.ns/=2) error stop 'Only collinear nspin=1,2 supported'
  if(nb/=1.and.nb/=4) error stop 'Only nbeta=1 or 4 supported'
  if(nb==1.and.ns/=1) error stop 'Scalar valence requires nspin=1'
  if(vskip<0.or.vskip+nv>nvs) error stop 'Valence slice outside projector file'
  read(u,*,iostat=ios) bax,cax
  call check(ios,'axes')
  read(u,*,iostat=ios) bands,mesh,znl,site,element
  call check(ios,'metadata')
  if(product(mesh)/=nk.or.znl(3)/=l) error stop 'Metadata dimensions disagree'
  read(u,*,iostat=ios) eta,q
  call check(ios,'energy/q metadata')
  read(u,'(A)',iostat=ios) cksfile
  call check(ios,'projector path')
  read(u,'(A)',iostat=ios) melfile
  call check(ios,'mels path')
  nm=2*l+1
  na=4*nm

  ! Native OCEAN parcksv: three int32 followed by complex128 payload.
  open(newunit=p,file=trim(cksfile),status='old',access='stream',form='unformatted',iostat=ios)
  call check(ios,'open parcksv')
  read(p,iostat=ios) header
  call check(ios,'parcksv header')
  nproj=header(1); ntot=header(2); nspin=header(3)
  if(nproj<1.or.int(ntot,int64)/=int(nvs,int64)*nk.or.nspin/=ns) &
    error stop 'Projector dimensions disagree'
  expected=12_int64+16_int64*nproj*nvs*nk*ns
  call bytes(cksfile,expected)
  allocate(proj(nproj,nvs,nk,ns),mel(nproj,nm),d(nv,nk,nm,ns))
  read(p,iostat=ios) proj
  call check(ios,'projector payload')
  close(p)
  if(.not.all(ieee_is_finite(real(proj))).or..not.all(ieee_is_finite(aimag(proj)))) &
    error stop 'Nonfinite projectors'
  open(newunit=p,file=trim(melfile),status='old',iostat=ios)
  call check(ios,'open mels')
  do m=1,nm
    do j=1,nproj
      read(p,*,iostat=ios) re,im
      call check(ios,'mels entry')
      if(.not.ieee_is_finite(re).or..not.ieee_is_finite(im)) error stop 'Nonfinite mels'
      mel(j,m)=cmplx(re,im,dp)
    enddo
  enddo
  do
    read(p,'(A)',iostat=ios) line
    if(ios<0) exit
    call check(ios,'mels trailing data')
    if(len_trim(line)>0) error stop 'Extra entries in mels'
  enddo
  close(p)
  ! OCEAN cksv_read: xes=sum_p mel(p,m)*proj(p,v,k,spin).
  ! SUM is intentional: Fortran DOT_PRODUCT conjugates its first argument.
  ! d=<v|r_b|g>=conjg(xes). No legacy RIXS intensity factor nspin.
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
  allocate(x(nc,nv,nk,nb),y(nc,nk,na))
  open(newunit=p,file='shg_rbxc_cv.info',status='unknown',iostat=ios)
  call check(ios,'create manifest (use fresh output directory)')
  write(p,'(A)') 'status RUNNING'
  write(p,'(A)') 'sector CORE_FROM_VALENCE'
  write(p,'(A)') 'operator R_CV=-conjg(xes)_hole'
  write(p,'(A)') 'binary_layout raw_complex128_Fortran_ncon_nk_nalpha_no_padding'
  write(p,'(A)') 'units bohr_squared_per_hartree'
  write(p,'(A)') 'electron_frame unchanged_from_valence_k_plus_q'
  write(p,'(A)') 'spin_convention nbeta1_per_spin_diagonal_replication_nbeta4_explicit'
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
  do j=1,ne
    read(u,*,iostat=ios) omega,omegaev
    call check(ios,'energy entry')
    read(u,'(A)',iostat=ios) infile
    call check(ios,'input file entry')
    read(u,'(A)',iostat=ios) outfile
    call check(ios,'output file entry')
    call bytes(infile,16_int64*nc*nv*nk*nb)
    open(newunit=e,file=trim(infile),status='old',access='stream',form='unformatted',iostat=ios)
    call check(ios,'open X')
    read(e,iostat=ios) x
    call check(ios,'read X')
    close(e)
    if(.not.all(ieee_is_finite(real(x))).or..not.all(ieee_is_finite(aimag(x)))) error stop 'Nonfinite X'
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
          do k=1,nk
            y(:,k,a)=-matmul(x(:,:,k,b),d(:,k,m,min(h,ns)))
          enddo
        enddo
      enddo
    enddo
    if(.not.all(ieee_is_finite(real(y))).or..not.all(ieee_is_finite(aimag(y)))) error stop 'Nonfinite RX'
    ynorm=norm2(reshape(real(y),[size(y)]))
    ynorm=hypot(ynorm,norm2(reshape(aimag(y),[size(y)])))
    open(newunit=e,file=trim(outfile),status='new',access='stream',form='unformatted',iostat=ios)
    call check(ios,'create RX')
    write(e,iostat=ios) y
    call check(ios,'write RX')
    close(e,iostat=ios)
    call check(ios,'close RX')
    write(p,'(A,1X,I0,3(1X,ES24.16),1X,A)') 'RbXc',j,omega,omegaev,ynorm,trim(outfile)
    flush(p)
  enddo
  write(p,'(A)') 'completion RBXC_CV_COMPLETE'
  close(p)
  close(u)
contains
  subroutine check(code,where)
    integer,intent(in)::code
    character(*),intent(in)::where
    if(code/=0) then
      write(*,*) trim(where),': iostat=',code
      error stop 1
    endif
  end subroutine
  subroutine bytes(path,want)
    character(*),intent(in)::path
    integer(int64),intent(in)::want
    integer(int64)::actual
    integer::stat
    inquire(file=trim(path),size=actual,iostat=stat)
    call check(stat,'file size '//trim(path))
    if(actual/=want) then
      write(*,*) trim(path),' bytes actual/expected ',actual,want
      error stop 'File size mismatch'
    endif
  end subroutine
end program
