! Numerical checks independent of MPI, wavefunctions and the BSE kernel.
program test_occupations
  use AI_kinds, only: DP
  use OCEAN_val_occupations
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  real(DP) :: e(4,2,2), w(2,2), mu, residual, top, kt, target, counted, df, nan
  real(DP), parameter :: targets(3) = [0.125_DP, 2.0_DP, 3.375_DP]
  integer :: ierr, s, k, t

  w(:,1) = [0.25_DP, 0.75_DP]
  w(:,2) = w(:,1)
  e(:,1,1) = [-0.1_DP, 0.1_DP, 0.1_DP, 2.0_DP]
  e(:,2,1) = e(:,1,1)
  ! An intrinsic gap with two degenerate conduction edges gives
  ! mu = midgap - kT*log(2)/2 in the dilute-carrier limit.
  do t = 1, 3
    kt = kB_hartree * real(10**t, DP)
    mu = 0.0_DP
    call solve_chemical_potential(e(:,:,:1), w(:,:1), 2.0_DP, kt, .false., mu, residual, top, ierr)
    call check(ierr == 0, 'intrinsic solve')
    call check(abs(mu + kt*log(2.0_DP)/2) < 1.0e-12_DP, 'intrinsic mu including underflow')
    call check(abs(residual) < 1.0e-12_DP, 'intrinsic count')
  enddo

  ! Distinct spin spectra and nonuniform k weights; compare against a direct
  ! sum at temperatures where the ordinary sum is well conditioned.
  e(:,1,1) = [-0.07_DP, -0.01_DP, 0.03_DP, 1.0_DP]
  e(:,2,1) = [-0.05_DP, 0.00_DP, 0.08_DP, 1.2_DP]
  e(:,:,2) = e(:,:,1) + 0.02_DP
  kt = 0.015_DP
  do s = 1, 2
    do t = 1, 3
      target = targets(t)
      mu = 0.0_DP
      call solve_chemical_potential(e(:,:,:s), w(:,:s), target, kt, .false., mu, residual, top, ierr)
      call check(ierr == 0, 'fractional/spin solve')
      counted = 0.0_DP
      do k = 1, 2
        counted = counted + sum(w(k,:s) * sum(fermi_occupation((e(:,k,:s)-mu)/kt), dim=1))
      enddo
      counted = counted * (3-s)
      call check(abs(counted-target) < 1.0e-11_DP, 'independent electron count')
      call solve_chemical_potential(e(:,:,:s), w(:,:s), target+0.1_DP, kt, .true., mu, residual, top, ierr)
      call check(ierr == 0 .and. abs(residual+0.1_DP) < 1.0e-11_DP, 'fixed mu residual')
    enddo
  enddo
  df = pair_occupation(-41.0_DP, -40.0_DP, 0.0_DP, 1.0_DP)
  call check(abs(df/(exp(-40.0_DP)-exp(-41.0_DP))-1.0_DP) < 1.0e-12_DP, 'occupied-tail cancellation')
  call check(pair_occupation(0.1_DP, -0.1_DP, 0.0_DP, kt) < 0, 'reverse pair sign')

  nan = ieee_value(0.0_DP, ieee_quiet_nan)
  mu = nan
  call solve_chemical_potential(e, w, 2.0_DP, kt, .true., mu, residual, top, ierr)
  call check(ierr /= 0, 'reject nonfinite fixed mu')
  mu = 0.0_DP
  w = 0.0_DP
  call solve_chemical_potential(e, w, 2.0_DP, kt, .false., mu, residual, top, ierr)
  call check(ierr /= 0, 'reject invalid weights')
  print *, 'Occupation numerics passed: low T, fractional count, spin, tails and invalid inputs.'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print *, 'FAILED: ', label
      stop 1
    endif
  end subroutine
end program
