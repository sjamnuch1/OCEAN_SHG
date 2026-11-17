! This file is part of OCEAN and is distributed under the University of
! Illinois/NCSA Open Source License. See License in the root directory.
module OCEAN_val_occupations
  use AI_kinds, only: DP
  use OCEAN_constants, only: eV2Hartree
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: fermi_occupation, pair_occupation, solve_chemical_potential, read_counting_spectrum, kB_hartree
  real(DP), parameter :: kB_hartree = 8.617333262145e-5_DP * eV2Hartree

contains

  elemental real(DP) function fermi_occupation(x) result(f)
    real(DP), intent(in) :: x  ! (energy - mu) / kT
    real(DP) :: z
    z = exp(-abs(x))
    if (x >= 0.0_DP) then
      f = z / (1.0_DP + z)
    else
      f = 1.0_DP / (1.0_DP + z)
    endif
  end function fermi_occupation

  elemental real(DP) function pair_occupation(ev, ec, mu, kt) result(df)
    real(DP), intent(in) :: ev, ec, mu, kt
    real(DP) :: xv, xc
    xv = (ev - mu) / kt
    xc = (ec - mu) / kt
    if (xc < 0.0_DP) then
      ! Subtract holes if both electronic occupations would be close to one.
      df = fermi_occupation(-xc) - fermi_occupation(-xv)
    else
      df = fermi_occupation(xv) - fermi_occupation(xc)
    endif
  end function pair_occupation

  subroutine read_counting_spectrum(nk, ns, energies, weights, ierr)
    integer, intent(in) :: nk, ns
    real(DP), allocatable, intent(out) :: energies(:,:,:), weights(:,:)
    integer, intent(out) :: ierr
    integer :: unit, nb, file_nk, file_ns, k, s, ios

    open(newunit=unit, file='val_occupation_energies.dat', status='old', action='read', iostat=ierr)
    if (ierr /= 0) return
    read(unit, *, iostat=ierr) nb, file_nk, file_ns
    if (ierr == 0) then
      if (nb < 1 .or. file_nk /= nk .or. file_ns /= ns) ierr = 1
    endif
    if (ierr == 0) allocate(energies(nb,nk,ns), weights(nk,ns), stat=ierr)
    if (ierr == 0) then
      do s = 1, ns
        do k = 1, nk
          read(unit, *, iostat=ierr) weights(k,s), energies(:,k,s)
          if (ierr /= 0) exit
        enddo
        if (ierr /= 0) exit
      enddo
    endif
    close(unit, iostat=ios)
    if (ierr /= 0) return
    if (.not. all(ieee_is_finite(energies)) .or. .not. all(ieee_is_finite(weights))) then
      ierr = 1
      return
    endif
    if (any(weights <= 0.0_DP) .or. any(abs(sum(weights, dim=1) - 1.0_DP) > 1.0e-10_DP)) ierr = 1
  end subroutine read_counting_spectrum

  subroutine solve_chemical_potential(energies, weights, nelec, kt, fixed, mu, residual, top, ierr)
    real(DP), intent(in) :: energies(:,:,:), weights(:,:), nelec, kt
    logical, intent(in) :: fixed
    real(DP), intent(inout) :: mu
    real(DP), intent(out) :: residual, top
    integer, intent(out) :: ierr
    real(DP) :: degeneracy, target, delta, lo, hi, mid, cmp, balance
    integer :: nb, ns, nref, iter

    ierr = 1
    if (.not. ieee_is_finite(nelec) .or. .not. ieee_is_finite(kt)) return
    if (kt <= 0.0_DP) return
    ns = size(energies, 3)
    nb = size(energies, 1)
    degeneracy = real(3-ns, DP) ! one spin channel represents two collinear spin states
    target = nelec / degeneracy
    if (target <= 0.0_DP .or. target >= real(nb, DP)*ns) return

    ! Use a filled-band reference to avoid subtracting two large counts in a gap.
    nref = nint(target/ns)
    delta = target - real(nref, DP)*ns
    if (.not. fixed) then
      lo = minval(energies) - 64.0_DP * kt
      hi = maxval(energies) + 64.0_DP * kt
      do iter = 1, 200
        mid = lo + (hi-lo)*0.5_DP
        if (mid <= lo .or. mid >= hi) exit
        if (delta == 0.0_DP .and. nref > 0 .and. nref < nb) then
          ! Logarithms preserve intrinsic carrier balance when both tails underflow.
          cmp = log_carriers((energies(nref+1:,:,:) - mid)/kt, weights) &
              - log_carriers((mid - energies(:nref,:,:))/kt, weights)
        else
          cmp = count_residual(mid)
        endif
        if (cmp == 0.0_DP) exit
        if (cmp < 0.0_DP) then
          lo = mid
        else
          hi = mid
        endif
      enddo
      mu = mid
    endif
    balance = count_residual(mu)
    residual = degeneracy * balance
    top = maxval(fermi_occupation((energies(nb,:,:) - mu)/kt))
    if (.not. fixed .and. abs(residual) > 1.0e-9_DP) return
    ierr = 0

  contains

    real(DP) function count_residual(chemical_potential) result(value)
      real(DP), intent(in) :: chemical_potential
      integer :: k, s
      value = -delta
      do s = 1, ns
        do k = 1, size(weights, 1)
          value = value + weights(k,s) * ( &
              sum(fermi_occupation((energies(nref+1:,k,s) - chemical_potential)/kt)) &
            - sum(fermi_occupation((chemical_potential - energies(:nref,k,s))/kt)))
        enddo
      enddo
    end function count_residual

  end subroutine solve_chemical_potential

  real(DP) function log_carriers(x, weights) result(value)
    real(DP), intent(in) :: x(:,:,:), weights(:,:)
    real(DP) :: logf(size(x,1),size(x,2),size(x,3)), shift
    integer :: k, s
    logf = -max(x, 0.0_DP) - log(1.0_DP + exp(-abs(x)))
    shift = maxval(logf)
    value = 0.0_DP
    do s = 1, size(x,3)
      do k = 1, size(x,2)
        value = value + weights(k,s)*sum(exp(logf(:,k,s) - shift))
      enddo
    enddo
    value = shift + log(value)
  end function log_carriers
end module OCEAN_val_occupations
