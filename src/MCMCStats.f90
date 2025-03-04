module MCMCStats
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  implicit none
contains

  !----------------------------------------------------------
  ! Function: ComputeAutocorrelation
  ! Purpose: Compute the autocorrelation function (ACF) of a time series x
  !          up to a maximum lag maxlag.
  !
  ! Input:
  !   x      - Real vector of MCMC samples.
  !   maxlag - Maximum lag to compute the ACF.
  !
  ! Output:
  !   acf    - Real vector of length (maxlag+1) containing ACF values,
  !            with acf(0) = 1.
  !----------------------------------------------------------
  function ComputeAutocorrelation(x, maxlag_in) result(acf)
    implicit none
    real(fp_kind), intent(in) :: x(:)
    integer, intent(in) :: maxlag_in
    real(fp_kind), allocatable, dimension(:) :: acf
    integer :: n, k, t, maxlag
    real(fp_kind) :: mean_x, var_x, numerator

    n = size(x)
    maxlag = min(maxlag_in, n-1)
    allocate(acf(0:maxlag))
    mean_x = sum(x) / n
    var_x = sum((x - mean_x)**2) / n

    acf(0) = 1.0d0
    do k = 1, maxlag
      numerator = 0.0d0
      do t = 1, n - k
        numerator = numerator + (x(t) - mean_x) * (x(t+k) - mean_x)
      end do
      !print *,k
      !print *,n-k 
      !print *,numerator
      !print *, var_x
      acf(k) = numerator / real(n - k, fp_kind) / var_x
    end do
  end function ComputeAutocorrelation

  !----------------------------------------------------------
  ! Function: ComputeIntegratedAutocorrelationTime
  ! Purpose: Compute the integrated autocorrelation time tau_int from an ACF vector.
  !
  ! Input:
  !   acf - Autocorrelation function vector, with acf(0)=1.
  !
  ! Output:
  !   tau_int - Integrated autocorrelation time.
  !----------------------------------------------------------
  function ComputeIntegratedAutocorrelationTime(acf) result(tau_int)
    implicit none
    real(fp_kind), intent(in) :: acf(:)
    real(fp_kind) :: tau_int
    integer :: k, maxlag

    maxlag = size(acf) - 1
    tau_int = 1.0d0
    do k = 1, maxlag
       tau_int = tau_int + 2.0d0 * acf(k)
    end do
  end function ComputeIntegratedAutocorrelationTime

  !----------------------------------------------------------
  ! Function: ComputeEffectiveSampleSize
  ! Purpose: Compute the effective sample size (ESS) for a time series x,
  !          given its integrated autocorrelation time tau_int.
  !
  ! Input:
  !   x       - Time series (vector of MCMC samples).
  !   tau_int - Integrated autocorrelation time.
  !
  ! Output:
  !   ess - Effective sample size, N / tau_int.
  !----------------------------------------------------------
  function ComputeEffectiveSampleSize(x, tau_int) result(ess)
    implicit none
    real(fp_kind), intent(in) :: x(:)
    real(fp_kind), intent(in) :: tau_int
    real(fp_kind) :: ess

    ess = size(x) / tau_int
  end function ComputeEffectiveSampleSize

end module MCMCStats
