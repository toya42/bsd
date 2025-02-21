!=====================================================================
! Module: LikelihoodAndPrior
! Purpose: Compute likelihood, prior, and log-posterior.
!          Uses fixed-size arrays (dimension(N)) for speed.
!=====================================================================
module LikelihoodAndPrior
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only : fp_kind
  use ModelFunctions
  use GlobalData
  implicit none
contains

  !----------------------------------------------------------
  ! ComputeLambda:
  ! For each energy point, compute λ = I_inc * f_ratio(E, theta)
  ! Here, lambda is an array of fixed size N.
  !----------------------------------------------------------
  subroutine ComputeLambda(theta, lambda)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), dimension(N), intent(out) :: lambda
    integer(int32) :: i

    do i = 1, N
       lambda(i) = I_inc(i) * f_ratio(E(i), theta)
    end do
  end subroutine ComputeLambda

  !----------------------------------------------------------
  ! ComputeLogLikelihood:
  ! Compute the log-likelihood assuming Poisson statistics.
  ! For each measurement point i, the contribution is:
  !     I_ab(i)*log(lambda(i)) - lambda(i)
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogLikelihood(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), dimension(N) :: lambda
    integer(int32) :: i

    ! Compute the expected counts (lambda) for all energy points.
    call ComputeLambda(theta, lambda)

    ComputeLogLikelihood = 0.0d0
    do i = 1, N
       if (lambda(i) > 0.0d0) then
          ComputeLogLikelihood = ComputeLogLikelihood + I_ab(i)*log(lambda(i)) - lambda(i)
       end if
    end do
  end function ComputeLogLikelihood

  !----------------------------------------------------------
  ! ComputeLogPrior:
  ! Evaluate the log-prior probability for the model parameters.
  ! In this example, we assume a flat prior with constraints:
  !    - For low-energy peaks: peak center must be less than E0.
  !    - For high-energy peaks: peak center must be greater than E0.
  ! If a constraint is violated, a very large negative value is returned.
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogPrior(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    integer(int32) :: k, j

    ComputeLogPrior = 0.0d0
    ! Check low-energy peak constraints.
    do k = 1, size(theta%low)
       if (theta%low(k)%mu >= theta%step%E0) then
          ComputeLogPrior = -1.0d300
          return
       end if
    end do
    ! Check high-energy peak constraints.
    do j = 1, size(theta%high)
       if (theta%high(j)%mu <= theta%step%E0) then
          ComputeLogPrior = -1.0d300
          return
       end if
    end do
    ! Flat prior otherwise.
    ComputeLogPrior = 0.0d0
  end function ComputeLogPrior

  !----------------------------------------------------------
  ! ComputeLogPosterior:
  ! Calculate the log-posterior for a given set of model parameters,
  ! taking into account an inverse temperature (beta) for replica exchange.
  ! The log-posterior is defined as:
  !    logPosterior = beta * log-likelihood + log-prior
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogPosterior(theta, beta_val)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), intent(in) :: beta_val

    ComputeLogPosterior = beta_val * ComputeLogLikelihood(theta) + ComputeLogPrior(theta)
  end function ComputeLogPosterior

end module LikelihoodAndPrior
