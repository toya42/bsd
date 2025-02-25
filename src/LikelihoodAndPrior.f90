!=====================================================================
! Module: LikelihoodAndPrior
! Purpose: Compute likelihood, normal prior, and log-posterior.
!=====================================================================
module LikelihoodAndPrior
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only : fp_kind
  use ModelFunctions
  use GlobalData
  use PriorProposal
  implicit none
contains

  !----------------------------------------------------------
  ! ComputeLambda: For each energy point, compute λ = I_inc * f_ratio(E, theta)
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
  ! ComputeLogLikelihood: Compute log-likelihood assuming Poisson statistics.
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogLikelihood(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), dimension(N) :: lambda
    integer(int32) :: i
    call ComputeLambda(theta, lambda)
    ComputeLogLikelihood = 0.0d0
    do i = 1, N
       if (lambda(i) > 0.0d0) then
          ComputeLogLikelihood = ComputeLogLikelihood + I_ab(i)*log(lambda(i)) - lambda(i)
       end if
    end do
  end function ComputeLogLikelihood

  !----------------------------------------------------------
  ! ComputeLogPrior: Evaluate log-prior
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogPrior(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind) :: log_prior
    integer :: k, j

    log_prior = 0.0d0
    ! Step function priors:
    !log_prior = log_prior + GammaLogPDF (theta%step%Ba,    alpha_Ba,    beta_Ba)
    log_prior = log_prior + GammaLogPDF(theta%step%Bb, alpha_Bb, beta_Bb)
    log_prior = log_prior + GammaLogPDF(theta%step%H,  alpha_H,  beta_H)
    log_prior = log_prior + NormalLogPDF(theta%step%E0, mu_E0, sigma_E0)
    log_prior = log_prior + UniformLogPDF(theta%step%Gamma, a_Gamma, b_Gamma)

    ! White-line priors:
    log_prior = log_prior + GammaLogPDF(theta%WL%A_WL, alpha_A_WL, beta_A_WL)
    log_prior = log_prior + NormalLogPDF(theta%WL%mu_WL, mu_mu_WL, sigma_mu_WL)
    log_prior = log_prior + UniformLogPDF(theta%WL%sigma_G_WL, a_sigma_G_WL, b_sigma_G_WL)
    log_prior = log_prior + UniformLogPDF(theta%WL%gamma_L_WL, a_gamma_L_WL, b_gamma_L_WL)

    ! Low-energy peak priors:
    do k = 1, size(theta%low)
       log_prior = log_prior + GammaLogPDF(theta%low(k)%A, alpha_low_A, beta_low_A)
       log_prior = log_prior + UniformLogPDF(theta%low(k)%mu, a_low_mu, b_low_mu)
       log_prior = log_prior + UniformLogPDF(theta%low(k)%sigma_G, a_low_sigma_G, b_low_sigma_G)
       log_prior = log_prior + UniformLogPDF(theta%low(k)%gamma_L, a_low_gamma_L, b_low_gamma_L)
       !! Optionally enforce constraint: low-energy mu < E0.
       !if (theta%low(k)%mu >= theta%step%E0) then
       !   ComputeLogPrior = -1.0d300
       !   return
       !end if
    end do

    ! High-energy peak priors:
    do j = 1, size(theta%high)
       log_prior = log_prior + GammaLogPDF(theta%high(j)%A, alpha_high_A, beta_high_A)
       log_prior = log_prior + UniformLogPDF(theta%high(j)%mu, a_high_mu, b_high_mu)
       log_prior = log_prior + UniformLogPDF(theta%high(j)%sigma_G, a_high_sigma_G, b_high_sigma_G)
       log_prior = log_prior + UniformLogPDF(theta%high(j)%gamma_L, a_high_gamma_L, b_high_gamma_L)
       ! Optionally enforce constraint: high-energy mu > E0.
       !if (theta%high(j)%mu <= theta%step%E0) then
       !   ComputeLogPrior = -1.0d300
       !   return
       !end if
    end do

    ComputeLogPrior = log_prior
  end function ComputeLogPrior

  !----------------------------------------------------------
  ! ComputeLogPosterior: Calculate log-posterior for given theta and replica inverse temperature beta.
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogPosterior(theta, beta_val)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), intent(in) :: beta_val
    ComputeLogPosterior = beta_val * ComputeLogLikelihood(theta) + ComputeLogPrior(theta)
  end function ComputeLogPosterior

end module LikelihoodAndPrior
