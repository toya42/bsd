!=====================================================================
! Module: PriorProposal
! Purpose: Define normal prior distributions for each parameter,
!          and compute proposal sigma as a fraction of the prior sigma.
!=====================================================================
module PriorProposal
  use, intrinsic :: iso_fortran_env
  use precision, only: fp_kind
  use GlobalData, only: pi, K1, K2,L_rep_max
  implicit none
  
  ! Prior parameters for step function:
  real(fp_kind), parameter :: mu_Ba = 0.0d0, sigma_Ba = 2.0d-9
  real(fp_kind), parameter :: mu_Bb = 0.0d0, sigma_Bb = 2.0d-9
  real(fp_kind), parameter :: alpha_H  = 10.0d0, beta_H  = 5.0d-7, sigma_H  = 1.0d-6
  real(fp_kind), parameter :: a_E0 = 2400.0d0, b_E0 = 2560.0d0, sigma_E0 = 2.0d0
  real(fp_kind), parameter :: a_Gamma = 1.0d-5, b_Gamma = 1.0d1, sigma_Gamma = 5.0d-1
  
  ! Prior parameters for white-line (WL) function:
  real(fp_kind), parameter :: alpha_A_WL = 10.0d0, beta_A_WL = 5.0d-7, sigma_A_WL = 1.0d-5
  real(fp_kind), parameter :: mu_mu_WL = 0.0d0, sigma_mu_WL = 2.0d0
  !real(fp_kind), parameter :: a_mu_WL = -20.0d0, b_mu_WL = 20.0d0, sigma_mu_WL      = 5.0d0
  real(fp_kind), parameter :: a_sigma_G_WL = 1.0d-5, b_sigma_G_WL = 1.0d2, sigma_sigma_G_WL = 1.0d-0
  real(fp_kind), parameter :: a_gamma_L_WL = 1.0d-5, b_gamma_L_WL = 1.0d2, sigma_gamma_L_WL = 5.0d-1
  
  ! Prior parameters for each low-energy peak:
  real(fp_kind), parameter :: alpha_low_A = 10.0d0, beta_low_A = 5.0d-7, sigma_low_A = 1.0d-5
  real(fp_kind), parameter :: a_low_mu      = 0.0d0, b_low_mu = 30.0d0, sigma_low_mu      = 5.0d0
  real(fp_kind), parameter :: a_low_sigma_G = 1.0d-5, b_low_sigma_G = 1.0d2, sigma_low_sigma_G = 1.0d-0
  real(fp_kind), parameter :: a_low_gamma_L = 1.0d-5, b_low_gamma_L = 1.0d2, sigma_low_gamma_L = 5.0d-1
  
  ! Prior parameters for each high-energy peak:
  real(fp_kind), parameter :: alpha_high_A = 10.0d0, beta_high_A = 5.0d-7, sigma_high_A = 4.0d-6
  real(fp_kind), parameter :: a_high_mu      = 0.0d0, b_high_mu = 200.0d0, sigma_high_mu      = 1.0d1
  real(fp_kind), parameter :: a_high_sigma_G = 1.0d-5, b_high_sigma_G = 1.0d2, sigma_high_sigma_G = 1.0d-0
  real(fp_kind), parameter :: a_high_gamma_L = 1.0d-5, b_high_gamma_L = 1.0d2, sigma_high_gamma_L = 5.0d-1
  
  ! Proposal scaling factor.
  real(fp_kind),allocatable, dimension(:,:) :: c_proposal
contains

  !----------------------------------------------------------
  ! Function: NormalLogPDF
  ! Purpose: Compute the log density of a normal distribution at x.
  !----------------------------------------------------------
  real(fp_kind) function NormalLogPDF(x, mu, sigma)
    implicit none
    real(fp_kind), intent(in) :: x, mu, sigma
    NormalLogPDF = -0.5d0 * ((x - mu)/sigma)**2 - log(sigma * sqrt(2.0d0*pi))
  end function NormalLogPDF

  !----------------------------------------------------------
  ! Function: GammaLogPDF
  ! Purpose: Compute the log density of a Gamma distribution at x.
  !          Parameterization: shape (alpha) and rate (beta).
  !          Valid for x > 0; returns a very low value if x <= 0.
  !----------------------------------------------------------
  real(fp_kind) function GammaLogPDF(x, alpha, beta)
    implicit none
    real(fp_kind), intent(in) :: x, alpha, beta
    if (x <= 0.0d0) then
       GammaLogPDF = -huge(1.0d0)
    else
       ! Using: log f(x) = alpha*log(beta) - log Gamma(alpha) + (alpha-1)*log x - beta*x
       GammaLogPDF = alpha * log(beta) - log_gamma(alpha) + (alpha - 1.0d0) * log(x) - beta * x
    end if
  end function GammaLogPDF

  !----------------------------------------------------------
  ! Function: UniformLogPDF
  ! Purpose: Compute the log density of a Uniform distribution over [a, b].
  !          Returns -log(b - a) if x is within [a, b]; otherwise returns a very low value.
  !----------------------------------------------------------
  real(fp_kind) function UniformLogPDF(x, a, b)
    implicit none
    real(fp_kind), intent(in) :: x, a, b
    if (x < a .or. x > b) then
       UniformLogPDF = -huge(1.0d0)
    else
       UniformLogPDF = -log(b - a)
    end if
  end function UniformLogPDF


  !----------------------------------------------------------
  ! Function: GetProposalSigma
  ! Purpose: Given the prior sigma for a parameter, return the proposal sigma.
  !----------------------------------------------------------
  real(fp_kind) function GetProposalSigma(sigma_prior,b,l)
    implicit none
    integer(int32), intent(in) :: b,l
    real(fp_kind), intent(in) :: sigma_prior
    GetProposalSigma = c_proposal(b,l) * sigma_prior
  end function GetProposalSigma


  !----------------------------------------------------------
  ! Function: GetPriorSigma
  ! Purpose: Return an array of prior sigma values for a given parameter block.
  ! Input:
  !   block_id : integer identifier for the block.
  !              Block 1: Step parameters [Ba, Bb, H, E0, Gamma] (5 values)
  !              Block 2: WL parameters [A_WL, mu_WL, sigma_G_WL, gamma_L_WL] (4 values)
  !              Blocks 3 to 2+K1: Low-energy peak parameters (4 values each)
  !              Blocks (3+K1) to (2+K1+K2): High-energy peak parameters (4 values each)
  ! Output:
  !   prior_sigma : allocatable real array containing the prior sigma values.
  !----------------------------------------------------------
  function GetPriorSigma(block_id) result(prior_sigma)
    implicit none
    integer, intent(in) :: block_id
    real(fp_kind), allocatable :: prior_sigma(:)
    integer :: idx

    select case (block_id)
    case (1)
       allocate(prior_sigma(5))
       prior_sigma = [ sigma_Ba, sigma_Bb, sigma_H, sigma_E0, sigma_Gamma ]
    case (2)
       allocate(prior_sigma(4))
       prior_sigma = [ sigma_A_WL, sigma_mu_WL, sigma_sigma_G_WL, sigma_gamma_L_WL ]
    case default
       if (block_id <= 2 + K1) then
          allocate(prior_sigma(4))
          prior_sigma = [ sigma_low_A, sigma_low_mu, sigma_low_sigma_G, sigma_low_gamma_L ]
       else if (block_id <= 2 + K1 + K2) then
          allocate(prior_sigma(4))
          prior_sigma = [ sigma_high_A, sigma_high_mu, sigma_high_sigma_G, sigma_high_gamma_L ]
       else
          print *, "Error: block_id out of range in GetPriorSigma."
          stop
       end if
    end select
  end function GetPriorSigma

  function GetPriorSigmaFull() result(prior_sigma)
    implicit none
    real(fp_kind), allocatable :: prior_sigma(:)
    integer :: idx,i

    allocate(prior_sigma(5+4+4*K1+4*K2))

    prior_sigma(1) = sigma_Ba
    prior_sigma(2) = sigma_Bb
    prior_sigma(3) = sigma_H
    prior_sigma(4) = sigma_E0
    prior_sigma(5) = sigma_Gamma

    prior_sigma(6) = sigma_A_WL
    prior_sigma(7) = sigma_mu_WL
    prior_sigma(8) = sigma_sigma_G_WL
    prior_sigma(9) = sigma_gamma_L_WL

    i = 10
    do idx=1,K1
      prior_sigma(i)   = sigma_low_A
      prior_sigma(i+1) = sigma_low_mu
      prior_sigma(i+2) = sigma_low_sigma_G
      prior_sigma(i+3) = sigma_low_gamma_L
      i=i+4
    end do
    do idx=1,K2
      prior_sigma(i)   = sigma_high_A
      prior_sigma(i+1) = sigma_high_mu
      prior_sigma(i+2) = sigma_high_sigma_G
      prior_sigma(i+3) = sigma_high_gamma_L
      i=i+4
    end do

  end function GetPriorSigmaFull


end module PriorProposal
