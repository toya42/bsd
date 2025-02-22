!=====================================================================
! Module: PriorProposal
! Purpose: Define normal prior distributions for each parameter,
!          and compute proposal sigma as a fraction of the prior sigma.
!=====================================================================
module PriorProposal
  use, intrinsic :: iso_fortran_env, only: real64
  use precision, only: fp_kind
  use GlobalData, only: pi
  implicit none
  
  ! Prior parameters for step function:
  real(fp_kind), parameter :: mu_Ba    = 0.0d0, sigma_Ba    = 1.0d-3
  real(fp_kind), parameter :: mu_Bb    = 0.0d0, sigma_Bb    = 1.0d0
  real(fp_kind), parameter :: mu_H     = 3.2d4, sigma_H     = 1.0d4
  real(fp_kind), parameter :: mu_E0    = 2470.0d0, sigma_E0    = 5.0d2
  real(fp_kind), parameter :: mu_Gamma = 5.0d0, sigma_Gamma = 1.0d0
  
  ! Prior parameters for white-line (WL) function:
  real(fp_kind), parameter :: mu_A_WL       = 5.0d4, sigma_A_WL       = 1.0d4
  real(fp_kind), parameter :: mu_mu_WL      = 2470.0d0, sigma_mu_WL      = 2.0d2
  real(fp_kind), parameter :: mu_sigma_G_WL = 1.0d1, sigma_sigma_G_WL = 0.5d1
  real(fp_kind), parameter :: mu_gamma_L_WL = 1.0d1, sigma_gamma_L_WL = 0.5d1
  
  ! Prior parameters for each low-energy peak:
  real(fp_kind), parameter :: mu_low_A       = 3.0d4, sigma_low_A       = 1.0d4
  real(fp_kind), parameter :: mu_low_mu      = 2470.0d0, sigma_low_mu      = 2.0d2
  real(fp_kind), parameter :: mu_low_sigma_G = 1.0d1, sigma_low_sigma_G = 0.5d1
  real(fp_kind), parameter :: mu_low_gamma_L = 1.0d1, sigma_low_gamma_L = 0.5d1
  
  ! Prior parameters for each high-energy peak:
  real(fp_kind), parameter :: mu_high_A       = 1.0d4, sigma_high_A       = 0.1d4
  real(fp_kind), parameter :: mu_high_mu      = 2480.0d0, sigma_high_mu      = 2.0d2
  real(fp_kind), parameter :: mu_high_sigma_G = 1.0d1, sigma_high_sigma_G = 0.5d1
  real(fp_kind), parameter :: mu_high_gamma_L = 1.0d1, sigma_high_gamma_L = 0.5d1
  
  ! Proposal scaling factor.
  real(fp_kind), parameter :: c_proposal = 0.1d0
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
  ! Function: GetProposalSigma
  ! Purpose: Given the prior sigma for a parameter, return the proposal sigma.
  !----------------------------------------------------------
  real(fp_kind) function GetProposalSigma(sigma_prior)
    implicit none
    real(fp_kind), intent(in) :: sigma_prior
    GetProposalSigma = c_proposal * sigma_prior
  end function GetProposalSigma

end module PriorProposal