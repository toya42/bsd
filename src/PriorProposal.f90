!=====================================================================
! Module: PriorProposal
! Purpose: Define normal prior distributions for each parameter,
!          and compute proposal sigma as a fraction of the prior sigma.
!=====================================================================
module PriorProposal
  use, intrinsic :: iso_fortran_env, only: real64
  use precision, only: fp_kind
  use GlobalData, only: pi, K1, K2
  implicit none
  
  ! Prior parameters for step function:
  real(fp_kind), parameter :: mu_Ba    = 0.0d0, sigma_Ba    = 1.0d-8
  real(fp_kind), parameter :: mu_Bb    = 2000.0d0*2.0d-10, sigma_Bb    = 1.0d-8
  real(fp_kind), parameter :: mu_H     = 0.0d0, sigma_H     = 1.0d-8
  real(fp_kind), parameter :: mu_E0    = 2480.0d0, sigma_E0    = 5.0d0
  real(fp_kind), parameter :: mu_Gamma = 1.0d-2, sigma_Gamma = 1.0d-4
  
  ! Prior parameters for white-line (WL) function:
  real(fp_kind), parameter :: mu_A_WL       = 0.0d0, sigma_A_WL       = 0.1d-8
  real(fp_kind), parameter :: mu_mu_WL      = 2480.0d0, sigma_mu_WL      = 5.0d0
  real(fp_kind), parameter :: mu_sigma_G_WL = 2.0d-1, sigma_sigma_G_WL = 0.5d-3
  real(fp_kind), parameter :: mu_gamma_L_WL = 2.0d-1, sigma_gamma_L_WL = 0.5d-3
  
  ! Prior parameters for each low-energy peak:
  real(fp_kind), parameter :: mu_low_A       = 1.5d4, sigma_low_A       = 0.1d3
  real(fp_kind), parameter :: mu_low_mu      = 2480.0d0, sigma_low_mu      = 2.0d1
  real(fp_kind), parameter :: mu_low_sigma_G = 5.0d-1, sigma_low_sigma_G = 0.5d-3
  real(fp_kind), parameter :: mu_low_gamma_L = 5.0d-1, sigma_low_gamma_L = 0.5d-3
  
  ! Prior parameters for each high-energy peak:
  real(fp_kind), parameter :: mu_high_A       = 1.0d4, sigma_high_A       = 0.1d3
  real(fp_kind), parameter :: mu_high_mu      = 2480.0d0, sigma_high_mu      = 2.0d1
  real(fp_kind), parameter :: mu_high_sigma_G = 1.0d1, sigma_high_sigma_G = 0.5d-2
  real(fp_kind), parameter :: mu_high_gamma_L = 1.0d1, sigma_high_gamma_L = 0.5d-2
  
  ! Proposal scaling factor.
  real(fp_kind), parameter :: c_proposal = 1.0d0
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

end module PriorProposal
