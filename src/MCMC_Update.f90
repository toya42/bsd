!=====================================================================
! Module: MCMC_Update
! Purpose: Perform blockwise Metropolis-Hastings updates on ModelParameters.
!          Uses helper routines to extract and update parameter blocks,
!          and adapts proposal scale using the prior sigma.
!=====================================================================
module MCMC_Update
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use GlobalData
  use LikelihoodAndPrior
  use RNG
  use PriorProposal   ! Provides GetProposalSigma and GetPriorSigma
  use MultivariateProposal  ! New module for multivariate proposals.

  implicit none
contains

  !----------------------------------------------------------
  ! Function: TotalBlocks
  ! Returns the total number of parameter blocks.
  ! Block mapping:
  !   Block 1: Step parameters: [Ba, Bb, H, E0, Gamma] (5 parameters)
  !   Block 2: WL parameters: [A_WL, mu_WL, sigma_G_WL, gamma_L_WL] (4 parameters)
  !   Blocks 3 to 2+K1: Each low-energy peak (4 parameters each)
  !   Blocks (3+K1) to (2+K1+K2): Each high-energy peak (4 parameters each)
  !----------------------------------------------------------
  integer function TotalBlocks() result(total)
    implicit none
    total = 2 + K1 + K2
  end function TotalBlocks

  !----------------------------------------------------------
  ! Subroutine: ExtractBlock
  ! Extracts a block of parameters from the full ModelParameters.
  ! Input:
  !   theta     : full model parameters.
  !   block_id  : integer in 1 .. TotalBlocks()
  ! Output:
  !   block     : real array containing the parameters of the block.
  !----------------------------------------------------------
  subroutine ExtractBlock(theta, block_id, block)
    implicit none
    type(ModelParameters), intent(in) :: theta
    integer(int32), intent(in) :: block_id
    real(fp_kind), allocatable, intent(out) :: block(:)
    integer :: idx

    select case (block_id)
    case (1)
       allocate(block(5))
       block = [ theta%step%Ba, theta%step%Bb, theta%step%H, theta%step%E0, theta%step%Gamma ]
    case (2)
       allocate(block(4))
       block = [ theta%WL%A_WL, theta%WL%mu_WL, theta%WL%sigma_G_WL, theta%WL%gamma_L_WL ]
    case default
       if (block_id <= 2 + K1) then
          idx = block_id - 2
          allocate(block(4))
          block = [ theta%low(idx)%A, theta%low(idx)%mu, theta%low(idx)%sigma_G, theta%low(idx)%gamma_L ]
       else if (block_id <= 2 + K1 + K2) then
          idx = block_id - (2 + K1)
          allocate(block(4))
          block = [ theta%high(idx)%A, theta%high(idx)%mu, theta%high(idx)%sigma_G, theta%high(idx)%gamma_L ]
       else
          print *, "Error: block_id out of range in ExtractBlock."
          stop
       end if
    end select
  end subroutine ExtractBlock

  !----------------------------------------------------------
  ! Subroutine: UpdateBlock
  ! Writes new values from a candidate block into the full ModelParameters.
  ! Input:
  !   block_id  : integer block identifier.
  !   new_block : real array with new parameter values.
  ! In/out:
  !   theta     : full model parameters to be updated.
  !----------------------------------------------------------
  subroutine UpdateBlock(theta, block_id, new_block)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    integer(int32), intent(in) :: block_id
    real(fp_kind), intent(in) :: new_block(:)
    integer :: idx

    select case (block_id)
    case (1)
       if (size(new_block) /= 5) then
          print *, "Error: Block 1 size mismatch in UpdateBlock."
          stop
       end if
       theta%step%Ba    = new_block(1)
       theta%step%Bb    = new_block(2)
       theta%step%H     = new_block(3)
       theta%step%E0    = new_block(4)
       theta%step%Gamma = new_block(5)
    case (2)
       if (size(new_block) /= 4) then
          print *, "Error: Block 2 size mismatch in UpdateBlock."
          stop
       end if
       theta%WL%A_WL       = new_block(1)
       theta%WL%mu_WL       = new_block(2)
       theta%WL%sigma_G_WL  = new_block(3)
       theta%WL%gamma_L_WL  = new_block(4)
    case default
       if (block_id <= 2 + K1) then
          idx = block_id - 2
          if (size(new_block) /= 4) then
             print *, "Error: Low-energy block size mismatch in UpdateBlock."
             stop
          end if
          theta%low(idx)%A       = new_block(1)
          theta%low(idx)%mu      = new_block(2)
          theta%low(idx)%sigma_G = new_block(3)
          theta%low(idx)%gamma_L = new_block(4)
       else if (block_id <= 2 + K1 + K2) then
          idx = block_id - (2 + K1)
          if (size(new_block) /= 4) then
             print *, "Error: High-energy block size mismatch in UpdateBlock."
             stop
          end if
          theta%high(idx)%A       = new_block(1)
          theta%high(idx)%mu      = new_block(2)
          theta%high(idx)%sigma_G = new_block(3)
          theta%high(idx)%gamma_L = new_block(4)
       else
          print *, "Error: block_id out of range in UpdateBlock."
          stop
       end if
    end select
  end subroutine UpdateBlock

  !----------------------------------------------------------
  ! Subroutine: ProposeNew
  ! Generate a new candidate for a block by adding a normally distributed
  ! random perturbation to each parameter.
  ! Inputs:
  !   current_block : current parameter values.
  !   prior_sigma   : prior sigma values for each parameter in the block.
  ! Output:
  !   proposed_block: newly proposed parameter values.
  !----------------------------------------------------------
  subroutine ProposeNew(current_block, proposed_block, prior_sigma)
    implicit none
    real(fp_kind), intent(in) :: current_block(:)
    real(fp_kind), intent(in) :: prior_sigma(:)
    real(fp_kind), intent(out), allocatable :: proposed_block(:)
    integer(int32) :: i, block_size
    real(fp_kind) :: perturb, sigma_i

    block_size = size(current_block)
    if (size(prior_sigma) /= block_size) then
       print *, "Error: Size mismatch in ProposeNew between current_block and prior_sigma."
       stop
    end if

    allocate(proposed_block(block_size))
    proposed_block = current_block
    do i = 1, block_size
       call NormalRandom(perturb)  ! Draw perturbation from N(0,1)
       sigma_i = GetProposalSigma(prior_sigma(i))
       proposed_block(i) = current_block(i) + sigma_i * perturb
    end do
  end subroutine ProposeNew

  !----------------------------------------------------------
  ! Subroutine: BlockwiseMHUpdate
  ! Purpose: Perform a Metropolis-Hastings update on a given block 
  !          using a multivariate normal proposal with covariance.
  !----------------------------------------------------------
  subroutine BlockwiseMHUpdate(theta, block_id, beta_val)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    integer(int32), intent(in) :: block_id
    real(fp_kind), intent(in) :: beta_val
    real(fp_kind), allocatable :: current_block(:), proposed_block(:)
    real(fp_kind), allocatable :: prior_sigma(:)
    type(ModelParameters) :: theta_candidate
    real(fp_kind) :: logPost_current, logPost_proposed, delta, u

    ! Extract the current block from theta.
    call ExtractBlock(theta, block_id, current_block)
    
    ! Get the prior sigma vector for this block.
    ! Here we define a helper: GetPriorSigma (already defined in PriorProposal).
    prior_sigma = GetPriorSigma(block_id)
    
    ! Generate a candidate update using a multivariate normal proposal.
    call ProposeNew(current_block, proposed_block, prior_sigma)
    !call ProposeNewMV(current_block, prior_sigma, proposed_block)

    ! Create a candidate copy of theta.
    theta_candidate = theta
    call UpdateBlock(theta_candidate, block_id, proposed_block)

    ! Compute the log-posterior for current and candidate theta.
    logPost_current = ComputeLogPosterior(theta, beta_val)
    logPost_proposed = ComputeLogPosterior(theta_candidate, beta_val)
    delta = logPost_proposed - logPost_current

    call RandomUniform(u)
    if (u < exp(delta)) then
       ! Accept the candidate update.
       call UpdateBlock(theta, block_id, proposed_block)
    end if

    deallocate(current_block, proposed_block, prior_sigma)
  end subroutine BlockwiseMHUpdate

  !----------------------------------------------------------
  ! Subroutine: MCMC_UpdateReplica remains unchanged.
  !----------------------------------------------------------
  subroutine MCMC_UpdateReplica(theta, beta_val)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    real(fp_kind), intent(in) :: beta_val
    integer :: b, total_blocks
    total_blocks = TotalBlocks()
    do b = 1, total_blocks
       call BlockwiseMHUpdate(theta, b, beta_val)
    end do
  end subroutine MCMC_UpdateReplica

end module MCMC_Update
