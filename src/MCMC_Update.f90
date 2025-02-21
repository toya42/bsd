!=====================================================================
! Module: MCMC_Update
! Purpose: Perform blockwise Metropolis-Hastings updates on ModelParameters.
!          Uses helper routines to extract and update parameter blocks.
!=====================================================================
module MCMC_Update
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use GlobalData
  use LikelihoodAndPrior
  use RNG
  implicit none
contains

  !----------------------------------------------------------
  ! Function: TotalBlocks
  ! Returns the total number of parameter blocks.
  ! Block mapping:
  !   Block 1: Step (5 parameters)
  !   Block 2: WL (4 parameters)
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
    real(fp_kind), intent(out), allocatable, dimension(:) :: block
    integer :: total, idx

    total = TotalBlocks()
    select case (block_id)
    case (1)
       ! Block 1: Step parameters: [Ba, Bb, H, E0, Gamma]
       allocate(block(5))
       block = [ theta%step%Ba, theta%step%Bb, theta%step%H, theta%step%E0, theta%step%Gamma ]
    case (2)
       ! Block 2: WL parameters: [A_WL, mu_WL, sigma_G_WL, gamma_L_WL]
       allocate(block(4))
       block = [ theta%WL%A_WL, theta%WL%mu_WL, theta%WL%sigma_G_WL, theta%WL%gamma_L_WL ]
    case default
       if (block_id <= 2 + K1) then
          ! Blocks 3 to 2+K1: low-energy peaks.
          idx = block_id - 2
          allocate(block(4))
          block = [ theta%low(idx)%A, theta%low(idx)%mu, theta%low(idx)%sigma_G, theta%low(idx)%gamma_L ]
       else if (block_id <= 2 + K1 + K2) then
          ! Blocks (3+K1) to (2+K1+K2): high-energy peaks.
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
       theta%step%Ba = new_block(1)
       theta%step%Bb = new_block(2)
       theta%step%H  = new_block(3)
       theta%step%E0 = new_block(4)
       theta%step%Gamma = new_block(5)
    case (2)
       if (size(new_block) /= 4) then
          print *, "Error: Block 2 size mismatch in UpdateBlock."
          stop
       end if
       theta%WL%A_WL = new_block(1)
       theta%WL%mu_WL = new_block(2)
       theta%WL%sigma_G_WL = new_block(3)
       theta%WL%gamma_L_WL = new_block(4)
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
  ! Generate a new candidate for a block by adding a small random perturbation.
  ! Uses the intrinsic random_number routine.
  !----------------------------------------------------------
  subroutine ProposeNew(current_block, proposed_block)
    implicit none
    real(fp_kind), intent(inout) :: current_block(:)
    real(fp_kind), intent(out), allocatable, dimension(:) :: proposed_block
    integer(int32) :: i, block_size
    real(fp_kind) :: perturb, sigma

    block_size = size(current_block)
    allocate(proposed_block(block_size))
    sigma = 0.1d0   ! Tuning parameter; adjust as needed.
    proposed_block = current_block
    do i = 1, block_size
       call RandomUniform(perturb)
       ! Map perturbation from [0,1] to [-sigma, sigma]
       proposed_block(i) = current_block(i) + sigma*(2.0d0*perturb - 1.0d0)
    end do
  end subroutine ProposeNew

  !----------------------------------------------------------
  ! Subroutine: BlockwiseMHUpdate
  ! Perform a Metropolis-Hastings update on a given block of parameters.
  !----------------------------------------------------------
  subroutine BlockwiseMHUpdate(theta, block_id, beta_val)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    integer(int32), intent(in) :: block_id
    real(fp_kind), intent(in) :: beta_val
    real(fp_kind), allocatable :: current_block(:), proposed_block(:)
    type(ModelParameters) :: theta_candidate
    real(fp_kind) :: logPost_current, logPost_proposed, delta, u

    ! Extract the current block from theta.
    call ExtractBlock(theta, block_id, current_block)

    ! Generate a candidate update.
    call ProposeNew(current_block, proposed_block)

    ! Copy current theta into theta_candidate.
    theta_candidate = theta

    ! Update theta_candidate with the proposed block.
    call UpdateBlock(theta_candidate, block_id, proposed_block)

    ! Compute log-posterior for current and candidate theta.
    logPost_current = ComputeLogPosterior(theta, beta_val)
    logPost_proposed = ComputeLogPosterior(theta_candidate, beta_val)
    delta = logPost_proposed - logPost_current

    call RandomUniform(u)
    !print *,delta
    if(0.0<delta) then
       ! Accept proposal: update theta with proposed block.
       call UpdateBlock(theta, block_id, proposed_block)
    else if(delta<-1000) then
       continue
    else if (u < exp(delta)) then
       ! Accept proposal: update theta with proposed block.
       call UpdateBlock(theta, block_id, proposed_block)
    else
       ! Reject proposal: do nothing.
    end if

    deallocate(current_block, proposed_block)
  end subroutine BlockwiseMHUpdate

  !----------------------------------------------------------
  ! Subroutine: MCMC_UpdateReplica
  ! Update all blocks of theta for a single replica.
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
