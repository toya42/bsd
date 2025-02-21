!=====================================================================
! Module: ReplicaExchange
! Purpose: Implement odd–even replica exchange between neighboring replicas.
!          Exchanges are performed between adjacent replicas, alternating
!          odd-indexed and even-indexed pairs depending on the exchange_step.
!=====================================================================
module ReplicaExchange
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only : fp_kind
  use LikelihoodAndPrior    ! Provides ComputeLogLikelihood routine.
  use GlobalData            ! Contains beta array and other global parameters.
  use RNG                   ! Provides RandomUniform subroutine.
  implicit none
contains

  !----------------------------------------------------------
  ! Subroutine: SwapModelParameters
  ! Purpose: Swap two ModelParameters structures.
  !----------------------------------------------------------
  subroutine SwapModelParameters(theta1, theta2)
    implicit none
    type(ModelParameters), intent(inout) :: theta1, theta2
    type(ModelParameters) :: temp

    ! Simple swap using a temporary variable.
    temp = theta1
    theta1 = theta2
    theta2 = temp
  end subroutine SwapModelParameters

  !----------------------------------------------------------
  ! Subroutine: OddEvenExchange
  ! Purpose: Perform replica exchanges between neighboring replicas.
  !          Based on the exchange_step, alternate between:
  !            - Odd-indexed exchanges: pairs (1-2, 3-4, etc.)
  !            - Even-indexed exchanges: pairs (2-3, 4-5, etc.)
  !
  ! Input:
  !   theta_array   : array of ModelParameters for all replicas.
  !   beta_array    : array of inverse temperatures for each replica.
  !   exchange_step : an integer counter to determine odd/even exchange.
  !----------------------------------------------------------
  subroutine OddEvenExchange(theta_array, beta_array, exchange_step)
    implicit none
    type(ModelParameters), dimension(:), intent(inout) :: theta_array
    real(fp_kind), dimension(:), intent(in) :: beta_array
    integer(int32), intent(in) :: exchange_step  ! Exchange step counter
    real(fp_kind) :: logL_lower, logL_higher, delta_swap, r_swap, u
    integer(int32) :: l, num_replica

    ! Get the number of replicas in the array.
    num_replica = size(theta_array)

    ! Determine which pairs to attempt to swap based on exchange_step parity.
    if (mod(exchange_step, 2) == 1) then
       ! Odd-indexed exchange: swap replicas 1-2, 3-4, etc.
       do l = 1, num_replica - 1, 2
          logL_lower  = ComputeLogLikelihood(theta_array(l))
          logL_higher = ComputeLogLikelihood(theta_array(l+1))
          ! Compute the swap acceptance criterion:
          delta_swap = (beta_array(l+1) - beta_array(l)) * (logL_lower - logL_higher)
          if (0<delta_swap) then
            call SwapModelParameters(theta_array(l), theta_array(l+1))
          else if(delta_swap<-1000) then
            continue
          else 
            r_swap = exp(delta_swap)
            call RandomUniform(u)
            if (u < min(1.0d0, r_swap)) then
              call SwapModelParameters(theta_array(l), theta_array(l+1))
            end if
          end if
       end do
    else
       ! Even-indexed exchange: swap replicas 2-3, 4-5, etc.
       do l = 2, num_replica - 1, 2
          logL_lower  = ComputeLogLikelihood(theta_array(l))
          logL_higher = ComputeLogLikelihood(theta_array(l+1))
          delta_swap = (beta_array(l+1) - beta_array(l)) * (logL_lower - logL_higher)
          if (0<delta_swap) then
             call SwapModelParameters(theta_array(l), theta_array(l+1))
          else if(delta_swap<-1000) then
            continue
          else
            r_swap = exp(delta_swap)
            call RandomUniform(u)
            if (u < min(1.0d0, r_swap)) then
              call SwapModelParameters(theta_array(l), theta_array(l+1))
            end if
          end if
       end do
    end if

  end subroutine OddEvenExchange

end module ReplicaExchange
