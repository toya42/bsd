!=====================================================================
! Module: ThermodynamicIntegration
! Purpose: Accumulate log-likelihoods and compute the marginal likelihood.
!=====================================================================
module ThermodynamicIntegration
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use GlobalData
  implicit none
  ! Fixed-size accumulators for each replica (size = L_rep_max).
  real(fp_kind), dimension(L_rep_max) :: SumLogL = 0.0d0  ! Sum of log-likelihoods per replica
  integer(int32), dimension(L_rep_max) :: Counts = 0       ! Counts of evaluations per replica
contains

  !----------------------------------------------------------
  ! Subroutine: UpdateLikelihoodAccumulator
  ! Purpose: Update the accumulator for a given replica.
  ! Inputs:
  !   replica_index : Index of the replica (1 to L_rep)
  !   currentLogL   : Log-likelihood value from the current MCMC update
  !----------------------------------------------------------
  subroutine UpdateLikelihoodAccumulator(replica_index, currentLogL)
    implicit none
    integer(int32), intent(in) :: replica_index
    real(fp_kind), intent(in) :: currentLogL
    SumLogL(replica_index) = SumLogL(replica_index) + currentLogL
    Counts(replica_index) = Counts(replica_index) + 1
  end subroutine UpdateLikelihoodAccumulator

  !----------------------------------------------------------
  ! Subroutine: ComputeAverageLogLikelihood
  ! Purpose: Compute the average log-likelihood for each replica.
  ! Output:
  !   avgLogL : Array of average log-likelihoods (size must equal L_rep)
  !----------------------------------------------------------
  subroutine ComputeAverageLogLikelihood(avgLogL)
    implicit none
    real(fp_kind), dimension(:), intent(out) :: avgLogL
    integer(int32) :: ll
    if (size(avgLogL) /= L_rep) then
       print *, "Error: avgLogL size mismatch. Expected size:", L_rep
       stop
    end if
    do ll = 1, L_rep
       if (Counts(ll) > 0) then
          avgLogL(ll) = SumLogL(ll) / real(Counts(ll), fp_kind)
       else
          avgLogL(ll) = -1.0d300
       end if
    end do
  end subroutine ComputeAverageLogLikelihood

  !----------------------------------------------------------
  ! Function: ComputeMarginalLikelihood
  ! Purpose: Compute the marginal likelihood logZ using the trapezoidal rule.
  ! Inputs:
  !   avgLogL    : Array of average log-likelihoods for each replica.
  !   beta_array : Array of inverse temperatures for the replicas.
  ! Returns:
  !   logZ       : Approximated marginal likelihood (log Z).
  !----------------------------------------------------------
  real(fp_kind) function ComputeMarginalLikelihood(avgLogL, beta_array)
    implicit none
    real(fp_kind), dimension(:), intent(in) :: avgLogL
    real(fp_kind), dimension(:), intent(in) :: beta_array
    integer(int32) :: ll, nBeta
    nBeta = size(beta_array)
    ComputeMarginalLikelihood = 0.0d0
    do ll = 1, nBeta - 1
       ComputeMarginalLikelihood = ComputeMarginalLikelihood + &
            (beta_array(ll+1) - beta_array(ll)) * 0.5d0 * (avgLogL(ll) + avgLogL(ll+1))
    end do
  end function ComputeMarginalLikelihood

  !----------------------------------------------------------
  ! Function: ComputeBayesFreeEnergy
  ! Purpose: Compute the Bayes free energy F = -logZ.
  !----------------------------------------------------------
  real(fp_kind) function ComputeBayesFreeEnergy(logZ)
    implicit none
    real(fp_kind), intent(in) :: logZ
    ComputeBayesFreeEnergy = - logZ
  end function ComputeBayesFreeEnergy

end module ThermodynamicIntegration
