module precision
  use, intrinsic :: iso_fortran_env
  integer(int32), parameter :: fp_kind = real64 ! or real32
end module precision
!=====================================================================
! Module: GlobalData
! Purpose: Store global input arrays and configuration parameters.
!=====================================================================
module GlobalData
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  implicit none
  ! Parameters for the data size and MCMC settings
  integer(int32), parameter :: N = 100            ! Number of energy points
  integer(int32), parameter :: T_iter = 10000     ! Total number of MCMC iterations
  integer(int32), parameter :: T_burn = 2000      ! Burn-in iterations
  integer(int32), parameter :: m  = 4             ! Number of parameter blocks
  integer(int32), parameter :: L_rep = 4          ! Number of replicas

  ! Input arrays: energy and photon counts
  real(fp_kind), dimension(N) :: E       ! Energy array
  real(fp_kind), dimension(N) :: I_inc   ! Incident counts at each energy
  real(fp_kind), dimension(N) :: I_ab    ! Observed absorbed counts

  ! Model-related parameters: numbers of peaks in different domains
  integer(int32), parameter :: K1 = 4   ! Number of low-energy peaks
  integer(int32), parameter :: K2 = 4   ! Number of high-energy peaks

  ! Replica exchange settings: inverse temperatures (beta)
  real(fp_kind), dimension(L_rep) :: beta  ! For example, beta(1) < beta(2) < ... < beta(L_rep)=1

end module GlobalData

!=====================================================================
! Module: ModelFunctions
! Purpose: Define model functions for the absorption spectrum.
!=====================================================================
module ModelFunctions
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  implicit none

  ! Derived type for the step (edge) function parameters
  type :: ParametersStep
     real(fp_kind) :: Ba,Bb  ! Baseline intensity (assume linear function) Ba*E+Bb
     real(fp_kind) :: H      ! Step height
     real(fp_kind) :: E0     ! Edge energy (absorption-edge position)
     real(fp_kind) :: Gamma  ! Broadening parameter
  end type ParametersStep

  ! Derived type for the white-line (peak) parameters
  type :: ParametersWL
     real(fp_kind) :: A_WL       ! Peak amplitude
     real(fp_kind) :: mu_WL      ! Peak center (near E0)
     real(fp_kind) :: sigma_G_WL ! Gaussian width
     real(fp_kind) :: gamma_L_WL ! Lorentzian half-width
  end type ParametersWL

  ! Derived type for individual Voigt peak parameters (used in low and high energy domains)
  type :: ParametersPeak
     real(fp_kind) :: A       ! Peak amplitude
     real(fp_kind) :: mu      ! Peak center energy
     real(fp_kind) :: sigma_G ! Gaussian width
     real(fp_kind) :: gamma_L ! Lorentzian half-width
  end type ParametersPeak

  ! Combined type for all model parameters
  type :: ModelParameters
     type(ParametersStep) :: step
     type(ParametersWL)   :: WL
     type(ParametersPeak), dimension(K1) :: low   ! Low-energy peaks
     type(ParametersPeak), dimension(K2) :: high  ! High-energy peaks
  end type ModelParameters

contains

  !----------------------------------------------------------
  ! Dummy Voigt profile function (to be replaced by an accurate routine)
  !----------------------------------------------------------
  real(fp_kind) function Voigt(x, sigma, gamma)
    implicit none
    real(fp_kind), intent(in) :: x, sigma, gamma
    ! Simple approximation: average of Gaussian and Lorentzian components.
    Voigt = 0.5d0 * exp( - (x**2) / (2.0d0 * sigma**2) ) / ( sigma * sqrt(2.0d0*acos(-1.0d0)) ) &
            + 0.5d0 * ( gamma / (acos(-1.0d0) * (x**2 + gamma**2)) )
  end function Voigt

  !----------------------------------------------------------
  ! f_step: smoothed step function modeling the absorption edge.
  ! f_step = Ba*E+Bb + H * Φ((E - E0)/Gamma), where Φ(x) = 0.5 + (1/π) atan(x)
  !----------------------------------------------------------
  real(fp_kind) function f_step(E, theta_step)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ParametersStep), intent(in) :: theta_step
    real(fp_kind) :: x
    x = (E - theta_step%E0) / theta_step%Gamma
    f_step = theta_step%Ba*x + theta_step%Bb &
            + theta_step%H * (0.5d0 + (1.0d0/acos(-1.0d0))*atan(x))
  end function f_step

  !----------------------------------------------------------
  ! f_WL: white line peak using a Voigt profile.
  !----------------------------------------------------------
  real(fp_kind) function f_WL(E, theta_WL)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ParametersWL), intent(in) :: theta_WL
    f_WL = theta_WL%A_WL * Voigt(E - theta_WL%mu_WL, theta_WL%sigma_G_WL, theta_WL%gamma_L_WL)
  end function f_WL

  !----------------------------------------------------------
  ! f_low: Sum of Voigt profiles for low-energy peaks.
  !----------------------------------------------------------
  real(fp_kind) function f_low(E, theta_low)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ParametersPeak), dimension(:), intent(in) :: theta_low
    integer(int32) :: k
    f_low = 0.0d0
    do k = 1, size(theta_low)
       f_low = f_low + theta_low(k)%A * Voigt(E - theta_low(k)%mu, theta_low(k)%sigma_G, theta_low(k)%gamma_L)
    end do
  end function f_low

  !----------------------------------------------------------
  ! f_high: Sum of Voigt profiles for high-energy peaks.
  !----------------------------------------------------------
  real(fp_kind) function f_high(E, theta_high)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ParametersPeak), dimension(:), intent(in) :: theta_high
    integer(int32) :: j
    f_high = 0.0d0
    do j = 1, size(theta_high)
       f_high = f_high + theta_high(j)%A * Voigt(E - theta_high(j)%mu, theta_high(j)%sigma_G, theta_high(j)%gamma_L)
    end do
  end function f_high

  !----------------------------------------------------------
  ! f_ratio: Full absorption ratio function as sum of components.
  !----------------------------------------------------------
  real(fp_kind) function f_ratio(E, theta)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ModelParameters), intent(in) :: theta
    f_ratio = f_step(E, theta%step) + f_WL(E, theta%WL) + f_low(E, theta%low) + f_high(E, theta%high)
  end function f_ratio

end module ModelFunctions

!=====================================================================
! Module: LikelihoodAndPrior
! Purpose: Compute likelihood, prior, and log-posterior.
!=====================================================================
module LikelihoodAndPrior
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use ModelFunctions
  use GlobalData
  implicit none
contains

  !----------------------------------------------------------
  ! ComputeLambda: For each energy point, compute λ = I_inc * f_ratio.
  !----------------------------------------------------------
  subroutine ComputeLambda(theta, lambda)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), dimension(:), intent(out) :: lambda
    integer(int32) :: i
    do i = 1, size(lambda)
       lambda(i) = I_inc(i) * f_ratio(E(i), theta)
    end do
  end subroutine ComputeLambda

  !----------------------------------------------------------
  ! ComputeLogLikelihood: Compute log-likelihood assuming Poisson statistics.
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogLikelihood(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), allocatable :: lambda(:)
    integer(int32) :: i
    allocate(lambda(N))
    call ComputeLambda(theta, lambda)
    ComputeLogLikelihood = 0.0d0
    do i = 1, N
       if (lambda(i) > 0.0d0) then
          ComputeLogLikelihood = ComputeLogLikelihood + I_ab(i)*log(lambda(i)) - lambda(i)
       end if
    end do
    deallocate(lambda)
  end function ComputeLogLikelihood

  !----------------------------------------------------------
  ! ComputeLogPrior: Evaluate log-prior.
  ! In this example, we assume flat priors with constraints:
  !    Low-energy peak centers must be less than E0.
  !    High-energy peak centers must be greater than E0.
  !----------------------------------------------------------
  real(fp_kind) function ComputeLogPrior(theta)
    implicit none
    type(ModelParameters), intent(in) :: theta
    integer(int32) :: k, j
    ComputeLogPrior = 0.0d0
    do k = 1, size(theta%low)
       if (theta%low(k)%mu >= theta%step%E0) then
          ComputeLogPrior = -1.0d300
          return
       end if
    end do
    do j = 1, size(theta%high)
       if (theta%high(j)%mu <= theta%step%E0) then
          ComputeLogPrior = -1.0d300
          return
       end if
    end do
    ! Otherwise, use a flat prior (i.e., add 0)
    ComputeLogPrior = 0.0d0
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

!=====================================================================
! Module: MCMC_Update
! Purpose: Perform blockwise Metropolis-Hastings updates.
!=====================================================================
module MCMC_Update
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use LikelihoodAndPrior
  use RNG
  implicit none
contains

  !----------------------------------------------------------
  ! ProposeNew: Generate a new candidate for a block by adding a small random perturbation.
  ! Uses the intrinsic random_number (or external RNG if needed).
  !----------------------------------------------------------
  subroutine ProposeNew(current_block, proposed_block)
    implicit none
    real(fp_kind), intent(in) :: current_block(:)
    real(fp_kind), intent(out) :: proposed_block(:)
    integer(int32) :: i, n
    real(fp_kind) :: perturb, sigma
    n = size(current_block)
    sigma = 0.1d0   ! Tuning parameter; adjust as needed.
    do i = 1, n
       call RandomUniform(perturb)
       ! Map perturbation from [0,1] to [-sigma, sigma]
       proposed_block(i) = current_block(i) + sigma*(2.0d0*perturb - 1.0d0)
    end do
  end subroutine ProposeNew

  !----------------------------------------------------------
  ! BlockwiseMHUpdate: Perform a Metropolis-Hastings update on a block of parameters.
  ! For demonstration, we assume the block corresponds to a 1D array extracted from theta.
  ! In practice, you must extract and update the relevant parts of the ModelParameters.
  !----------------------------------------------------------
  subroutine BlockwiseMHUpdate(theta, block_indices, beta_val)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    integer(int32), intent(in) :: block_indices(:)  ! Placeholder: indices for block in flattened parameter vector.
    real(fp_kind), intent(in) :: beta_val
    real(fp_kind), allocatable :: current_block(:), proposed_block(:)
    real(fp_kind) :: logPost_current, logPost_proposed, delta, u
    integer(int32) :: i, n

    ! Placeholder: assume each block has one parameter.
    n = size(block_indices)
    allocate(current_block(n), proposed_block(n))
    do i = 1, n
       ! For demonstration, assume we extract a dummy value.
       current_block(i) = 0.0d0  ! Replace with proper extraction from theta.
    end do

    call ProposeNew(current_block, proposed_block)

    logPost_current = ComputeLogPosterior(theta, beta_val)
    ! Here, form a theta_proposed by replacing the dummy block.
    ! In this placeholder, we do not change theta; in a full implementation, update theta accordingly.
    logPost_proposed = ComputeLogPosterior(theta, beta_val)

    delta = logPost_proposed - logPost_current
    call RandomUniform(u)
    if (u < exp(delta)) then
       ! Accept proposal: update the block in theta (placeholder).
       ! In practice, assign proposed_block into theta at indices block_indices.
    else
       ! Reject proposal: do nothing.
    end if

    deallocate(current_block, proposed_block)
  end subroutine BlockwiseMHUpdate

  !----------------------------------------------------------
  ! MCMC_UpdateReplica: Update all blocks of theta for a single replica.
  !----------------------------------------------------------
  subroutine MCMC_UpdateReplica(theta, beta_val)
    implicit none
    type(ModelParameters), intent(inout) :: theta
    real(fp_kind), intent(in) :: beta_val
    integer(int32) :: b
    integer(int32), allocatable :: block_indices(:)
    do b = 1, m
       ! Placeholder: determine block_indices for block b.
       allocate(block_indices(1))
       block_indices(1) = b   ! This is a dummy index.
       call BlockwiseMHUpdate(theta, block_indices, beta_val)
       deallocate(block_indices)
    end do
  end subroutine MCMC_UpdateReplica

end module MCMC_Update

!=====================================================================
! Module: ReplicaExchange
! Purpose: Implement odd–even replica exchange between neighboring replicas.
!=====================================================================
module ReplicaExchange
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use LikelihoodAndPrior
  use GlobalData
  use RNG
  implicit none
contains

  !----------------------------------------------------------
  ! SwapModelParameters: Swap two ModelParameters structures.
  !----------------------------------------------------------
  subroutine SwapModelParameters(theta1, theta2)
    implicit none
    type(ModelParameters), intent(inout) :: theta1, theta2
    type(ModelParameters) :: temp
    temp = theta1
    theta1 = theta2
    theta2 = temp
  end subroutine SwapModelParameters

  !----------------------------------------------------------
  ! OddEvenExchange: Perform replica exchanges between neighboring replicas.
  ! Alternates odd-indexed and even-indexed exchanges based on exchange_step.
  !----------------------------------------------------------
  subroutine OddEvenExchange(theta_array, beta_array, exchange_step)
    implicit none
    type(ModelParameters), dimension(:), intent(inout) :: theta_array
    real(fp_kind), dimension(:), intent(in) :: beta_array
    integer(int32), intent(in) :: exchange_step  ! Exchange step counter
    real(fp_kind) :: logL_lower, logL_higher, delta_swap, r_swap, u
    integer(int32) :: l, num_replica

    num_replica = size(theta_array)

    if (mod(exchange_step, 2) == 1) then
       ! Odd-indexed exchange: swap replicas 1-2, 3-4, ...
       do l = 1, num_replica - 1, 2
          logL_lower = ComputeLogLikelihood(theta_array(l))
          logL_higher = ComputeLogLikelihood(theta_array(l+1))
          delta_swap = (beta_array(l+1) - beta_array(l)) * (logL_lower - logL_higher)
          r_swap = exp(delta_swap)
          call RandomUniform(u)
          if (u < min(1.0d0, r_swap)) then
             call SwapModelParameters(theta_array(l), theta_array(l+1))
          end if
       end do
    else
       ! Even-indexed exchange: swap replicas 2-3, 4-5, ...
       do l = 2, num_replica - 1, 2
          logL_lower = ComputeLogLikelihood(theta_array(l))
          logL_higher = ComputeLogLikelihood(theta_array(l+1))
          delta_swap = (beta_array(l+1) - beta_array(l)) * (logL_lower - logL_higher)
          r_swap = exp(delta_swap)
          call RandomUniform(u)
          if (u < min(1.0d0, r_swap)) then
             call SwapModelParameters(theta_array(l), theta_array(l+1))
          end if
       end do
    end if

  end subroutine OddEvenExchange

end module ReplicaExchange

!=====================================================================
! Module: ThermodynamicIntegration
! Purpose: Accumulate log-likelihoods and compute the marginal likelihood.
!=====================================================================
module ThermodynamicIntegration
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  real(fp_kind), dimension(L_rep) :: SumLogL = 0.0d0  ! Accumulator for each replica
  integer(int32), dimension(L_rep) :: Count = 0           ! Count for each replica
contains

  !----------------------------------------------------------
  ! UpdateLikelihoodAccumulator: Update the log-likelihood accumulator for a given replica.
  !----------------------------------------------------------
  subroutine UpdateLikelihoodAccumulator(replica_index, currentLogL)
    implicit none
    integer(int32), intent(in) :: replica_index
    real(fp_kind), intent(in) :: currentLogL
    SumLogL(replica_index) = SumLogL(replica_index) + currentLogL
    Count(replica_index) = Count(replica_index) + 1
  end subroutine UpdateLikelihoodAccumulator

  !----------------------------------------------------------
  ! ComputeAverageLogLikelihood: Compute average log-likelihood for each replica.
  !----------------------------------------------------------
  subroutine ComputeAverageLogLikelihood(avgLogL)
    implicit none
    real(fp_kind), dimension(:), intent(out) :: avgLogL
    integer(int32) :: l
    if (size(avgLogL) /= L) then
       print *, "Error: avgLogL size mismatch"
       stop
    end if
    do l = 1, L
       if (Count(l) > 0) then
          avgLogL(l) = SumLogL(l) / real(Count(l), real64)
       else
          avgLogL(l) = -1.0d300
       end if
    end do
  end subroutine ComputeAverageLogLikelihood

  !----------------------------------------------------------
  ! ComputeMarginalLikelihood: Use the trapezoidal rule over beta.
  !----------------------------------------------------------
  real(fp_kind) function ComputeMarginalLikelihood(avgLogL, beta_array)
    implicit none
    real(fp_kind), dimension(:), intent(in) :: avgLogL
    real(fp_kind), dimension(:), intent(in) :: beta_array
    integer(int32) :: l, nBeta
    nBeta = size(beta_array)
    ComputeMarginalLikelihood = 0.0d0
    do l = 1, nBeta - 1
       ComputeMarginalLikelihood = ComputeMarginalLikelihood + &
            (beta_array(l+1) - beta_array(l)) * 0.5d0 * (avgLogL(l) + avgLogL(l+1))
    end do
  end function ComputeMarginalLikelihood

  !----------------------------------------------------------
  ! ComputeBayesFreeEnergy: Compute F = -logZ.
  !----------------------------------------------------------
  real(fp_kind) function ComputeBayesFreeEnergy(logZ)
    implicit none
    real(fp_kind), intent(in) :: logZ
    ComputeBayesFreeEnergy = - logZ
  end function ComputeBayesFreeEnergy

end module ThermodynamicIntegration

!=====================================================================
! Module: RNG
! Purpose: Random number generation wrapper.
!=====================================================================
module RNG
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  implicit none
contains
  !----------------------------------------------------------
  ! RandomUniform: Returns a random number uniformly distributed in [0,1].
  !----------------------------------------------------------
  subroutine RandomUniform(u)
    implicit none
    real(fp_kind), intent(out) :: u
    call random_number(u)
  end subroutine RandomUniform
end module RNG

!=====================================================================
! Main Program: BayesianDeconvolution
! Purpose: Integrate all modules to perform Bayesian deconvolution with replica exchange.
!=====================================================================
program BayesianDeconvolution
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use ModelFunctions
  use LikelihoodAndPrior
  use MCMC_Update
  use ReplicaExchange
  use ThermodynamicIntegration
  use RNG
  implicit none
  integer(int32) :: t, l, exchange_step, i
  type(ModelParameters), dimension(L_rep) :: theta_array
  real(fp_kind) :: currentLogL
  real(fp_kind), dimension(L_rep) :: avgLogL

  !--------------------------
  ! Initialize Input Data
  !--------------------------
  do i = 1, N
     E(i) = 1.0d0 * i                  ! Dummy energy values
     I_inc(i) = 100.0d0                ! Dummy incident counts
     I_ab(i) = 50.0d0                  ! Dummy absorbed counts
  end do

  ! Initialize beta (inverse temperatures) for replicas.
  beta(1) = 0.1d0
  beta(2) = 0.4d0
  beta(3) = 0.7d0
  beta(4) = 1.0d0

  !--------------------------
  ! Initialize Model Parameters for each replica.
  ! For demonstration, we use the same initial guess for each replica.
  !--------------------------
  do l = 1, L
     theta_array(l)%step%Ba = 0.0d0
     theta_array(l)%step%Bb = 0.0d0
     theta_array(l)%step%H = 1.0d0
     theta_array(l)%step%E0 = 50.0d0
     theta_array(l)%step%Gamma = 5.0d0

     theta_array(l)%WL%A_WL = 1.0d0
     theta_array(l)%WL%mu_WL = 52.0d0
     theta_array(l)%WL%sigma_G_WL = 2.0d0
     theta_array(l)%WL%gamma_L_WL = 1.0d0

     do i = 1, K1
         theta_array(l)%low(i)%A = 0.5d0
         theta_array(l)%low(i)%mu = 40.0d0 - i    ! Ensure low peaks are below E0
         theta_array(l)%low(i)%sigma_G = 1.0d0
         theta_array(l)%low(i)%gamma_L = 0.5d0
     end do
     do i = 1, K2
         theta_array(l)%high(i)%A = 0.3d0
         theta_array(l)%high(i)%mu = 60.0d0 + i   ! Ensure high peaks are above E0
         theta_array(l)%high(i)%sigma_G = 1.5d0
         theta_array(l)%high(i)%gamma_L = 0.7d0
     end do
  end do

  exchange_step = 0

  !--------------------------
  ! Main MCMC Loop with Replica Exchange
  !--------------------------
  do t = 1, T_iter
     do l = 1, L
        call MCMC_UpdateReplica(theta_array(l), beta(l))
        currentLogL = ComputeLogLikelihood(theta_array(l))
        call UpdateLikelihoodAccumulator(l, currentLogL)
     end do

     exchange_step = exchange_step + 1
     call OddEvenExchange(theta_array, beta, exchange_step)

     ! Optionally: Save samples for replica with beta=1 (typically theta_array(L)) after T_burn.
  end do

  !--------------------------
  ! Post-processing: Compute average log-likelihood and marginal likelihood.
  !--------------------------
  call ComputeAverageLogLikelihood(avgLogL)
  print *, "Average log-likelihood for each replica:"
  do l = 1, L
     print *, "Replica ", l, " beta = ", beta(l), " avgLogL = ", avgLogL(l)
  end do

  print *, "Marginal likelihood (log Z): ", ComputeMarginalLikelihood(avgLogL, beta)
  print *, "Bayes free energy F: ", ComputeBayesFreeEnergy( ComputeMarginalLikelihood(avgLogL, beta) )

  !--------------------------
  ! Output: Posterior samples from replica with beta = 1 (theta_array(L))
  !--------------------------
  print *, "Posterior sample for beta=1 chain:"
  print *, "Step function E0: ", theta_array(L)%step%E0

end program BayesianDeconvolution
