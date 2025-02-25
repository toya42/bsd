!=====================================================================
! Main Program: BayesianDeconvolution
! Purpose: Integrate all modules to perform Bayesian deconvolution
!          with replica exchange.
!=====================================================================
program BayesianDeconvolution
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only : fp_kind
  use GlobalData
  use DataInput
  use DataOutput
  use ModelFunctions
  use LikelihoodAndPrior
  use MCMC_Update
  use ReplicaExchange
  use ThermodynamicIntegration
  use PosteriorStatistics
  use RNG
  implicit none

  integer(int32) :: t, l, exchange_step, i
  type(ModelParameters), allocatable, dimension(:) :: theta_array
  type(ModelParameters) :: theta_mode
  real(fp_kind) :: currentLogL
  real(fp_kind), allocatable, dimension(:) :: avgLogL

  !-----------------------------------------------------------
  ! Read experimental data from CSV file.
  !-----------------------------------------------------------
  call ReadExperimentalData()

  !-----------------------------------------------------------
  ! Read MCMC parameters from parameter file.
  !-----------------------------------------------------------
  call ReadMCMCParameters()


  allocate(theta_array(L_rep),avgLogL(L_rep))
  do l=1,L_rep
    allocate(theta_array(l)%low(K1))
    allocate(theta_array(l)%high(K2))
  end do


  call InitializeHistoryOutput("log")

  !-----------------------------------------------------------
  ! Initialize Beta (Inverse Temperatures) for Replicas
  ! In practice, beta might be read from a parameter file.
  !-----------------------------------------------------------
  !beta(1) = 0.1d0
  !beta(2) = 0.4d0
  !beta(3) = 0.7d0
  !beta(4) = 1.0d0
  beta(1) = 1e-2
  do l=2,L_rep
    beta(l) = 1.225**(l-L_rep)
    !print *,beta(l)
   end do
  !-----------------------------------------------------------
  ! Initialize Model Parameters for Each Replica.
  ! For demonstration, we use the same initial guess for each replica.
  !-----------------------------------------------------------
  do l = 1, L_rep
     theta_array(l)%step%Ba = 0.0d0
     theta_array(l)%step%Bb = 5.0d-8
     theta_array(l)%step%H  = 5.0d-6
     theta_array(l)%step%E0 = 2480.0d0
     theta_array(l)%step%Gamma = 1.0d-1

     theta_array(l)%WL%A_WL = 2.0d-5
     theta_array(l)%WL%mu_WL = 0.0d0
     theta_array(l)%WL%sigma_G_WL = 2.0d-1
     theta_array(l)%WL%gamma_L_WL = 2.0d-1

     do i = 1, K1
         theta_array(l)%low(i)%A = 0.0d0
         theta_array(l)%low(i)%mu = real((K1+1-i)*5)    ! Ensure low-energy peaks are below E0
         theta_array(l)%low(i)%sigma_G = 2.0d-1
         theta_array(l)%low(i)%gamma_L = 2.0d-1
     end do
     do i = 1, K2
         theta_array(l)%high(i)%A = 0.0d0
         theta_array(l)%high(i)%mu = real(i*10)   ! Ensure high-energy peaks are above E0
         theta_array(l)%high(i)%sigma_G = 2.0d-1
         theta_array(l)%high(i)%gamma_L = 2.0d-1
     end do
  end do

  !call InitializeSpectrumOutput("spec0")
  !print *,3
  !call ExportRestoredSpectrum(theta_array(L_rep))
  !stop

  !-----------------------------------------------------------
  ! Set the replica exchange counter to zero.
  !-----------------------------------------------------------
  exchange_step = 0

  !-----------------------------------------------------------
  ! Main MCMC Loop with Replica Exchange.
  !-----------------------------------------------------------
  do t = 1, T_iter
     if(mod(t,500)==0) then
      print *,'iteration:',t
     end if  
     do l = 1, L_rep
        ! Update the parameters for replica l using blockwise MH updates.
        call MCMC_UpdateReplica(theta_array(l), beta(l))
        call MCMC_block_UpdateReplica(theta_array(l), beta(l))
        ! Compute the current log-likelihood for replica l.
        currentLogL = ComputeLogLikelihood(theta_array(l))
        ! Update the likelihood accumulator for replica l.
        call UpdateLikelihoodAccumulator(l, currentLogL)
     end do

     ! Increment the exchange step counter.
    if(mod(t,50)==0) then
        exchange_step = exchange_step + 1
        ! Perform odd–even replica exchange across the replicas.
        call OddEvenExchange(theta_array, beta, exchange_step)
    end if

     if(t > T_burn) then 
          call AppendHistoryEntry(t, theta_array(L_rep))
     end if
     !call AppendHistoryEntry(t, theta_array(L_rep))
  end do

  call FinalizeHistoryOutput()

  !-----------------------------------------------------------
  ! Post-processing: Compute Average Log-Likelihoods and Marginal Likelihood.
  !-----------------------------------------------------------
  call ComputeAverageLogLikelihood(avgLogL)
  print *, "Average log-likelihood for each replica:"
  do l = 1, L_rep
     print *, "Replica ", l, " beta = ", beta(l), " avgLogL = ", avgLogL(l)
  end do

  print *, "Marginal likelihood (log Z): ", ComputeMarginalLikelihood(avgLogL, beta)
  print *, "Bayes free energy F: ", ComputeBayesFreeEnergy( ComputeMarginalLikelihood(avgLogL, beta) )

  !-----------------------------------------------------------
  ! Output: Posterior Samples from the Replica with beta = 1.
  !-----------------------------------------------------------
  print *, "Posterior sample for beta=1 chain:"
  print *, "Step function E0: ", theta_array(L_rep)%step%E0

  !print *,1
  call ComputePosteriorMode("log", theta_mode)
  !print *,2
  call InitializeSpectrumOutput("spec")
  !print *,3
  call ExportRestoredSpectrum(theta_mode)
  !print *,4
  call FinalizeSpectrumOutput()
  !print *,5

end program BayesianDeconvolution
