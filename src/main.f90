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
  use PriorProposal, only : c_proposal
  use FaddeevaTable
  use SortPeaks
  implicit none

  integer(int32) :: t, l, exchange_step, i, b, exchange_unit, idx_unit, now_unit
  integer(int32),allocatable,dimension(:) :: idx_now
  type(ModelParameters), allocatable, dimension(:) :: theta_array
  type(ModelParameters) :: theta_mode
  real(fp_kind) :: currentLogL,accept_ratio
  real(fp_kind), allocatable, dimension(:) :: avgLogL
  character(len=20) :: fmt_exchange, fmt_idx, fmt_now

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


  call InitializeHistoryOutput("log.txt")
  call InitializeFaddeevaTable
  !-----------------------------------------------------------
  ! Initialize Beta (Inverse Temperatures) for Replicas
  ! In practice, beta might be read from a parameter file.
  !-----------------------------------------------------------
  !beta(1) = 0.1d0
  !beta(2) = 0.4d0
  !beta(3) = 0.7d0
  !beta(4) = 1.0d0
  do l=1,L_rep
    !beta(l) = 1.2**(l-L_rep)
    !beta(l) = log(l+0.1d0)/log(L_rep+0.1)
    !beta(l) = real(l)/real(L_rep)
    !beta(l) = (real(l)/real(L_rep)+1.2d0**(l-L_rep))*0.5d0
    beta(l) = 0.5d0*tanh((l-L_rep*0.5)/(L_rep*0.2))+0.5
    idx_exchange(l) = l
    !print *,beta(l)
  end do
  beta(L_rep) = 1.0d0
  !beta(1) = beta(2)*0.5
  !-----------------------------------------------------------
  ! Initialize Model Parameters for Each Replica.
  ! For demonstration, we use the same initial guess for each replica.
  !-----------------------------------------------------------
  do l = 1, L_rep
     theta_array(l)%step%Ba = 0.0d0
     theta_array(l)%step%Bb = 0.0d0
     theta_array(l)%step%H  = 5.0d-6
     theta_array(l)%step%E0 = 2480.0d0
     theta_array(l)%step%Gamma = 1.0d0

     theta_array(l)%WL%A_WL = 2.0d-6
     theta_array(l)%WL%mu_WL = 0.0d0
     theta_array(l)%WL%sigma_G_WL = 2.0d-1
     theta_array(l)%WL%gamma_L_WL = 2.0d-1

     do i = 1, K1
         theta_array(l)%low(i)%A = 2.0d-6
         theta_array(l)%low(i)%mu = real((K1+1-i))    ! Ensure low-energy peaks are below E0
         theta_array(l)%low(i)%sigma_G = 2.0d-1
         theta_array(l)%low(i)%gamma_L = 2.0d-1
     end do
     do i = 1, K2
         theta_array(l)%high(i)%A = 2.0d-6
         theta_array(l)%high(i)%mu = real(i*10)   ! Ensure high-energy peaks are above E0
         theta_array(l)%high(i)%sigma_G = 2.0d-1
         theta_array(l)%high(i)%gamma_L = 2.0d-1
     end do
  end do

  !call InitializeSpectrumOutput("spec0")
  !print *,3
  !call ExportRestoredSpectrum(theta_array(L_rep))
  !stop

  call InitializeCounters(0)
  allocate(c_proposal((2+K1+K2),L_rep))
  c_proposal = 0.5d-2

  cnt_exchange = 0.0d0
  write(fmt_exchange,'(I0)') L_rep-1
  fmt_exchange = '(i8,'//trim(fmt_exchange)//'(f8.2))'
  !print *,fmt_exchange
  exchange_unit=19
  open(exchange_unit,file="exchange.txt",status='replace', action='write', form='formatted')

  write(fmt_idx,'(I0)') L_rep
  fmt_idx = '(i8,'//trim(fmt_idx)//'(i8))'
  !print *,fmt_exchange
  idx_unit=20
  open(idx_unit,file="index.txt",status='replace', action='write', form='formatted')

  write(fmt_now,'(I0)') L_rep
  fmt_now = '(i8,'//trim(fmt_now)//'(i8))'
  !print *,fmt_exchange
  now_unit=21
  open(now_unit,file="now.txt",status='replace', action='write', form='formatted')

  allocate(idx_now(L_rep))

  !-----------------------------------------------------------
  ! Set the replica exchange counter to zero.
  !-----------------------------------------------------------
  exchange_step = 0


  !-----------------------------------------------------------
  ! Main MCMC Loop with Replica Exchange.
  !-----------------------------------------------------------
  do t = 1, T_iter
    !print *,t

    if(mod(t,100)==0) then
      print *,'iteration:',t
    end if

    do l = 1, L_rep
        ! Update the parameters for replica l using blockwise MH updates.
        !call MCMC_UpdateReplica(theta_array(l), beta(l),l)
        call MCMC_block_UpdateReplica(theta_array(l), beta(l),l)
        ! Compute the current log-likelihood for replica l.
        currentLogL = ComputeLogLikelihood(theta_array(l))
        ! Update the likelihood accumulator for replica l.
        call UpdateLikelihoodAccumulator(l, currentLogL)
    end do


    ! Increment the exchange step counter.
    if(mod(t,100)==0) then
      exchange_step = exchange_step + 1
      ! Perform odd–even replica exchange across the replicas.
      call OddEvenExchange(theta_array, beta, exchange_step)
      write(idx_unit,fmt_idx) t,idx_exchange(1:L_rep)
      do l=1,L_rep
        do b=1,L_rep
          if(l==idx_exchange(b)) then
            idx_now(l) = b
          end if
        end do
      end do
      write(now_unit,fmt_now) t,idx_now(1:L_rep)
    end if

    ! tune c_proposal
    if(mod(t,100)==0 .and. t<=T_burn/2) then
      do l=1,L_rep
        do b=1,(2+K1+K2)
          accept_ratio = real(accepted_proposals(b,l))/real(total_proposals(b,l))*100
          if(accept_ratio<20.0) then
            c_proposal(b,l) = c_proposal(b,l)*0.9
          else if(accept_ratio>50.0) then
            c_proposal(b,l) = c_proposal(b,l)*1.1
          end if
          !if(t==T_burn/2) then
          !  print *,'iteration:',t
          !  print *,'Replica:',l
          !  print *,'block:',b
          !  print '("accept ratio(%) = ",f6.2)',accept_ratio
          !  print '(i5,"/",i5)', accepted_proposals(b,l),total_proposals(b,l)
          !  print '("c_proposal = ",f9.5)',c_proposal(b,l)
          !end if
        end do
      end do
      call InitializeCounters(1)
      if(t==T_burn/2) then
        print *,'c_proposal tuning step completed'
      end if
    end if

    if(t==T_burn) then
      print *,'burn-in step completed'
    end if

    ! output exchange history
    if(mod(t,2000)==0) then
      write(exchange_unit,fmt_exchange) t,real(cnt_exchange(1:L_rep-1))/real(total_exchange(1:L_rep-1))*1.0d2
      !print *, 'exchange output'
      !print fmt_exchange,t, real(cnt_exchange(1:L_rep))/real(50)*1.0d2
      block
        real(fp_kind) :: exr1, exr2
        real(fp_kind),parameter :: c_delta = 0.2d0
        if(t<=T_burn/2) then
          do l=1,L_rep-2
            exr1 = real(cnt_exchange(l  ))/real(total_exchange(l  ))*1.0d2
            exr2 = real(cnt_exchange(l+1))/real(total_exchange(l+1))*1.0d2
            if(exr1>=60.0 .and. exr2>=60.0) then
              !print *, 'tune beta (case1)'
              !print *, 'l=',l+1,l+2
              !print *, '<beta_before>',beta(l),beta(l+1)
              beta(l+1) = beta(l+1)-(beta(l+1)-beta(l))*c_delta
              beta(l+2) = beta(l+2)-(beta(l+2)-beta(l+1))*c_delta
              !print *, '<beta_after >',beta(l+1),beta(l+2)
            else if(exr1<exr2 .and. exr2>=60.0) then
              !print *, 'tune beta (case2)'
              !print *, 'l=',l+1
              !print *, '<beta_before>',beta(l+1)
              beta(l+1) = beta(l+1)-(beta(l+1)-beta(l))*c_delta
              !print *, '<beta_after >',beta(l+1)
            else if(exr1>exr2 .and. exr1>=60.0) then
              !print *, 'tune beta (case3)'
              !print *, 'l=',l+1
              !print *, '<beta_before>',beta(l+1)
              beta(l+1) = beta(l+1)+(beta(l+2)-beta(l+1))*c_delta
              !print *, '<beta_after >',beta(l+1)
            end if
          end do
          beta(L_rep) = 1.0d0
          call BubbleSortBetas(beta, .true.)
          if(t==T_burn/2) then
            do l=1,L_rep
              write(31,*) l,beta(l)
            end do
          end if
        end if
      end block
      cnt_exchange = 0
      total_exchange = 0
    end if

    ! store and output data
    if(t > T_burn) then 
      call AppendHistoryEntry(t, theta_array(L_rep))
    end if
     !call AppendHistoryEntry(t, theta_array(L_rep))
    

  end do

  print *,'mcmc step completed'

  close(idx_unit)
  close(exchange_unit)
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
  call ComputePosteriorMode("log.txt", theta_mode)
  !print *,2
  call InitializeSpectrumOutput("spec.txt")
  !print *,3
  call ExportRestoredSpectrum(theta_mode)
  !print *,4
  call FinalizeSpectrumOutput()
  !print *,5

end program BayesianDeconvolution
