module MCMCStats
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use ModelFunctions
  use PriorProposal, only : c_proposal
  implicit none
  integer(int32), parameter :: length_stats_cal=1000
  integer(int32),private, parameter :: stats_tau_unit = 63
  integer(int32),private, parameter :: stats_ess_unit = 64
  integer(int32),private, parameter :: c_proposal_unit = 65
  real(fp_kind), parameter :: rate_cp = 0.005d0
  real(fp_kind), parameter :: tau_target = 10.0
  integer(int32) :: stats_count = 0
  integer(int32) :: stats_output_count = 0

  real(fp_kind), private, allocatable, dimension(:,:,:) :: history_buffer
  integer(int32) :: parameter_length
  character(len=20), private :: fmt

contains

   !----------------------------------------------------------
   ! Subroutine: InitializeHistoryOutput
   ! Purpose: Open the history output file and initialize the buffer.
   ! Input:
   !   filename: Name of the output CSV file for history.
   !----------------------------------------------------------
  subroutine InitializeHistoryOutputForStats()
    implicit none
    !character(len=*), intent(in) :: filename

    parameter_length = 1 + 5 + 4 + 4*K1 + 4*K2
    allocate(history_buffer(parameter_length, L_rep, length_stats_cal))
    stats_count = 0
    stats_output_count = 0

     write(fmt,'(I0)') parameter_length
     fmt = '('//trim(fmt)//'(E20.8e3))'


    open(unit=stats_tau_unit, file='tau_tune_history.txt',form='formatted', status='replace', action='write')
    open(unit=stats_ess_unit, file='ess_tune_history.txt',form='formatted', status='replace', action='write')
    open(unit=c_proposal_unit, file='c_proposal_tune_history.txt',form='formatted', status='replace', action='write')


  end subroutine InitializeHistoryOutputForStats

   !----------------------------------------------------------
   ! Subroutine: AppendHistoryEntry
   ! Purpose: Append a history entry (current theta parameters) for a given iteration.
   !----------------------------------------------------------
   subroutine AppendHistoryForStats(iteration, l, theta)
     implicit none
     integer, intent(in) :: iteration,l
     type(ModelParameters), intent(in) :: theta
     real(fp_kind), dimension(parameter_length) :: params_one_d
     integer :: pos, i
 
     pos = 1
     params_one_d(pos) = real(iteration, fp_kind)
     pos = pos + 1
 
     ! Step parameters: [Ba, Bb, H, E0, Gamma]
     params_one_d(pos:pos+4) = [ theta%step%Ba, theta%step%Bb, theta%step%H, theta%step%E0, theta%step%Gamma ]
     pos = pos + 5
 
     ! White-line parameters: [A_WL, mu_WL, sigma_G_WL, gamma_L_WL]
     params_one_d(pos:pos+3) = [ theta%WL%A_WL, theta%WL%mu_WL, theta%WL%sigma_G_WL, theta%WL%gamma_L_WL ]
     pos = pos + 4
 
     ! Low-energy peaks: each has 4 parameters.
     do i = 1, K1
        params_one_d(pos:pos+3) = [ theta%low(i)%A, theta%low(i)%mu, theta%low(i)%sigma_G, theta%low(i)%gamma_L ]
        pos = pos + 4
     end do
 
     ! High-energy peaks: each has 4 parameters.
     do i = 1, K2
        params_one_d(pos:pos+3) = [ theta%high(i)%A, theta%high(i)%mu, theta%high(i)%sigma_G, theta%high(i)%gamma_L ]
        pos = pos + 4
     end do
 
     ! Add the entry to the buffer.
     !if(l==L_rep) stats_count = stats_count + 1
     history_buffer(:,l,stats_count) = params_one_d
 
     ! If buffer is full, flush it.
     if (stats_count == length_stats_cal) then
        !print *,'append l',l
        call CalcStats(l)
     end if
   end subroutine AppendHistoryForStats
 
   !----------------------------------------------------------
   ! Subroutine: FlushHistoryBuffer
   ! Purpose: Write the buffered history entries to file and reset the buffer.
   !----------------------------------------------------------
  subroutine CalcStats(l)
    implicit none
    integer(int32), intent(in) :: l
    integer :: i,acf_size
    !real(fp_kind) :: tau_int,ess
    real(fp_kind), allocatable,dimension(:) :: acf, p_temp
    real(fp_kind),dimension(parameter_length-1) :: tau_int, ess

    acf_size = min(length_stats_cal/5, length_stats_cal-1)
    allocate(acf(0:acf_size), p_temp(length_stats_cal))

    !print *,'l',l

    do i=1,parameter_length-1
      p_temp = history_buffer(i+1,l,:)
      acf = ComputeAutocorrelation(p_temp,acf_size)
      tau_int(i) = ComputeIntegratedAutocorrelationTime(acf)
      ess(i) = ComputeEffectiveSampleSize(p_temp,tau_int(i))

      if(tau_int(i)==401.0) then
        c_proposal(i,l) = c_proposal(i,l)*0.5d0
        cycle
      end if

      if(tau_int(i)>tau_target) then
        !print *,'l,i:',l,i
        !print *,'before',c_proposal(i,l)
        c_proposal(i,l) = c_proposal(i,l)*exp(rate_cp*(tau_int(i)-tau_target)/tau_target)
        !print *,'after',c_proposal(i,l)
      end if

    end do

    if(l==L_rep) then
      stats_output_count = stats_output_count+1
      write(stats_tau_unit,fmt) real(stats_output_count*length_stats_cal,fp_kind),tau_int(:)
      write(stats_ess_unit,fmt) real(stats_output_count*length_stats_cal,fp_kind),    ess(:)
      write(c_proposal_unit,fmt) real(stats_output_count*length_stats_cal,fp_kind),c_proposal(:,l)
      stats_count = 0
    end if


 

   end subroutine CalcStats


  !----------------------------------------------------------
  ! Function: ComputeAutocorrelation
  ! Purpose: Compute the autocorrelation function (ACF) of a time series x
  !          up to a maximum lag maxlag.
  !
  ! Input:
  !   x      - Real vector of MCMC samples.
  !   maxlag - Maximum lag to compute the ACF.
  !
  ! Output:
  !   acf    - Real vector of length (maxlag+1) containing ACF values,
  !            with acf(0) = 1.
  !----------------------------------------------------------
  function ComputeAutocorrelation(x, maxlag_in) result(acf)
    implicit none
    real(fp_kind), intent(in) :: x(:)
    integer, intent(in) :: maxlag_in
    real(fp_kind), allocatable, dimension(:) :: acf
    integer :: n, k, t, maxlag
    real(fp_kind) :: mean_x, var_x, numerator

    n = size(x)
    maxlag = min(maxlag_in, n-1)
    allocate(acf(0:maxlag))
    mean_x = sum(x) / n
    var_x = sum((x - mean_x)**2) / n

    acf(0) = 1.0d0
    do k = 1, maxlag
      numerator = 0.0d0
      do t = 1, n - k
        numerator = numerator + (x(t) - mean_x) * (x(t+k) - mean_x)
      end do
      !print *,k
      !print *,n-k 
      !print *,numerator
      !print *, var_x
      acf(k) = numerator / real(n - k, fp_kind) / var_x
    end do
  end function ComputeAutocorrelation

  !----------------------------------------------------------
  ! Function: ComputeIntegratedAutocorrelationTime
  ! Purpose: Compute the integrated autocorrelation time tau_int from an ACF vector.
  !
  ! Input:
  !   acf - Autocorrelation function vector, with acf(0)=1.
  !
  ! Output:
  !   tau_int - Integrated autocorrelation time.
  !----------------------------------------------------------
  function ComputeIntegratedAutocorrelationTime(acf) result(tau_int)
    implicit none
    real(fp_kind), intent(in) :: acf(:)
    real(fp_kind) :: tau_int
    integer :: k, maxlag

    maxlag = size(acf) - 1
    tau_int = 1.0d0
    do k = 1, maxlag
       tau_int = tau_int + 2.0d0 * acf(k)
    end do
  end function ComputeIntegratedAutocorrelationTime

  !----------------------------------------------------------
  ! Function: ComputeEffectiveSampleSize
  ! Purpose: Compute the effective sample size (ESS) for a time series x,
  !          given its integrated autocorrelation time tau_int.
  !
  ! Input:
  !   x       - Time series (vector of MCMC samples).
  !   tau_int - Integrated autocorrelation time.
  !
  ! Output:
  !   ess - Effective sample size, N / tau_int.
  !----------------------------------------------------------
  function ComputeEffectiveSampleSize(x, tau_int) result(ess)
    implicit none
    real(fp_kind), intent(in) :: x(:)
    real(fp_kind), intent(in) :: tau_int
    real(fp_kind) :: ess

    ess = size(x) / tau_int
  end function ComputeEffectiveSampleSize

end module MCMCStats
