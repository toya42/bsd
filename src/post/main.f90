!=====================================================================
! Main Program: BayesianDeconvolution
! Purpose: Integrate all modules to perform Bayesian deconvolution
!          with replica exchange.
!=====================================================================
program PostProcessing
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use DataInput
  use DataOutput
  use ModelFunctions
  use PosteriorStatistics
  use FaddeevaTable
  implicit none

  integer(int32) :: l
  type(ModelParameters), allocatable, dimension(:) :: theta_array
  type(ModelParameters) :: theta_mode

  !-----------------------------------------------------------
  ! Read experimental data from CSV file.
  !-----------------------------------------------------------
  call ReadExperimentalData()

  !-----------------------------------------------------------
  ! Read MCMC parameters from parameter file.
  !-----------------------------------------------------------
  call ReadMCMCParameters()

  print *,'bin:'
  read *, NBINS
  print *,NBINS

  print *,'nhist:'
  read *,nhist
  print *,nhist

  allocate(theta_array(L_rep))
  do l=1,L_rep
    allocate(theta_array(l)%low(K1))
    allocate(theta_array(l)%high(K2))
  end do


  call InitializeHistoryOutput("log.txt")
  call InitializeFaddeevaTable
  !print *,1
  call ComputePosteriorMode("log.txt", theta_mode)
  !print *,2
  call InitializeSpectrumOutput("spec.txt")
  !print *,3
  call ExportRestoredSpectrum(theta_mode)
  !print *,4
  call FinalizeSpectrumOutput()
  !print *,5

end program PostProcessing
