program BayesianDeconvolution
  use GlobalData
  use DataInput
  use ModelFunctions
  use LikelihoodAndPrior
  use MCMC_Update
  use ReplicaExchange
  use ThermodynamicIntegration
  use RNG
  implicit none
  integer :: t, l, exchange_step, i
  type(ModelParameters), allocatable, dimension(:) :: theta_array
  real(8) :: currentLogL
  real(8), allocatable, dimension(:) :: avgLogL

  !-----------------------------------------------------------
  ! Step 1: Read experimental data from CSV file.
  !-----------------------------------------------------------
  call ReadExperimentalData()

  !-----------------------------------------------------------
  ! Step 2: Read MCMC parameters from parameter file.
  !-----------------------------------------------------------
  call ReadMCMCParameters()

end program BayesianDeconvolution