!=====================================================================
! Module: GlobalData
! Purpose: Store global input arrays and configuration parameters.
! Note: Variables that used to be hard-coded parameters (N, T_iter, etc.)
!       are now read from input files.
!=====================================================================
module GlobalData
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  implicit none

  ! mathematical constant
  real(fp_kind), parameter :: pi=4.0_real64*atan(1.0_real64)

  integer(int32), parameter :: L_rep_max = 20
  ! These will be set by reading the experimental data file and the MCMC parameter file.
  integer(int32) :: N = 0          ! Number of measurement points (to be read from CSV)
  integer(int32) :: T_iter = 0     ! Total number of MCMC iterations
  integer(int32) :: T_burn = 0     ! Burn-in iterations
  integer(int32) :: L_rep = 0      ! Number of replicas (for replica exchange)
  integer(int32) :: K1 = 0         ! Number of low-energy peaks
  integer(int32) :: K2 = 0         ! Number of high-energy peaks

  ! Input arrays: these will be allocated after reading the CSV file.
  real(fp_kind), allocatable, dimension(:) :: E      ! Energy array
  real(fp_kind), allocatable, dimension(:) :: I_inc  ! Incident counts (scaled)
  real(fp_kind), allocatable, dimension(:) :: I_ab   ! Absorbed counts (scaled)

  ! Replica exchange settings: inverse temperatures (beta)
  real(fp_kind), allocatable, dimension(:) :: beta  ! Size will be L_rep

end module GlobalData
