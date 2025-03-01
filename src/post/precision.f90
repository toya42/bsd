!=====================================================================
! Module: precision
! Purpose: Define the floating-point kind (e.g. real64) for use throughout.
!=====================================================================
module precision
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  integer(int32), parameter :: fp_kind = real64  ! Change to real32 if desired.
end module precision
