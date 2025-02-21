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
