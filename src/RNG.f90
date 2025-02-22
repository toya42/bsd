!=====================================================================
! Module: RNG
! Purpose: Random number generation wrapper.
!=====================================================================
module RNG
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData, only: pi
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

  !----------------------------------------------------------
  ! Subroutine: NormalRandom
  ! Purpose: Generate a normally distributed random number (mean=0, variance=1)
  !          using the Box-Muller transform.
  !----------------------------------------------------------
  subroutine NormalRandom(x)
    implicit none
    real(fp_kind), intent(out) :: x
    real(fp_kind) :: u1, u2, r, theta
    call RandomUniform(u1)
    call RandomUniform(u2)
    ! Avoid taking log(0)
    if (u1 == 0.0_fp_kind) then
       u1 = 1.0e-10_fp_kind
    end if
    r = sqrt(-2.0_fp_kind * log(u1))
    theta = 2.0_fp_kind * pi * u2
    x = r * cos(theta)
  end subroutine NormalRandom
end module RNG
