!=====================================================================
! Module: ModelFunctions
! Purpose: Define model functions for the absorption spectrum.
!          Implements a Voigt profile function using the Faddeeva method.
!=====================================================================
module ModelFunctions
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use precision, only: fp_kind
  use GlobalData
  use FaddeevaTable
  implicit none

  ! Derived type for the step (edge) function parameters.
  type :: ParametersStep
     real(fp_kind) :: Ba, Bb  ! Baseline: modeled as a linear function: Ba*E + Bb.
     real(fp_kind) :: H       ! Step height.
     real(fp_kind) :: E0      ! Edge energy (absorption-edge position).
     real(fp_kind) :: Gamma   ! Broadening parameter.
  end type ParametersStep

  ! Derived type for the white-line (peak) parameters.
  type :: ParametersWL
     real(fp_kind) :: A_WL       ! Peak amplitude.
     real(fp_kind) :: mu_WL      ! Peak center (near E0).
     real(fp_kind) :: sigma_G_WL ! Gaussian width.
     real(fp_kind) :: gamma_L_WL ! Lorentzian half-width.
  end type ParametersWL

  ! Derived type for individual Voigt peak parameters (for low and high energy domains).
  type :: ParametersPeak
     real(fp_kind) :: A       ! Peak amplitude.
     real(fp_kind) :: mu      ! Peak center energy.
     real(fp_kind) :: sigma_G ! Gaussian width.
     real(fp_kind) :: gamma_L ! Lorentzian half-width.
  end type ParametersPeak

  ! Combined type for all model parameters.
  type :: ModelParameters
     type(ParametersStep) :: step
     type(ParametersWL)   :: WL
     type(ParametersPeak), allocatable, dimension(:) :: low   ! Low-energy peaks.
     type(ParametersPeak), allocatable, dimension(:) :: high  ! High-energy peaks.
  end type ModelParameters

contains

  !---------------------------------------------------------------------
  ! Function: Voigt
  ! Purpose: Compute the Voigt profile V(x; sigma, gamma)
  !          using the Faddeeva function.
  ! Formula: V(x; sigma, gamma) = Re[w(z)] / (sigma * sqrt(2*pi)),
  ! where z = (x + i*gamma) / (sigma * sqrt(2)).
  !---------------------------------------------------------------------
  real(fp_kind) function Voigt(x, sigma, gamma)
    implicit none
    real(fp_kind), intent(in) :: x, sigma, gamma
    real(fp_kind) :: scale
    real(fp_kind) :: xx, yy
    complex(fp_kind) :: w_val

    !if(sigma<0.0d0 .or. gamma<0.0d0) then
    !  print *,'sigma',sigma
    !  print *,'gamma',gamma
    !end if

    ! Compute scaling factor: sigma * sqrt(2)
    scale = sigma * sqrt(2.0d0)
    ! Compute real and imaginary parts of z.
    xx = x / scale
    yy = gamma / scale
    ! Compute Faddeeva function w(z).
    !print *,xx,yy
    !w_val = Faddeeva(xx, yy)
    w_val = FastFaddeeva(abs(xx), yy)
    ! Voigt profile is the real part divided by (sigma * sqrt(2*pi)).
    Voigt = real(w_val) / ( sigma * sqrt(2.0d0 * pi) )
    !if(Voigt<0.0) then
    !    print *,"Voigt=",Voigt
    !end if
  end function Voigt

  !---------------------------------------------------------------------
  ! f_step: Smoothed step function modeling the absorption edge.
  ! f_step =  Bb + H * Φ((E - E0)/(0.5*Gamma)), (Ba:not used)
  ! where Φ(x) = 0.5 + (1/π) atan(x)
  !---------------------------------------------------------------------
  real(fp_kind) function f_step(E, theta_step)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ParametersStep), intent(in) :: theta_step
    real(fp_kind) :: x
    ! Scale the argument by 2*Gamma as specified.
    x = (E - theta_step%E0) / (0.5 * theta_step%Gamma)
    f_step = theta_step%Ba * (E - theta_step%E0) + theta_step%Bb + theta_step%H * (0.5d0 + (1.0d0/pi)*atan(x))
    !f_step = theta_step%Ba * E + theta_step%Bb + theta_step%H * (0.5d0 + (1.0d0/pi)*atan(x))
    !f_step = theta_step%Bb + theta_step%H * (0.5d0 + (1.0d0/pi)*atan(x))
  end function f_step

  !---------------------------------------------------------------------
  ! f_WL: White-line peak using a Voigt profile.
  !---------------------------------------------------------------------
  real(fp_kind) function f_WL(E, theta_WL, E0)
    implicit none
    real(fp_kind), intent(in) :: E, E0
    type(ParametersWL), intent(in) :: theta_WL
    real(fp_kind) :: mu
    mu = E0+theta_WL%mu_WL
    !f_WL = theta_WL%A_WL * Voigt(E - theta_WL%mu_WL, theta_WL%sigma_G_WL, theta_WL%gamma_L_WL)
    f_WL = theta_WL%A_WL * Voigt(E - mu, theta_WL%sigma_G_WL, theta_WL%gamma_L_WL)

  end function f_WL

  !---------------------------------------------------------------------
  ! f_low: Sum of Voigt profiles for low-energy peaks.
  !---------------------------------------------------------------------
  real(fp_kind) function f_low(E, theta_low, E0)
    implicit none
    real(fp_kind), intent(in) :: E, E0
    type(ParametersPeak), dimension(:), intent(in) :: theta_low
    integer(int32) :: k
    real(fp_kind) :: mu

    f_low = 0.0d0
    do k = 1, size(theta_low)
       mu = E0-theta_low(k)%mu
       f_low = f_low + theta_low(k)%A * Voigt(E - mu, theta_low(k)%sigma_G, theta_low(k)%gamma_L)
    end do
  end function f_low

  !---------------------------------------------------------------------
  ! f_high: Sum of Voigt profiles for high-energy peaks.
  !---------------------------------------------------------------------
  real(fp_kind) function f_high(E, theta_high, E0)
    implicit none
    real(fp_kind), intent(in) :: E, E0
    type(ParametersPeak), dimension(:), intent(in) :: theta_high
    integer(int32) :: j
    real(fp_kind) :: mu

    f_high = 0.0d0
    do j = 1, size(theta_high)
       mu = E0+theta_high(j)%mu
       f_high = f_high + theta_high(j)%A * Voigt(E - mu, theta_high(j)%sigma_G, theta_high(j)%gamma_L)
    end do
  end function f_high

  !---------------------------------------------------------------------
  ! f_ratio: Full absorption ratio function as the sum of components.
  !---------------------------------------------------------------------
  real(fp_kind) function f_ratio(E, theta)
    implicit none
    real(fp_kind), intent(in) :: E
    type(ModelParameters), intent(in) :: theta

    !print *,"f_ratio",f_ratio
    !print *, theta%step
    !print *,"step",f_step(E, theta%step)
    !print *,"WL",f_WL(E, theta%WL)
    !print *,"Low",f_low(E, theta%low)
    !print *,"High",f_high(E, theta%high)


    f_ratio = f_step(E, theta%step) + f_WL(E, theta%WL, theta%step%E0) &
                                   + f_low(E, theta%low, theta%step%E0)&
                                   + f_high(E, theta%high, theta%step%E0)

    !if(f_ratio<0.0) then
    !    print *,"negative f_ratio"
    !    print *,"f_ratio",f_ratio
    !    print *,"step",f_step(E, theta%step)
    !    print *,"WL",f_WL(E, theta%WL)
    !    print *,"Low",f_low(E, theta%low)
    !    print *,"High",f_high(E, theta%high)
    !end if

  end function f_ratio

end module ModelFunctions
