module FaddeevaTable
  use, intrinsic :: iso_fortran_env
  use precision, only: fp_kind
  implicit none
  ! Define grid parameters
  double precision, parameter :: xcut = 7.77d0, ycut = 7.46d0
  double precision, parameter :: h = 1.d0/63.d0
  integer, parameter :: nx = 490, ny = 470
  integer, parameter :: idim = (nx+2)*(ny+2)
  
  ! Arrays to store the precomputed table
  double precision, allocatable, dimension(:) :: wtreal, wtimag
  double precision :: hrecip
  integer :: kstep
  
contains

  !---------------------------------------------------------------------
  ! Function: Faddeeva
  ! Purpose: Compute the Faddeeva function w(z) = exp(-z^2)*erfc(-i*z)
  !          for z = xx + i*yy using a continued-fraction approximation.
  !          This routine is based on the Cernlib code provided.
  ! Inputs:
  !   xx, yy - real components of the complex argument z.
  ! Output:
  !   Returns a complex number w = wx + i*wy.
  !---------------------------------------------------------------------
  function Faddeeva(xx, yy) result(w)
    implicit none
    real(fp_kind), intent(in) :: xx, yy
    complex(fp_kind) :: w
    real(fp_kind) :: wx, wy
    integer :: n, nc, nu
    real(fp_kind) :: x, y, q, h, xl, xh, yh, tx, ty, tn, sx, sy, saux
    real(fp_kind), parameter :: cc = 1.12837916709551d0
    real(fp_kind), parameter :: xlim = 5.33d0, ylim = 4.29d0
    real(fp_kind), parameter :: fac1 = 3.2d0, fac2 = 23.0d0, fac3 = 21.0d0
    real(fp_kind), dimension(34) :: rx, ry

    x = abs(xx)
    y = abs(yy)

    if ( (y < ylim) .and. (x < xlim) ) then
       q  = (1.0d0 - y / ylim) * sqrt(1.0d0 - (x/xlim)**2)
       h  = 1.0d0 / (fac1 * q)
       nc = 7 + int(fac2 * q)
       xl = h**(1 - nc)
       xh = y + 0.5d0/h
       yh = x
       nu = 10 + int(fac3 * q)
       rx(nu+1) = 0.0d0
       ry(nu+1) = 0.0d0
       do n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = 0.5d0 * tx / tn
          ry(n) = 0.5d0 * ty / tn
       end do

       sx = 0.0d0
       sy = 0.0d0
       do n = nc, 1, -1
          saux = sx + xl
          sx = rx(n) * saux - ry(n) * sy
          sy = rx(n) * sy + ry(n) * saux
          xl = h * xl
       end do
       wx = cc * sx
       wy = cc * sy
    else
       xh = y
       yh = x
       rx(1) = 0.0d0
       ry(1) = 0.0d0
       do n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = 0.5d0 * tx / tn
          ry(1) = 0.5d0 * ty / tn
       end do
       wx = cc * rx(1)
       wy = cc * ry(1)
    end if

    if (yy < 0.0d0) then
       wx = 2.0d0 * exp(yy*yy - xx*xx) * cos(2.0d0*xx*yy) - wx
       wy = -2.0d0 * exp(yy*yy - xx*xx) * sin(2.0d0*xx*yy) - wy
       if (xx > 0.0d0) wy = -wy
    else
       if (xx < 0.0d0) wy = -wy
    end if

    w = cmplx(wx, wy, kind=fp_kind)
  end function Faddeeva



  !---------------------------------------------------------------------
  ! Subroutine: InitializeFaddeevaTable
  ! Purpose: Precompute and store values of w(z) on a uniform grid.
  ! The grid spans x=0 to (nx+1)*h and y=0 to (ny+1)*h.
  !---------------------------------------------------------------------
  subroutine InitializeFaddeevaTable()
    implicit none
    integer :: i, j, k
    double precision :: x, y
    complex(fp_kind) :: w
    ! Set reciprocal of h and grid stride
    hrecip = 1.d0 / h
    kstep = nx + 2
    allocate(wtreal(idim))
    allocate(wtimag(idim))
    k = 0
    do j = 0, ny+1
       do i = 0, nx+1
          k = k + 1
          x = dble(i) * h
          y = dble(j) * h
          w = Faddeeva(x, y)
          wtreal(k) = real(w,fp_kind)
          wtimag(k) = imag(w)
       end do
    end do
  end subroutine InitializeFaddeevaTable


  !---------------------------------------------------------------------
  ! Function: FastFaddeeva
  ! Purpose: Return the complex error function w(z) for given x and y
  ! using the precomputed table via nearest-neighbor lookup.
  !---------------------------------------------------------------------
  function FastFaddeeva(x, y) result(w)
    implicit none
    double precision, intent(in) :: x, y
    complex(8) :: w
    integer :: i, j, index
    ! Ensure x and y are nonnegative for the table (our grid is for x>=0, y>=0)
    if (x < 0.d0 .or. y < 0.d0) then
       print *, "FastFaddeeva: x and y must be nonnegative for table lookup."
       stop
    end if
    ! Determine the grid indices. Here we use nearest-neighbor.
    i = int(x/h + 0.5d0)
    j = int(y/h + 0.5d0)
    if (i < 0) i = 0
    if (i > nx+1) i = nx+1
    if (j < 0) j = 0
    if (j > ny+1) j = ny+1
    ! For Fortran arrays (1-indexed), compute the index.
    index = i + j*(nx+2) + 1
    w = cmplx(wtreal(index), wtimag(index))
  end function FastFaddeeva

end module FaddeevaTable
