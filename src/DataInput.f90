!=====================================================================
! Module: DataInput
! Purpose: Read experimental data from CSV and MCMC parameters from a text file.
!=====================================================================
module DataInput
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use GlobalData
  implicit none
contains

  !----------------------------------------------------------
  ! Subroutine: ReadExperimentalData
  ! Purpose: Read experimental data from a CSV file.
  !   The file has a header in the first line.
  !   Columns: 
  !       1: Energy, 2: Measurement time, 
  !       3: Measurement proportional to I_inc, 4: Measurement proportional to I_ab, etc.
  !   A proportionality constant of 1e5 is applied to the measured counts.
  !   The subroutine reads the file name from standard input.
  !----------------------------------------------------------
  subroutine ReadExperimentalData()
    implicit none
    character(len=256) :: filename
    character(len=256) :: line
    integer :: unit, ios, count_lines, i
    real(fp_kind) :: energy_val, time_val, meas_inc, meas_ab
    real(fp_kind), parameter :: scale = 1.0d3
    ! First pass: count the number of data lines (excluding header)
    unit = 10
    print *, "Enter experimental data CSV filename:"
    read(*, '(A)') filename
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file ", trim(filename)
       stop
    end if

    ! Skip the header line.
    read(unit, '(A)', iostat=ios) line
    count_lines = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       count_lines = count_lines + 1
    end do
    close(unit)

    ! Set N from the file.
    N = count_lines
    print *, "Number of measurement points (N) = ", N

    ! Allocate arrays
    allocate(E(N))
    allocate(I_inc(N))
    allocate(I_ab(N))

    ! Second pass: read the data.
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error reopening file ", trim(filename)
       stop
    end if

    ! Skip header
    read(unit, '(A)', iostat=ios) line

    i = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       i = i + 1
       ! Assume CSV format: use internal read with comma delimiter.
       ! For simplicity, we assume the columns are separated by commas.
       read(line, *) energy_val, time_val, meas_inc, meas_ab
       E(i) = energy_val
       I_inc(i) = meas_inc * scale
       I_ab(i) = meas_ab
       !print *, energy_val, time_val, meas_inc, meas_ab
    end do
    close(unit)
  end subroutine ReadExperimentalData

  !----------------------------------------------------------
  ! Subroutine: ReadMCMCParameters
  ! Purpose: Read MCMC parameters from a text file (param_mcmc.txt).
  !   The file format is as follows (one parameter per line):
  !       10000 T_iter
  !        2000 T_burn
  !          50 L_rep
  !           4 K1
  !           4 K2
  !----------------------------------------------------------
  subroutine ReadMCMCParameters()
    implicit none
    character(len=256) :: paramFilename
    character(len=256) :: dummy
    integer :: unit, ios

    unit = 20
    print *, "Enter MCMC parameters filename:"
    read(*, '(A)') paramFilename

    open(unit, file=trim(paramFilename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file ", trim(paramFilename)
       stop
    end if

    read(unit, *) T_iter, dummy  ! Expect line: e.g., 10000 T_iter
    read(unit, *) T_burn, dummy  ! e.g., 2000 T_burn
    read(unit, *) L_rep, dummy   ! e.g., 50 L_rep
    read(unit, *) K1, dummy      ! e.g., 4 K1
    read(unit, *) K2, dummy      ! e.g., 4 K2

    close(unit)

    if(L_rep>L_rep_max) then
       print *, "Error Replica size, L_rep must be smaller than ", L_rep_max+1
       stop
    end if


    print *, "MCMC parameters:"
    print *, "  T_iter =", T_iter
    print *, "  T_burn =", T_burn
    print *, "  L_rep =", L_rep
    print *, "  K1 =", K1
    print *, "  K2 =", K2

    ! Allocate beta array with size L_rep and initialize (e.g., linear or geometric progression)
    allocate(beta(L_rep))
    ! For example, a simple linear progression from 0.1 to 1.0:
    do unit = 1, L_rep
       beta(unit) = 0.1d0 + (1.0d0 - 0.1d0) * real(unit-1, fp_kind) / real(L_rep-1, fp_kind)
    end do
  end subroutine ReadMCMCParameters

end module DataInput
