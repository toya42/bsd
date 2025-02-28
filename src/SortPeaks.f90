!=====================================================================
! Module: SortPeaks
! Purpose: Provide a bubble sort routine to order an array of 
!          ParametersPeak (used for low and high peak groups) based on mu.
!=====================================================================
module SortPeaks
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use ModelFunctions, only: ParametersPeak
  implicit none
contains

  !----------------------------------------------------------
  ! Subroutine: BubbleSortPeaks
  ! Purpose: Sort an array of ParametersPeak in place according to the
  !          mu field. Optionally, sort in ascending (default) or descending order.
  !
  ! Inputs:
  !   peaks     - Array of ParametersPeak to be sorted.
  !   ascending - (Optional) Logical flag indicating sort order.
  !               Default is .true. for ascending.
  !----------------------------------------------------------
  subroutine BubbleSortPeaks(peaks, ascending)
    implicit none
    type(ParametersPeak), intent(inout) :: peaks(:)
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: n, i
    type(ParametersPeak) :: temp
    logical :: swapped

    if (present(ascending)) then
       asc = ascending
    else
       asc = .true.
    end if

    n = size(peaks)
    do
       swapped = .false.
       do i = 1, n - 1
          if (asc) then
             if (peaks(i)%mu > peaks(i+1)%mu) then
                temp = peaks(i)
                peaks(i) = peaks(i+1)
                peaks(i+1) = temp
                swapped = .true.
             end if
          else
             if (peaks(i)%mu < peaks(i+1)%mu) then
                temp = peaks(i)
                peaks(i) = peaks(i+1)
                peaks(i+1) = temp
                swapped = .true.
             end if
          end if
       end do
       if (.not. swapped) exit
    end do
  end subroutine BubbleSortPeaks

  subroutine BubbleSortBetas(betas, ascending)
    implicit none
    real(fp_kind), intent(inout) :: betas(:)
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: n, i
    real(real64) :: temp
    logical :: swapped

    if (present(ascending)) then
       asc = ascending
    else
       asc = .true.
    end if

    n = size(betas)
    do
       swapped = .false.
       do i = 1, n - 1
          if (asc) then
             if (betas(i) > betas(i+1)) then
                temp = betas(i)
                betas(i) = betas(i+1)
                betas(i+1) = temp
                swapped = .true.
             end if
          else
             if (betas(i) < betas(i+1)) then
                temp = betas(i)
                betas(i) = betas(i+1)
                betas(i+1) = temp
                swapped = .true.
             end if
          end if
       end do
       if (.not. swapped) exit
    end do
  end subroutine BubbleSortBetas


end module SortPeaks
