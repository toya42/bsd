!=====================================================================
! Module: PosteriorStatistics
! Purpose: Compute the mode of the posterior distribution from a
!          history CSV file (for the beta=1 chain) and reconstruct
!          a ModelParameters structure.
!=====================================================================
module PosteriorStatistics
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use GlobalData
  use ModelFunctions
  implicit none
  integer, parameter :: NBINS = 50  ! Number of bins for histogram mode calculation.
contains

  !----------------------------------------------------------
  ! Function: ComputeModeFromVector
  ! Purpose: Compute the mode (most frequent value) of a real vector
  !          using a fixed-bin histogram.
  !----------------------------------------------------------
  real(fp_kind) function ComputeModeFromVector(vec)
    implicit none
    real(fp_kind), intent(in) :: vec(:)
    integer :: i, nbins
    real(fp_kind) :: xmin, xmax, bin_width, value
    integer, allocatable :: counts(:)
    integer :: bin_idx, max_idx, count_max

    nbins = NBINS
    xmin = minval(vec)
    xmax = maxval(vec)
    if (xmax == xmin) then
       ComputeModeFromVector = xmin
       return
    end if
    bin_width = (xmax - xmin) / real(nbins, fp_kind)
    allocate(counts(nbins))
    counts = 0

    do i = 1, size(vec)
       value = vec(i)
       bin_idx = int((value - xmin) / bin_width) + 1
       if (bin_idx > nbins) bin_idx = nbins
       counts(bin_idx) = counts(bin_idx) + 1
    end do

    max_idx = 1
    count_max = counts(1)
    do i = 2, nbins
       if (counts(i) > count_max) then
          count_max = counts(i)
          max_idx = i
       end if
    end do

    ComputeModeFromVector = xmin + (real(max_idx, fp_kind) - 0.5d0)*bin_width
    deallocate(counts)
  end function ComputeModeFromVector

  !----------------------------------------------------------
  ! Subroutine: ReadHistoryFile
  ! Purpose: Read a CSV file (without header) containing the history 
  !          of parameters from the beta=1 chain into a 2D array.
  !----------------------------------------------------------
  subroutine ReadHistoryFile(filename, history_data, n_history, n_columns)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp_kind), allocatable, dimension(:,:) :: history_data
    integer, intent(out) :: n_history, n_columns
    character(len=256) :: line
    integer :: unit, ios, count, i, j
    character(len=20), allocatable, dimension(:) :: tokens
    real(fp_kind), allocatable :: temp_row(:)
    integer :: num_fields

    unit = 101
    ! First pass: count the number of data lines and determine n_columns.
    n_history = 0
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error reading history file: ", trim(filename)
       stop
    end if
    n_columns = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       n_history = n_history + 1
       if (n_history == 1) then
          call Tokenize(line, tokens, num_fields)
          n_columns = num_fields
          deallocate(tokens)
       end if
    end do
    close(unit)

    allocate(history_data(n_history, n_columns))
    ! Second pass: read the data.
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    count = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       count = count + 1
       call Tokenize(line, tokens, num_fields)
       if (num_fields /= n_columns) then
          print *, "Error: Inconsistent number of columns in history file."
          stop
       end if
       allocate(temp_row(n_columns))
       do j = 1, n_columns
          read(tokens(j), *) temp_row(j)
       end do
       history_data(count, :) = temp_row(:)
       deallocate(temp_row)
       deallocate(tokens)
    end do
    close(unit)
  end subroutine ReadHistoryFile

  !----------------------------------------------------------
  ! Subroutine: Tokenize
  ! Purpose: Split a string line into tokens separated by commas.
  !          Each token is limited to 20 characters.
  !----------------------------------------------------------
  subroutine Tokenize(line, tokens, num_tokens)
    implicit none
    character(len=*), intent(in) :: line
    character(len=20), allocatable, dimension(:), intent(out) :: tokens
    integer, intent(out) :: num_tokens
    integer :: pos, start, len_line, token_count, i
    character(len=20) :: token

    len_line = len_trim(line)
    token_count = 0
    start = 1
    pos = 1
    allocate(tokens(0))
    do while (pos <= len_line)
       if (line(pos:pos) == ',') then
          token = adjustl(line(start:pos-1))
          token = token(1:min(len_trim(token),20))  ! Limit to 20 characters.
          token_count = token_count + 1
          call append_token(tokens, token, token_count)
          start = pos + 1
       end if
       pos = pos + 1
    end do
    ! Add last token.
    if (start <= len_line) then
       token = adjustl(line(start:len_line))
       token = token(1:min(len_trim(token),20))
       token_count = token_count + 1
       call append_token(tokens, token, token_count)
    end if
    num_tokens = token_count
  end subroutine Tokenize

  !----------------------------------------------------------
  ! Subroutine: append_token
  ! Purpose: Append a new token to an allocatable array of tokens.
  !----------------------------------------------------------
  subroutine append_token(tokens, new_token, new_size)
    implicit none
    character(len=20), allocatable, dimension(:), intent(inout) :: tokens
    character(len=*), intent(in) :: new_token
    integer, intent(in) :: new_size
    character(len=20), allocatable, dimension(:) :: temp
    integer :: i

    if (.not. allocated(tokens)) then
       allocate(tokens(new_size))
       tokens(new_size) = new_token
    else
       allocate(temp(new_size))
       do i = 1, new_size - 1
          temp(i) = tokens(i)
       end do
       temp(new_size) = new_token
       deallocate(tokens)
       tokens = temp
    end if
  end subroutine append_token

  !----------------------------------------------------------
  ! Subroutine: ComputePosteriorMode
  ! Purpose: Compute the mode for each parameter from the history file
  !          (for the beta=1 chain) and reconstruct a ModelParameters structure.
  !----------------------------------------------------------
  subroutine ComputePosteriorMode(history_filename, theta_mode)
    implicit none
    character(len=*), intent(in) :: history_filename
    type(ModelParameters), intent(out) :: theta_mode
    real(fp_kind), allocatable, dimension(:,:) :: history_data
    integer :: n_history, n_columns, j, pos
    real(fp_kind), allocatable, dimension(:) :: mode_vector

    ! Expected n_columns = 1 + 5 + 4 + 4*K1 + 4*K2.
    n_columns = 1 + 5 + 4 + 4*K1 + 4*K2

    call ReadHistoryFile(history_filename, history_data, n_history, n_columns)

    allocate(mode_vector(n_columns))
    ! For each column, compute the mode.
    do j = 1, n_columns
       mode_vector(j) = ComputeModeFromVector(history_data(:, j))
    end do

    ! Reconstruct theta_mode from mode_vector.
    pos = 2   ! Skip iteration number in column 1.
    theta_mode%step%Ba    = mode_vector(pos); pos = pos + 1
    theta_mode%step%Bb    = mode_vector(pos); pos = pos + 1
    theta_mode%step%H     = mode_vector(pos); pos = pos + 1
    theta_mode%step%E0    = mode_vector(pos); pos = pos + 1
    theta_mode%step%Gamma = mode_vector(pos); pos = pos + 1

    theta_mode%WL%A_WL      = mode_vector(pos); pos = pos + 1
    theta_mode%WL%mu_WL      = mode_vector(pos); pos = pos + 1
    theta_mode%WL%sigma_G_WL = mode_vector(pos); pos = pos + 1
    theta_mode%WL%gamma_L_WL = mode_vector(pos); pos = pos + 1

    do j = 1, K1
       theta_mode%low(j)%A       = mode_vector(pos); pos = pos + 1
       theta_mode%low(j)%mu      = mode_vector(pos); pos = pos + 1
       theta_mode%low(j)%sigma_G = mode_vector(pos); pos = pos + 1
       theta_mode%low(j)%gamma_L = mode_vector(pos); pos = pos + 1
    end do

    do j = 1, K2
       theta_mode%high(j)%A       = mode_vector(pos); pos = pos + 1
       theta_mode%high(j)%mu      = mode_vector(pos); pos = pos + 1
       theta_mode%high(j)%sigma_G = mode_vector(pos); pos = pos + 1
       theta_mode%high(j)%gamma_L = mode_vector(pos); pos = pos + 1
    end do

    deallocate(mode_vector, history_data)
  end subroutine ComputePosteriorMode

end module PosteriorStatistics
