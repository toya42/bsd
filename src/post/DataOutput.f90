!=====================================================================
! Module: DataOutput
! Purpose: Record the history of model parameters during MCMC and
!          export the experimental and restored spectrum.
!=====================================================================
module DataOutput
  use, intrinsic :: iso_fortran_env
  use precision, only : fp_kind
  use GlobalData
  use ModelFunctions
  implicit none

  integer(int32), parameter :: BUFFER_SIZE = 1000
  integer(int32) :: history_length
  ! Explanation:
  ! 1: iteration number
  ! 5: step parameters [Ba, Bb, H, E0, Gamma]
  ! 4: WL parameters [A_WL, mu_WL, sigma_G_WL, gamma_L_WL]
  ! 4*K1: low-energy peaks (each 4 parameters)
  ! 4*K2: high-energy peaks (each 4 parameters)

  ! Buffer for storing history entries (each row is one entry)
  real(fp_kind), allocatable, dimension(:,:) :: history_buffer
  integer(int32) :: buffer_count = 0

  integer(int32), private :: history_unit = 99       ! File unit for history output.
  integer(int32), private :: spectrum_unit = 100       ! File unit for spectrum export.
  logical :: history_initialized = .false.
  logical :: spectrum_initialized = .false.
  
  character(len=64) :: fmt1,fmt2

contains

  !----------------------------------------------------------
  ! Subroutine: InitializeHistoryOutput
  ! Purpose: Open the history output file and initialize the buffer.
  ! Input:
  !   filename: Name of the output CSV file for history.
  !----------------------------------------------------------
  subroutine InitializeHistoryOutput(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer(int32) :: iostat_local,idx,i
    character(len=20),dimension(8) :: c_temp
    character(len=20),dimension(10+4*K1+4*K2) :: header
    character(len=50) :: c_idx, c_pn


    history_length = 1 + 5 + 4 + 4*K1 + 4*K2

    !open(unit = history_unit, file = trim(filename), status='replace', action='write', iostat=iostat_local)
    !if (iostat_local /= 0) then
    !   print *, "Error opening history output file: ", trim(filename)
    !   stop
    !end if
    !history_initialized = .true.
    !allocate(history_buffer(BUFFER_SIZE, history_length))
    !buffer_count = 0
    !! Write CSV header.

    write(fmt1,'(I0)') history_length
    fmt1 = '('//trim(fmt1)//'(a20))'
    write(fmt2,'(I0)') history_length
    fmt2 = '('//trim(fmt2)//'(E20.8e3))'
    !print *, fmt1,fmt2

    !header(1) = '1_iteration'
    !! step
    !header( 2) = '2_Step_Ba'
    !header( 3) = '3_Step_Bb'
    !header( 4) = '4_Step_H'
    !header( 5) = '5_Step_E0'
    !header( 6) = '6_Step_Gamma'
    !! wl
    !header( 7) = '7_WL_A'
    !header( 8) = '8_WL_mu'
    !header( 9) = '9_WL_Sigma'
    !header(10) = '10_WL_Gamma'

    !c_temp(1) = '_low_A'
    !c_temp(2) = '_low_mu'
    !c_temp(3) = '_low_Sigma'
    !c_temp(4) = '_low_Gamma'
    !c_temp(5) = '_high_A'
    !c_temp(6) = '_high_mu'
    !c_temp(7) = '_high_Sigma'
    !c_temp(8) = '_high_Gamma'


    !idx = 11
    !do i=1,K1
    !  write(c_pn,'(I0)') i
    !  write(c_idx,'(I0)') idx
    !  header(idx  ) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(1))
    !  write(c_idx,'(I0)') idx+1
    !  header(idx+1) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(2))
    !  write(c_idx,'(I0)') idx+2
    !  header(idx+2) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(3))
    !  write(c_idx,'(I0)') idx+3
    !  header(idx+3) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(4))
    !  idx = idx+4
    !end do
    !do i=1,K2
    !  write(c_pn,'(I0)') i
    !  write(c_idx,'(I0)') idx
    !  header(idx  ) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(5))
    !  write(c_idx,'(I0)') idx+1
    !  header(idx+1) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(6))
    !  write(c_idx,'(I0)') idx+2
    !  header(idx+2) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(7))
    !  write(c_idx,'(I0)') idx+3
    !  header(idx+3) = trim(c_idx)//"_"//trim(c_pn)//trim(c_temp(8))
    !  idx = idx+4
    !end do

    !!write(history_unit,fmt1) header(:)

  end subroutine InitializeHistoryOutput

  !----------------------------------------------------------
  ! Subroutine: AppendHistoryEntry
  ! Purpose: Append a history entry (current theta parameters) for a given iteration.
  !----------------------------------------------------------
  subroutine AppendHistoryEntry(iteration, theta)
    implicit none
    integer, intent(in) :: iteration
    type(ModelParameters), intent(in) :: theta
    real(fp_kind), dimension(history_length) :: entry
    integer :: pos, i

    pos = 1
    entry(pos) = real(iteration, fp_kind)
    pos = pos + 1

    ! Step parameters: [Ba, Bb, H, E0, Gamma]
    entry(pos:pos+4) = [ theta%step%Ba, theta%step%Bb, theta%step%H, theta%step%E0, theta%step%Gamma ]
    pos = pos + 5

    ! White-line parameters: [A_WL, mu_WL, sigma_G_WL, gamma_L_WL]
    entry(pos:pos+3) = [ theta%WL%A_WL, theta%WL%mu_WL, theta%WL%sigma_G_WL, theta%WL%gamma_L_WL ]
    pos = pos + 4

    ! Low-energy peaks: each has 4 parameters.
    do i = 1, K1
       entry(pos:pos+3) = [ theta%low(i)%A, theta%low(i)%mu, theta%low(i)%sigma_G, theta%low(i)%gamma_L ]
       pos = pos + 4
    end do

    ! High-energy peaks: each has 4 parameters.
    do i = 1, K2
       entry(pos:pos+3) = [ theta%high(i)%A, theta%high(i)%mu, theta%high(i)%sigma_G, theta%high(i)%gamma_L ]
       pos = pos + 4
    end do

    ! Add the entry to the buffer.
    buffer_count = buffer_count + 1
    history_buffer(buffer_count, :) = entry

    ! If buffer is full, flush it.
    if (buffer_count == BUFFER_SIZE) then
       call FlushHistoryBuffer()
    end if
  end subroutine AppendHistoryEntry

  !----------------------------------------------------------
  ! Subroutine: FlushHistoryBuffer
  ! Purpose: Write the buffered history entries to file and reset the buffer.
  !----------------------------------------------------------
  subroutine FlushHistoryBuffer()
    implicit none
    integer :: i

    if (.not. history_initialized) then
       print *, "History output file not initialized."
       stop
    end if

    do i = 1, buffer_count
       ! Write one history line as CSV.
       write(history_unit,fmt2) history_buffer(i, :)
    end do
    buffer_count = 0
  end subroutine FlushHistoryBuffer

  !----------------------------------------------------------
  ! Subroutine: FinalizeHistoryOutput
  ! Purpose: Flush any remaining history entries and close the file.
  !----------------------------------------------------------
  subroutine FinalizeHistoryOutput()
    implicit none
    if (history_initialized) then
       if (buffer_count > 0) call FlushHistoryBuffer()
       close(history_unit)
       history_initialized = .false.
    end if
  end subroutine FinalizeHistoryOutput

  !----------------------------------------------------------
  ! Subroutine: InitializeSpectrumOutput
  ! Purpose: Open the output file for spectrum export.
  ! Input:
  !   filename: Name of the spectrum output CSV file.
  !----------------------------------------------------------
  subroutine InitializeSpectrumOutput(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: iostat_local
    open(unit = spectrum_unit, file = trim(filename), status='replace', action='write', iostat=iostat_local)
    if (iostat_local /= 0) then
       print *, "Error opening spectrum output file: ", trim(filename)
       stop
    end if
    spectrum_initialized = .true.
    ! Write header for spectrum file.
    write(spectrum_unit, '(A)') 'Energy, I_inc, I_ab, f_ratio, Restored_Spectrum, step, wl, low, high'
  end subroutine InitializeSpectrumOutput

  !----------------------------------------------------------
  ! Subroutine: ExportRestoredSpectrum
  ! Purpose: For post-processing, export the experimental spectrum and 
  !          the restored spectrum computed using the estimated parameters.
  !          Columns: Energy, I_inc, I_ab, f_ratio(E, theta_est)
  ! Input:
  !   theta_est: Estimated ModelParameters (e.g., from beta=1 chain).
  !----------------------------------------------------------
  subroutine ExportRestoredSpectrum(theta_est)
    implicit none
    type(ModelParameters) :: theta_est
    integer(int32) :: i,enum
    real(fp_kind) :: restored_val,step,wl,low,high,fr,emin,emax,ene,I_inc_ave,e0
    if (.not. spectrum_initialized) then
       print *, "Spectrum output file not initialized."
       stop
    end if

    do i = 1, N
       step = f_step(E(i), theta_est%step) *I_inc(i)
       wl = f_WL(E(i), theta_est%WL, theta_est%step%E0)*I_inc(i)
       low = f_low(E(i), theta_est%low, theta_est%step%E0)*I_inc(i)
       high = f_high(E(i), theta_est%high, theta_est%step%E0)*I_inc(i)
       fr = f_ratio(E(i), theta_est)
       restored_val = fr*I_inc(i)
       !print *,E(i),restored_val
       write(spectrum_unit, '(10E18.8e3)') E(i), I_inc(i), I_ab(i), fr, restored_val,step,wl,low,high, step+wl
    end do

   e0 = theta_est%step%E0
   write(16,*) '==peak position=='
   write(16,*) 'step'
   write(16,*) e0
   write(16,*) 'WL'
   write(16,'(4E18.8e3)') e0+theta_est%WL%mu_WL, theta_est%WL%A_WL, theta_est%WL%sigma_G_WL, theta_est%WL%gamma_L_WL
   write(16,*) 'low'
   do i=1,K1
      write(16,*) i
      write(16,'(4E18.8e3)') e0-theta_est%low(i)%mu, theta_est%low(i)%A, theta_est%low(i)%sigma_G, theta_est%low(i)%gamma_L
   end do
   write(16,*) 'high'
   do i=1,K2
      write(16,*) i
      write(16,'(4E18.8e3)') e0+theta_est%high(i)%mu, theta_est%high(i)%A, theta_est%high(i)%sigma_G, theta_est%high(i)%gamma_L
   end do

   emin = E(1)
   emax = E(N)
   enum = 10000
   I_inc_ave = sum(I_inc)/size(I_inc)
   do i=1,enum
      ene = (emax-emin)/real(enum)*real(i)+emin
      step = f_step(ene, theta_est%step) *I_inc_ave
      wl = f_WL(ene, theta_est%WL, theta_est%step%E0)*I_inc_ave
      low = f_low(ene, theta_est%low, theta_est%step%E0)*I_inc_ave
      high = f_high(ene, theta_est%high, theta_est%step%E0)*I_inc_ave
      fr = f_ratio(ene, theta_est)
      restored_val = fr*I_inc_ave
      !print *,E(i),restored_val
      write(17, '(7E18.8e3)') ene, fr, restored_val,step,wl,low,high
   end do
  

  end subroutine ExportRestoredSpectrum

  !----------------------------------------------------------
  ! Subroutine: FinalizeSpectrumOutput
  ! Purpose: Close the spectrum output file.
  !----------------------------------------------------------
  subroutine FinalizeSpectrumOutput()
    implicit none
    if (spectrum_initialized) then
       close(spectrum_unit)
       spectrum_initialized = .false.
    end if
  end subroutine FinalizeSpectrumOutput

end module DataOutput
