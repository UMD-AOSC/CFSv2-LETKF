program sscycle2

  use sigio_module
  use sfcio_module

  implicit none

  type(sigio_head) :: headsig
  type(sfcio_head) :: headsfc
  type(sigio_data) :: datasig
  type(sfcio_data) :: datasfc

  integer :: iret
  character(len=32) :: arg

  !------------------------------------------------------------
  
  ! get the date from the command line argument
  call get_command_argument(1, arg)
  if (len_trim(arg) /= 10) then
     write (0,*) "[Error] expected YYYYMMDDHH as argument"
     stop 1
  end if

  ! write (*,*) "Changing gfs output headers to be initial conditions at", &
  !      arg

  !! read in the current data
  call sigio_srohdc(11, 'fort.11',  headsig, datasig, iret)
  if (iret /= 0) then
     write (0,*) '[Error] open and read sigma file: fort.11'
     stop 1
  end if
  call sfcio_srohdc(12, 'fort.12',  headsfc, datasfc, iret)
  if (iret /= 0) then
     write (0,*) '[Error] open and read surface file: fort.12'
     stop 1
  end if

  !! make changes
  headsig%fhour = 0
  read (arg(1:4), *) headsig%idate(4)
  read (arg(5:6), *) headsig%idate(2)
  read (arg(7:8), *) headsig%idate(3)
  read (arg(9:10),*) headsig%idate(1)
  headsfc%fhour = 0
  headsfc%idate = headsig%idate

  call sigio_swohdc(21, 'fort.21', headsig, datasig, iret)
  if (iret /= 0) then
     write (0,*) '[Error] cant write to sigma file: fort.21'
     stop 1
  end if
  call sfcio_swohdc(22, 'fort.22', headsfc, datasfc, iret)
  if (iret /= 0) then
     write (0,*) '[Error] cant write to surface file: fort.22'
     stop 1
  end if

end program sscycle2
