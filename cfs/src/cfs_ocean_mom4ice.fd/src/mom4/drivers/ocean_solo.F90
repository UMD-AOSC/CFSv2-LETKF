!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison
!</CONTACT>
!
!<REVIEWER EMAIL="V.Balaji@noaa.gov"> V. Balaji
!</REVIEWER>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Driver for ocean-only simulations.
!</OVERVIEW>
!
!<DESCRIPTION>
! Driver for the ocean-only simulations. Similar to the FMS coupler, but 
! allows one to run the ocean model without compiling  other models. 
! Much simpler than the full FMS coupler. 
! </DESCRIPTION>
!
! <NAMELIST NAME="ocean_solo_nml">
!
!   <DATA NAME="date_init"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. If the restart file
!      ocean_solo.res is present, date_init will be taken from there.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
!     'no_calendar' can not be used because the time_manager's date  function are used. 
!     
!   </DATA>
!   <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
!     The number of months that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
!     The number of days that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
!     The number of hours that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
!     The number of minutes that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
!     The number of seconds that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="dt_ocean"  TYPE="integer"  DEFAULT="0">
!     Ocean model time step in seconds. 
!   </DATA>
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds for coupling between ocean and atmospheric models: 
!     must be an integral multiple of dt_ocean. This is the "slow" timestep.
!     Note that for an ocean_solo model, the coupling to an "atmosphere" is the coupling 
!     to some data files.  In this case, dt_cpld represents the time which data is updated.
!     For example, if data is "daily", then dt_cpld=86400 should be chosen.  
!     If data is fixed, then dt_cpld of any integer of dt_ocean is fine, with
!     dt_cpld=86400 the default. 
!   </DATA>
!
! </NAMELIST>
!
!   <NOTE>
!     <PRE>
!     1.The actual run length will be the sum of months, 
!       days, hours, minutes, and seconds. A run length of zero
!       is not a valid option. 
!     2.The run length must be an integral multiple of the coupling 
!       timestep dt_cpld. 
!     </PRE>
!   </NOTE>
!
  use constants_mod,            only: constants_init
  use data_override_mod,        only: data_override_init, data_override
  use diag_manager_mod,         only: diag_manager_init, register_diag_field, diag_manager_end
  use field_manager_mod,        only: field_manager_init
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: mpp_error, FATAL, mpp_pe, mpp_npes, mpp_set_current_pelist
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_time, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
  use time_manager_mod,         only: JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )

  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_types_mod,          only: ocean_data_type, ice_ocean_boundary_type
#ifdef PRISM
  use ocean_prism_mod,          only: ocean_prism_init, localComm, ocean_prism_terminate
#endif
  implicit none

  type (ocean_data_type)                 :: Ocean_sfc          
  type(ice_ocean_boundary_type), target  :: Ice_ocean_boundary 

  ! define some time types 
  type(time_type) :: Time_init    ! initial time for experiment
  type(time_type) :: Time_start   ! start time for experiment
  type(time_type) :: Time_end     ! end time for experiment (as determined by dtts)
  type(time_type) :: Run_len      ! length of experiment 
  type(time_type) :: Time        
  type(time_type) :: Time_step_ocean    
  type(time_type) :: Time_step_coupled

  character(len=17) :: calendar = 'julian'

  integer :: dt_ocean = 3600
  integer :: dt_cpld  = 86400
  integer :: num_cpld_calls  = 0
  integer :: num_ocean_calls = 0
  integer :: nc, no
  integer :: calendar_type=-1

  integer :: date_init(6)=0, date(6)
  integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: yr, mon, day, hr, min, sec

  integer :: isc,iec,jsc,jec
  integer :: unit, log_unit, io_status, ierr

  integer :: memuse
  integer :: flags=0, override_clock
  integer :: nfields 
  
  character(len=256) :: version = ''
  character(len=256) :: tag = ''

  logical :: ocean_seg_start
  logical :: ocean_seg_end

  character(len=9) :: month
  
  namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds, dt_cpld, dt_ocean

! initialize shared modules
#ifdef PRISM
  call ocean_prism_init(ierr)
  call fms_init(localComm)
#else
  call fms_init()
#endif
  call constants_init

  flags = MPP_CLOCK_SYNC

  ! provide for namelist over-ride
  unit = open_namelist_file('input.nml')
  read  (unit, ocean_solo_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(),'(/47x,a/)') '======== MODEL BEING DRIVEN BY OCEAN_SOLO_MOD ========'
  write (stdout(), ocean_solo_nml)  
  write (stdlog(), ocean_solo_nml)
  ierr = check_nml_error(io_status,'ocean_solo_nml')
  call close_file (unit)

  write (stdlog(),'(/,80("="),/(a))') trim(version), trim(tag)

  ! set the calendar 

  select case( uppercase(trim(calendar)) )
  case( 'JULIAN' )
     calendar_type = JULIAN
  case( 'NOLEAP' )
     calendar_type = NOLEAP
  case( 'THIRTY_DAY' )
     calendar_type = THIRTY_DAY_MONTHS
  case( 'NO_CALENDAR' )
     calendar_type = NO_CALENDAR
  case default
     call mpp_error( FATAL, 'ocean_solo: ocean_solo_nml entry calendar must be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
  end select 

  ! get ocean_solo restart : this can override settings from namelist
  if (file_exist('INPUT/ocean_solo.res')) then
      call mpp_open(unit,'INPUT/ocean_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
      read(unit,*) calendar_type 
      read(unit,*) date_init
      read(unit,*) date
      call close_file(unit)
  endif
      
  call set_calendar_type (calendar_type)

!!$ initialize pelists for ocean ensembles set current pelist to ensemble member
!!$ need to call prior to diagnostics_init
!!$ code presently not supported (mjh)
!!$  call ocean_ensemble_init() 
                             
  call field_manager_init(nfields)

  call diag_manager_init()

  call time_interp_external_init()

  if (sum(date_init) <= 0) then
      call mpp_error(FATAL,'==>Error from ocean_solo_mod: date_init must be set either in ocean_solo.res or in ocean_solo_nml')
  else
      Time_init  = set_date(date_init(1),date_init(2), date_init(3), &
           date_init(4),date_init(5),date_init(6))
  endif

  if (file_exist('INPUT/ocean_solo.res')) then
      Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
  else
      Time_start = Time_init
      date = date_init
  endif

  Time_end          = increment_date(Time_start, years, months, days, hours, minutes, seconds)
  Run_len           = Time_end - Time_start
  Time_step_ocean   = set_time(dt_ocean, 0)
  Time_step_coupled = set_time(dt_cpld, 0)
  num_cpld_calls    = Run_len / Time_step_coupled
  num_ocean_calls   = Time_step_coupled / Time_step_ocean
  Time = Time_start

  call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_APPEND,threading=MPP_SINGLE)

  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call close_file (unit)  
  
  call ocean_model_init(Ocean_sfc, Time_init, Time, Time_step_ocean )

  call data_override_init(Ocean_domain_in = Ocean_sfc%domain)

  override_clock = mpp_clock_id('Override', flags=flags,grain=CLOCK_COMPONENT)
  
  call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)
  
  allocate ( Ice_ocean_boundary% u_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% v_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% t_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% q_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% salt_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% lw_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux (isc:iec,jsc:jec), &
             Ice_ocean_boundary% lprec (isc:iec,jsc:jec), &
             Ice_ocean_boundary% fprec (isc:iec,jsc:jec), &
             Ice_ocean_boundary% runoff (isc:iec,jsc:jec), &
             Ice_ocean_boundary% calving (isc:iec,jsc:jec), &
             Ice_ocean_boundary% p (isc:iec,jsc:jec))

  Ice_ocean_boundary%u_flux    = 0.0
  Ice_ocean_boundary%v_flux    = 0.0
  Ice_ocean_boundary%t_flux    = 0.0
  Ice_ocean_boundary%q_flux    = 0.0
  Ice_ocean_boundary%salt_flux = 0.0
  Ice_ocean_boundary%lw_flux   = 0.0
  Ice_ocean_boundary%sw_flux   = 0.0
  Ice_ocean_boundary%lprec     = 0.0
  Ice_ocean_boundary%fprec     = 0.0
  Ice_ocean_boundary%runoff    = 0.0
  Ice_ocean_boundary%calving   = 0.0
  Ice_ocean_boundary%p         = 0.0

  ! loop over the coupled calls 
  do nc=1, num_cpld_calls
     
     call mpp_clock_begin(override_clock)
     call ice_ocn_bnd_from_data(Ice_ocean_boundary)
     call mpp_clock_end(override_clock)
     
     ! loop over the ocean calls 
     do no=1, num_ocean_calls
        ocean_seg_start = ( no .eq. 1 )     
        ocean_seg_end   = ( no .eq. num_ocean_calls )
        call update_ocean_model(Ice_ocean_boundary, Ocean_sfc, ocean_seg_start, ocean_seg_end, num_ocean_calls)
        Time = Time + Time_step_ocean
     enddo

  enddo

  ! close some of the main components 

  call ocean_model_end(Ocean_sfc)

  call diag_manager_end(Time)

  ! need to reset pelist before calling mpp_clock_end
  call mpp_set_current_pelist() 

    ! write restart file
    call mpp_open( unit, 'RESTART/ocean_solo.res', nohdrs=.TRUE. )
    if ( mpp_pe().EQ.mpp_root_pe() )then
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        call get_date(Time_init,yr,mon,day,hr,min,sec)
        write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
             'Model start time:   year, month, day, hour, minute, second'
        call get_date(Time_end ,yr,mon,day,hr,min,sec)
        write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
             'Current model time: year, month, day, hour, minute, second'
    end if
    call mpp_close(unit)
    call fms_io_exit
    call fms_end    
#ifdef PRISM
    call ocean_prism_terminate
#endif
  contains


!====================================================================
! get forcing data from data_overide 
  subroutine ice_ocn_bnd_from_data(x)

      type (ice_ocean_boundary_type) :: x
      integer                        :: i,j
      type(time_type)                :: Time_next

      Time_next = Time + Time_step_ocean      
      call data_override('OCN', 't_flux',    x%t_flux   , Time_next)
      call data_override('OCN', 'u_flux',    x%u_flux   , Time_next)
      call data_override('OCN', 'v_flux',    x%v_flux   , Time_next)
      call data_override('OCN', 'q_flux',    x%q_flux   , Time_next)
      call data_override('OCN', 'salt_flux', x%salt_flux, Time_next)
      call data_override('OCN', 'lw_flux',   x%lw_flux  , Time_next)
      call data_override('OCN', 'sw_flux',   x%sw_flux  , Time_next)
      call data_override('OCN', 'lprec',     x%lprec    , Time_next)
      call data_override('OCN', 'fprec',     x%fprec    , Time_next)
      call data_override('OCN', 'runoff',    x%runoff   , Time_next)
      call data_override('OCN', 'calving',   x%calving  , Time_next)
      call data_override('OCN', 'p',         x%p        , Time_next)
            
  end subroutine ice_ocn_bnd_from_data


end program main
  
        

  





