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
! Modified: Xingren Wu
!           Xingren.Wu@noaa.gov
module ocean_tracer_util_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski
!</CONTACT>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov">
! S. M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! This module contains many routines of use for tracers in mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Tracer utility module for mom4. 
!</DESCRIPTION>
!
use mpp_mod,             only: stdout, stdlog, FATAL
use mpp_mod,             only: mpp_error, mpp_chksum, mpp_pe, mpp_min, mpp_max
use platform_mod,        only: i8_kind
  
use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,     only: ocean_time_type
use ocean_util_mod,      only: write_timestamp
use ocean_workspace_mod, only: wrk1 

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6


#include <ocean_memory.h>

logical :: module_is_initialized = .false.

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_tracer_util_init
public tracer_min_max
public tracer_prog_chksum 
public tracer_diag_chksum 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_util_init">
!
! <DESCRIPTION>
! Initialize mom4 tracer utilities.
! </DESCRIPTION>
!
subroutine ocean_tracer_util_init (Grid, Domain)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain

  if (module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_tracer_util_mod (ocean_tracer_util_init): module already initialized')
  endif 

  module_is_initialized = .true.

  write( stdlog(),'(/a/)') trim(version)

  Dom => Domain
  Grd => Grid

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif

end subroutine ocean_tracer_util_init
! </SUBROUTINE> NAME="ocean_tracer_util_init">



!#######################################################################
! <SUBROUTINE NAME="tracer_min_max">
!
! <DESCRIPTION>
! Compute the global min and max for tracers.            
!
! Vectorized using maxloc() and minloc() intrinsic functions by 
! Russell.Fiedler@csiro.au.
!
! </DESCRIPTION>
!
subroutine tracer_min_max(Time, Tracer)
  
  type(ocean_time_type), intent(in)        :: Time
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  real    :: tmax, tmin, xtmax, xtmin, ytmax, ytmin, ztmax, ztmin, tmax0, tmin0
  integer :: i,j,k, itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: tau 
  real    :: fudge
  
  ! arrays to enable vectorization
  integer :: iminarr(3),imaxarr(3)

  if (.not. module_is_initialized) then
     call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod (tracer_min_max): ocean_tracer_util_mod needs initialization ')
  endif 

  tmax=-1.e10;tmin=1.e10
  itmax=0;jtmax=0;ktmax=0
  itmin=0;jtmin=0;ktmin=0

  tau = Time%tau

  call write_timestamp(Time%model_time)
  
  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,tau)

  if(ANY(Grd%tmask(isc:iec,jsc:jec,:) > 0.)) then
     iminarr=minloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     imaxarr=maxloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     itmin=iminarr(1)+isc-1
     jtmin=iminarr(2)+jsc-1
     ktmin=iminarr(3)
     itmax=imaxarr(1)+isc-1
     jtmax=imaxarr(2)+jsc-1
     ktmax=imaxarr(3)
  end if
  
  if(ktmin >0) then
     tmin=wrk1(itmin,jtmin,ktmin)
     tmax=wrk1(itmax,jtmax,ktmax)
  end if

!!$  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,tau)

!!$  do k=1,nk
!!$    do j=jsc,jec
!!$       do i=isc,iec
!!$           if (wrk1(i,j,k) > tmax .and. Grd%tmask(i,j,k) > 0.) then
!!$              tmax = wrk1(i,j,k)
!!$              itmax = i;jtmax=j;ktmax=k
!!$           endif
!!$           if (wrk1(i,j,k) < tmin .and. Grd%tmask(i,j,k) > 0.) then
!!$              tmin = wrk1(i,j,k)
!!$              itmin = i;jtmin=j;ktmin=k
!!$           endif
!!$        enddo
!!$     enddo
!!$  enddo

  ! use "fudge" to distinguish processors when tracer extreme is independent of processor
  fudge = 1.0 + 1.e-12*mpp_pe() 
  tmax = tmax*fudge
  tmin = tmin*fudge
  if(tmax == 0.0) then 
    tmax = tmax + 1.e-12*mpp_pe() 
  endif 
  if(tmin == 0.0) then 
    tmin = tmin + 1.e-12*mpp_pe() 
  endif 
  

  tmax0=tmax;tmin0=tmin
  xtmax = Grd%xt(itmax,jtmax)
  xtmin = Grd%xt(itmin,jtmin)
  ytmax = Grd%yt(itmax,jtmax)
  ytmin = Grd%yt(itmin,jtmin)
  ztmax = Grd%zt(ktmax)
  ztmin = Grd%zt(ktmin)

  call mpp_max(tmax)
  call mpp_min(tmin)

  if (tmax0 == tmax) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The maximum T is ', tmax,&
               ' deg C at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod: The maximum temperature is outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The maximum S is ', tmax,&
               ' PSU at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod: The maximum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The maximum '//trim(Tracer%name)// ' is ', tmax,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,','&
               ,ktmax,'),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod: The maximum tracer is outside allowable range')
          endif 
      endif
  endif
  
  if (tmin0 == tmin) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The minimum T is ', tmin,&
               ' deg C at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod (tracer_min_max): minimum temp outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The minimum S is ', tmin,&
               ' PSU at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod: The minimum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es24.17,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') &
               ' The minimum '//trim(Tracer%name)// ' is ', tmin,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,','&
               ,ktmin,'),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod (tracer_min_max): minimum tracer outside allowable range')
          endif
      endif
  endif

  return


end subroutine tracer_min_max
! </SUBROUTINE>  NAME="tracer_min_max"


!#######################################################################
! <SUBROUTINE NAME="tracer_prog_chksum">
!
! <DESCRIPTION>
! Compute checksums for prognostic tracers 
! </DESCRIPTION>
subroutine tracer_prog_chksum(Time, Tracer, index, chksum)

  type(ocean_time_type), intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(in)  :: Tracer
  integer, intent(in)                       :: index
  integer(i8_kind), optional, intent(inout) :: chksum
  integer(i8_kind)                          :: chk_sum
  
  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod (tracer_prog_chksum): module not yet initialized ')
  endif 

  write(stdout(),*) '=== Prognostic tracer checksum follows ==='
  write(stdout(),*) 'Tracer name = ', Tracer%name

  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,index)*Grd%tmask(isc:iec,jsc:jec,:)

  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  
  write(stdout(),*) 'Tracer chksum = ',  chk_sum

  if (PRESENT(chksum)) chksum = chk_sum

end subroutine tracer_prog_chksum
! </SUBROUTINE>  NAME="tracer_prog_chksum"


!#######################################################################
! <SUBROUTINE NAME="tracer_diag_chksum">
!
! <DESCRIPTION>
! Compute checksums for diagnostic tracers 
! </DESCRIPTION>
subroutine tracer_diag_chksum(Time, Tracer, chksum)

  type(ocean_time_type), intent(in)         :: Time
  type(ocean_diag_tracer_type), intent(in)  :: Tracer
  integer(i8_kind), optional, intent(inout) :: chksum
  integer(i8_kind)                          :: chk_sum
  
  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_util_mod (tracer_diag_chksum): module not yet initialized ')
  endif 

  write(stdout(),*) '=== Diagnostic tracer checksum follows ==='
  write(stdout(),*) 'Tracer name = ', Tracer%name

  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)

  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  
  write(stdout(),*) 'Tracer chksum = ',  chk_sum

  if (PRESENT(chksum)) chksum = chk_sum

end subroutine tracer_diag_chksum
! </SUBROUTINE>  NAME="tracer_diag_chksum"


end module ocean_tracer_util_mod

