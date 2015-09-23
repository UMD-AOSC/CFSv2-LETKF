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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_spec_mod - sea ice and SST specified from data as per GFDL climate group !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_spec_mod

use fms_mod, only: open_namelist_file, check_nml_error, close_file, &
                   stdlog, mpp_pe, mpp_root_pe, write_version_number

use time_manager_mod, only: time_type
use data_override_mod,only: data_override
use constants_mod,    only: Tfreeze

implicit none
include 'netcdf.inc'
private
public :: get_sea_surface

character(len=128), parameter :: version = '$Id$'
character(len=128), parameter :: tagname = '$Name$'

logical :: module_is_initialized = .false.

logical :: mcm_ice = .false. ! When mcm_ice=.true., ice is handled as in supersource
namelist / ice_spec_nml / mcm_ice

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_sea_surface - get SST, ice concentration and thickness from data         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_sea_surface(time, ts, cn, iceh)
type (time_type),                         intent(in)  :: time
real, dimension(:, :),                    intent(out) :: ts
real, dimension(size(ts,1),size(ts,2),2), intent(out) :: cn
real, dimension(size(ts,1),size(ts,2)),   intent(out) :: iceh

real, dimension(size(ts,1),size(ts,2))                :: sst, icec

real ::  t_sw_freeze = -1.8
integer :: ierr, io, unit

if(.not.module_is_initialized) then
   unit = open_namelist_file()
   ierr=1
   do while (ierr /= 0)
     read(unit, nml=ice_spec_nml, iostat=io, end=20)
     ierr = check_nml_error (io, 'ice_spec_nml')
   enddo
20 call close_file (unit)
   call write_version_number(version, tagname)
   if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=ice_spec_nml)
   module_is_initialized = .true.
endif

  icec = 0.0; iceh = 0.0; sst = t_sw_freeze;
  call data_override('ICE', 'sic_obs', icec, time)
  call data_override('ICE', 'sit_obs', iceh, time)
  call data_override('ICE', 'sst_obs', sst, time)

  if(mcm_ice) then
    icec = 0.0
!   TK Mod: Limit minimum non-zero sea ice thickness to 0.01m.
!           This is to eliminate some very thin but non-zero
!           sea ice thickness values, where they really should be zero
!           but have become nonzero due to spatial interpolation
!           where the input grid and model grid are not
!           EXACTLY the same.  0.01 was obtained by trial and
!           error to roughly match supersource behavior.
!           5/22/01; 8/23/01

    where (iceh < 0.01) iceh=0.0
    where (iceh>0.0)
      icec = 1.0
      sst  = t_sw_freeze
    end where
  else
    where (icec >= 0.2)
      iceh = max(iceh, 1.0)
      sst = t_sw_freeze
    else where
      icec = 0.0
      iceh = 0.0
    end where
  endif

  where (icec==0.0 .and. sst<=t_sw_freeze) sst = t_sw_freeze+1e-10

  cn(:,:,2) = icec
  cn(:,:,1) = 1-cn(:,:,2)
  ts = sst+Tfreeze

  return
end subroutine get_sea_surface

end module ice_spec_mod
