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
module ocean_bbc_mod
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matthew Harrison
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Set bottom boundary conditions 
!</OVERVIEW>
! Set bottom boundary conditions 
!<DESCRIPTION>
! 
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bbc_nml">
!  <DATA NAME="cdbot" UNITS="dimensionless" TYPE="real">
! Dimensionless coefficient for quadratic bottom drag. 
!  </DATA> 
!  <DATA NAME="u2tides" UNITS="m^2/s^2" TYPE="real">
! Residual bottom velocity due to unresolved tidal fluctuations that contribute 
! to bottom dissipation.  Should be turned off when running with explicit
! representation of tides. 
!  </DATA> 
!</NAMELIST>
!
use fms_mod,         only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_mod,         only: mpp_error, FATAL, stdout, stdlog
use mpp_domains_mod, only: mpp_update_domains, BGRID_NE

use ocean_types_mod,   only: ocean_velocity_type, ocean_domain_type
use ocean_types_mod,   only: ocean_grid_type, ocean_prog_tracer_type
use ocean_types_mod,   only: ocean_time_type
use ocean_domains_mod, only: get_local_indices

implicit none

private

#include <ocean_memory.h>

integer :: num_prog_tracers

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

public :: ocean_bbc_init
public :: get_ocean_bbc

real :: cdbot = 1.e-3
real :: u2tides = 2.5e-3

namelist /ocean_bbc_nml/ cdbot, u2tides

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_bbc_init">
!
! <DESCRIPTION>
! Initialize the bottom boundary condition module
! </DESCRIPTION>
!
subroutine ocean_bbc_init(Grid, Domain, T_prog)

type(ocean_grid_type), target, intent(in)   :: Grid
type(ocean_domain_type), target, intent(in) :: Domain
type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

integer :: n, ioun, io_status, ierr

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

ioun = open_namelist_file()
read(ioun, ocean_bbc_nml, iostat=io_status)
write (stdout(),'(/)')
write (stdout(), ocean_bbc_nml)
write (stdlog(), ocean_bbc_nml)
ierr = check_nml_error(io_status,'ocean_bbc_nml')
call close_file(ioun)

#ifndef STATIC_MEMORY
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid

num_prog_tracers = size(T_prog(:))

do n=1, num_prog_tracers
#ifndef STATIC_MEMORY
   allocate(T_prog(n)%btf(isd:ied,jsd:jed))
#endif
   T_prog(n)%btf = 0.0
enddo

return

end subroutine ocean_bbc_init
! </SUBROUTINE> NAME="ocean_bbc_init"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_bbc">
!
! <DESCRIPTION>
! Set bottom boundary conditions for velocity and tracer.
! </DESCRIPTION>
!
subroutine get_ocean_bbc(Time, Velocity, T_prog)

type(ocean_time_type), intent(in)           :: Time
type(ocean_velocity_type), intent(inout)    :: Velocity
type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

integer :: n, taum1, i, j, kz
real    :: uvmag


if (.not. module_is_initialized) then 
   call mpp_error(FATAL,'==>Error from ocean_bbc_mod (get_ocean_bbc): module must be initialized ')
endif 

do n=1, num_prog_tracers
   T_prog(n)%btf(isd:ied,jsd:jed) = 0.0
enddo

taum1 = Time%taum1

do j=jsc,jec
   do i=isc,iec
      kz = Grd%kmu(i,j)
      if (kz /= 0) then
         uvmag    = sqrt(u2tides + Velocity%u(i,j,kz,1,taum1)**2 + Velocity%u(i,j,kz,2,taum1)**2)
         Velocity%bmf(i,j,1) = cdbot*Velocity%u(i,j,kz,1,taum1)*uvmag
         Velocity%bmf(i,j,2) = cdbot*Velocity%u(i,j,kz,2,taum1)*uvmag
      endif
   enddo
enddo

call mpp_update_domains(Velocity%bmf(:,:,1),Velocity%bmf(:,:,2),Dom%domain2d,gridtype=BGRID_NE)

return

end subroutine get_ocean_bbc
! </SUBROUTINE> NAME="get_ocean_bbc"

end module ocean_bbc_mod


