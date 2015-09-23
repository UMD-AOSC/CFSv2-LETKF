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
module ocean_velocity_diag_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Numerical diagnostics for velocity related quantities. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Numerical diagnostics for velocity related quantities.
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_velocity_diag_nml">
!  <DATA NAME="diag_freq" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!</NAMELIST>
!
use constants_mod,    only: grav
use diag_manager_mod, only: register_diag_field, need_data, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,          only: FATAL, stdout, stdlog
use mpp_domains_mod,  only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,          only: mpp_error, mpp_max, mpp_sum
use time_manager_mod, only: time_type, increment_time

use ocean_bih_friction_mod, only: horz_bih_viscosity_check, horz_bih_reynolds_check
use ocean_domains_mod,      only: get_local_indices
use ocean_lap_friction_mod, only: horz_lap_viscosity_check, horz_lap_reynolds_check
use ocean_operators_mod,    only: FAX, FAY
use ocean_types_mod,        only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,        only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_types_mod,        only: ocean_velocity_type, ocean_density_type
use ocean_util_mod,         only: write_timestamp

implicit none

private

#include <ocean_memory.h>

real    :: dtuv

logical :: used

! for kinetic energy 
integer :: id_ke_tot=-1
real    :: ke_tot

! for potential energy 
integer :: id_pe_tot=-1
logical :: potential_energy_first=.true.
real, dimension(:), allocatable :: potential_0 ! potential energy (joules) at initial timestep. 

! for output
integer :: unit=6

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

real :: test_nml=0.0
integer :: global_sum_flag      ! flag for mpp_global_sum

integer :: diag_freq = -1
logical :: do_bitwise_exact_sum = .false.

public  ocean_velocity_diag_init
public  ocean_velocity_diagnostics 
public  kinetic_energy
public  velocity_change
public  potential_energy
private velocity_land_cell_check

namelist /ocean_velocity_diag_nml/ diag_freq, do_bitwise_exact_sum


contains

!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_velocity_diag module containing subroutines
! diagnosing velocity related properties of the simulation.  These are 
! not terms in the equations, but rather they are diagnosed from 
! terms. 
! </DESCRIPTION>
!
subroutine ocean_velocity_diag_init(Grid, Domain, Time, Time_steps, Velocity)

type(ocean_grid_type), target, intent(in)   :: Grid
type(ocean_domain_type), target, intent(in) :: Domain
type(ocean_time_type), intent(in)           :: Time
type(ocean_time_steps_type), intent(in)     :: Time_steps
type(ocean_velocity_type), intent(in)       :: Velocity

integer :: n, ioun, io_status, ierr

if (module_is_initialized) return

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

ioun = open_namelist_file()
read(ioun, ocean_velocity_diag_nml, iostat=io_status)
write (stdlog(), ocean_velocity_diag_nml)
write (stdout(),'(/)')
write (stdout(), ocean_velocity_diag_nml)
ierr = check_nml_error(io_status,'ocean_velocity_diag_nml')
call close_file(ioun)

if (diag_freq == 0) diag_freq = 1

if(do_bitwise_exact_sum) then
   global_sum_flag = BITWISE_EXACT_SUM
else
   global_sum_flag = NON_BITWISE_EXACT_SUM
endif

#ifndef STATIC_MEMORY
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid

dtuv = Time_steps%dtuv

! kinetic energy 
id_ke_tot = register_diag_field ('ocean_model', 'ke_tot', Time%model_time, &
                                 'Total ke', '10^-15 joules',missing_value=0.0, range=(/0.0,1000.0/))
! potential energy 
id_pe_tot = register_diag_field ('ocean_model', 'pe_tot', &
            Time%model_time, 'pe(time)-pe(first)', '10^-15 joules',&
            missing_value=0.0, range=(/0.0,1000.0/))

allocate (potential_0(nk))

return

end subroutine ocean_velocity_diag_init
! </SUBROUTINE>  NAME="ocean_velocity_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the velocity. 
! </DESCRIPTION>
!
subroutine ocean_velocity_diagnostics(Time, Thickness, Velocity, Dens)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  type(ocean_density_type), intent(in)   :: Dens
  type(time_type)                        :: next_time 
  real                                   :: pe_tot

  ! diagnostics that do not send_data to diag_manager  
  if (diag_freq > 0) then 
     if(mod(Time%itt,diag_freq) == 0) then
        call horz_lap_viscosity_check
        call horz_lap_reynolds_check(Time, Velocity)
        call horz_bih_viscosity_check
        call horz_bih_reynolds_check(Time, Velocity)
        call velocity_land_cell_check(Time, Velocity)
        call write_timestamp(Time%model_time)
        call kinetic_energy(Time, Thickness, Velocity, Dens, ke_tot, .false., .true.)
        call potential_energy(Time, Thickness, Dens, pe_tot, .false., .true.)
     endif  
  endif

  ! diagnostics that send_data to diag_manager 
  next_time = increment_time(Time%model_time, int(dtuv), 0)
  if (need_data(id_ke_tot,next_time)) then 
    call kinetic_energy(Time, Thickness, Velocity, Dens, ke_tot, .true., .false.)
  endif 
  if (need_data(id_pe_tot,next_time)) then 
    call potential_energy(Time, Thickness, Dens, pe_tot, .true., .false.)
  endif 

end subroutine ocean_velocity_diagnostics
! </SUBROUTINE>  NAME="ocean_velocity_diagnostics"


!#######################################################################
! <SUBROUTINE NAME="potential_energy">
!
! <DESCRIPTION>
! Compute potential energy (Joules) relative to initial time step. 
! </DESCRIPTION>
!
subroutine potential_energy (Time, Thickness, Dens, pe_tot, diag_flag, write_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_density_type), intent(in)   :: Dens
  real, intent(out)                      :: pe_tot
  logical, intent(in), optional          :: diag_flag 
  logical, intent(in), optional          :: write_flag 
  logical                                :: send_diagnostics
  logical                                :: write_diagnostics

  real, dimension(nk) :: potential
  integer :: i, j, k, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_velocity_diag_mod(potential_energy): module needs initialization ')
  endif 

  tau = Time%tau

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! assume write_diagnostics=.true., unless write_flag says it is false.
  if (PRESENT(write_flag)) then 
    write_diagnostics = write_flag
  else 
    write_diagnostics = .true.
  endif 

  potential(:) = 0.0
  do k=1,nk
     potential(k) = mpp_global_sum(Dom%domain2d, &
                    Grd%dat(:,:)*Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k)*Dens%rho(:,:,k)*Thickness%ztp(:,:,k),&
                    global_sum_flag) 
     potential(k) = potential(k)*grav
  enddo
  if(potential_energy_first) then 
    potential_energy_first=.false.
    do k=1,nk
       potential_0(k) = potential(k)
    enddo
  endif 

  ! computing pe_tot relative to initial time step value reduces roundoff errors
  pe_tot = 0.0
  do k=1,nk
    pe_tot = pe_tot + (potential(k)-potential_0(k))
  enddo

  if(write_diagnostics) then 
     write(stdout(),'(/,1x,a,e20.12)') 'Potential energy relative to first time (Joules) = ', pe_tot
  endif 

  if(send_diagnostics) then
    if (id_pe_tot > 0) used = send_data (id_pe_tot, pe_tot/1.e15, Time%model_time)
  endif 

end subroutine potential_energy
! </SUBROUTINE>  NAME="potential_energy"


!#######################################################################
! <SUBROUTINE NAME="kinetic_energy">
!
! <DESCRIPTION>
! Compute global integrated horizontal kinetic energy.
! </DESCRIPTION>
!
subroutine kinetic_energy (Time, Thickness, Velocity, Dens, ke_tot, diag_flag, write_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  type(ocean_density_type), intent(in)   :: Dens
  real, intent(out)                      :: ke_tot
  logical, intent(in), optional          :: diag_flag 
  logical, intent(in), optional          :: write_flag 
  logical                                :: send_diagnostics
  logical                                :: write_diagnostics

  real, dimension(isd:ied,jsd:jed) :: rho_u_cell
  real                             :: mass_u_cell, vol_u_cell
  integer                          :: i, j, k, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_velocity_diag_mod(kinetic_energy): module needs initialization ')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! assume write_diagnostics=.true., unless write_flag says it is false.
  if (PRESENT(write_flag)) then 
    write_diagnostics = write_flag
  else 
    write_diagnostics = .true.
  endif 

  tau = Time%tau

  ! total horizontal kinetic energy(joules) at "tau"

  ke_tot = 0.0
  do k=1,nk
    rho_u_cell(:,:) = FAX(FAY(Dens%rho(:,:,k)))
    do j=jsc,jec
      do i=isc,iec
        vol_u_cell  = Grd%dau(i,j)*Thickness%dhu(i,j,k,tau)*Grd%umask(i,j,k)
        mass_u_cell = rho_u_cell(i,j)*vol_u_cell
        ke_tot      = ke_tot + 0.5*mass_u_cell*(Velocity%u(i,j,k,1,tau)**2 + Velocity%u(i,j,k,2,tau)**2)
      enddo
    enddo
  enddo
  call mpp_sum (ke_tot)

  if(write_diagnostics) then 
     write(stdout(),'(1x,a,e20.12)') 'Kinetic energy  (Joules)                         = ', ke_tot
  endif 

  if(send_diagnostics) then 
    if (id_ke_tot > 0) used = send_data (id_ke_tot, ke_tot/1.e15, Time%model_time)
  endif 

end subroutine kinetic_energy
! </SUBROUTINE> NAME="kinetic_energy"


!#######################################################################
! <SUBROUTINE NAME="velocity_land_cell_check">
!
! <DESCRIPTION>
! See if there are any points over land with nonzero ocean velocity 
! </DESCRIPTION>
!
subroutine velocity_land_cell_check(Time, Velocity)

  implicit none
  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  integer                               :: i, j, k, num, taup1

  taup1 = Time%taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_velocity_diag_mod(velocity_land_cell_check): module needs initialization ')
  endif 

  call write_timestamp(Time%model_time)
  write (stdout(),'(//60x,a/)') ' Land cell summary:'
  write (stdout(),'(1x,a/)')'Locations (if any) where land cell velocity is non-zero...'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) == 0 .and. (Velocity%u(i,j,k,1,taup1) /= 0.0 .or. Velocity%u(i,j,k,2,taup1) /= 0.0)) then
          num = num + 1
          write (unit,9000) i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), Velocity%u(i,j,k,1,taup1), Velocity%u(i,j,k,2,taup1)
        endif
      enddo
    enddo
  enddo

  call mpp_max(num)
  if (num > 0) call mpp_error(FATAL,'==>Error: found nonzero ocean velocity over land points. ')

9000  format(/' " =>Error: Land cell at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m) has u=',e10.3, ' v=',e10.3)

end subroutine velocity_land_cell_check
! </SUBROUTINE> NAME="velocity_land_cell_check"



!#######################################################################
! <SUBROUTINE NAME="velocity_change">
!
! <DESCRIPTION>
! Determine the number of points that have large single-time step 
! changes in the abs of the velocity.  
! </DESCRIPTION>
!
subroutine velocity_change(Time, Velocity,velocity_change_max,velocity_change_max_num)

  implicit none
  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  real, intent(in)                      :: velocity_change_max
  integer, intent(in)                   :: velocity_change_max_num
  real                                  :: abschange(2)
  integer                               :: i, j, k, n, num
  integer                               :: taum1, tau, taup1 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_velocity_diag_mod(velocity_change): module needs initialization ')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  call write_timestamp(Time%model_time)
  write (stdout(),'(//60x,a/)') ' Velocity change summary (test for leap-frog noise):'
  write (stdout(),'(1x,a,e12.6/)')'Locations (if any) where abs(0.5*(u(taup1)+u(taum1))-u(tau)) (m/s) > ',velocity_change_max
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%umask(i,j,k) > 0) then 
           do n=1,2
             abschange(n) = abs(0.5*(Velocity%u(i,j,k,n,taup1)+Velocity%u(i,j,k,n,taum1))-Velocity%u(i,j,k,n,tau))
           enddo
           if (abschange(1) > velocity_change_max .or. abschange(2) > velocity_change_max) then
             num = num + 1
             write (unit,9000) i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), abschange(1), abschange(2)
           endif 
        endif
      enddo
    enddo
  enddo

  call mpp_max(num)
  if (num > velocity_change_max_num) then 
    call mpp_error(FATAL,'Found too many points where velocity changed too much over single time step.')
  endif 

9000  format(/' " =>Error: Ocean at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m) has excessive change(u)=',e10.3, ' or excessive change(v)=',e10.3)

end subroutine velocity_change
! </SUBROUTINE> NAME="velocity_change"


end module ocean_velocity_diag_mod


