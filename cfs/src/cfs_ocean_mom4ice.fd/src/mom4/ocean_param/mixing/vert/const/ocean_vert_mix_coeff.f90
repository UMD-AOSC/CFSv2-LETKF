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
module ocean_vert_mix_coeff_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! Constant vertical viscosity and diffusivity 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes a time independent vertical viscosity and diffusivity. 
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_mix_coeff_const_nml">
!  <DATA NAME="kappa_h" UNITS="m^2/sec" TYPE="real">
!  The constant vertical diffusivity.  Used for cases when wanting a space-time
!  independent diffusivity.  The "h" is historical and stands for "heat".
!  </DATA> 
!  <DATA NAME="kappa_m" UNITS="m^2/sec" TYPE="real">
!  The constant vertical viscosity.  Used for cases when wanting a space-time
!  independent viscosity.  
!  </DATA> 
!  <DATA NAME="diff_cbt_limit" UNITS="m^2/sec" TYPE="real">
!  The largest allowable vertical diffusivity.  Of use for cases where vertically unstable
!  columns are stabilized with a large vertical diffusivity.  
!  </DATA> 
!</NAMELIST>

use constants_mod,      only: pi
use fms_mod,            only: write_version_number, FATAL, stdout, stdlog
use fms_mod,            only: open_namelist_file, check_nml_error, close_file
use mpp_io_mod,         only: mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII
use mpp_mod,            only: mpp_error

use ocean_density_mod,  only: density_delta_z
use ocean_domains_mod,  only: get_local_indices
use ocean_types_mod,    only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,    only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,    only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_vert_mix_mod, only: diff_cbt_back


implicit none

public ocean_vert_mix_coeff_init
public vertical_mix_coeff

private

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: const/ocean_vert_mix_coeff.F90 ($Id$)'
character (len=128) :: tagname = &
     '$Name$'

integer :: index_temp, index_salt, num_prog_tracers, num_diag_tracers
logical :: module_is_initialized = .FALSE.

real :: kappa_h        = 0.1e-4  ! constant vertical diffusivity (m^2/sec)
real :: kappa_m        = 1.0e-4  ! constant vertical viscosity (m^2/sec)
real :: diff_cbt_limit = 1.0     ! diffusivity used to vertically adjust (m^2/sec)

namelist /ocean_vert_mix_coeff_const_nml/ kappa_h, kappa_m, diff_cbt_limit

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_coeff_init">
!
! <DESCRIPTION>
! Initialize the constant vertical diffusivity module.
! </DESCRIPTION>
!
subroutine ocean_vert_mix_coeff_init (Grid, Domain, Time, Time_steps, T_prog, T_diag)

  type(ocean_grid_type), intent(in), target              :: Grid
  type(ocean_domain_type), intent(in), target            :: Domain
  type(ocean_time_type), intent(in)                      :: Time
  type(ocean_time_steps_type), intent(in)                :: Time_steps 
  type(ocean_prog_tracer_type), intent(in), dimension(:) :: T_prog
  type(ocean_diag_tracer_type), intent(in), dimension(:) :: T_diag

  integer :: k, n, ioun, ierr, io_status
  real    :: dzmin

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_mix_coeff_mod (ocean_vert_mix_coeff_init): module is already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_vert_mix_coeff_const_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_vert_mix_coeff_const_nml)  
  write (stdlog(),ocean_vert_mix_coeff_const_nml)
  ierr = check_nml_error(io_status,'ocean_vert_mix_coeff_const_nml')
  call close_file (ioun)

  if(kappa_h==0.0) then 
    write(stdout(),'(/a)')'==>USING constant vertical diffusivity with kappa_h = 0.0.'
  else 
    write(stdout(),'(/a)')'==>USING constant vertical diffusivity with kappa_h > 0.0.'
  endif
  if(kappa_m==0.0) then 
    write(stdout(),'(a/)')'==>USING constant vertical viscosity with kappa_m = 0.0.'
  else 
    write(stdout(),'(a/)')'==>USING constant vertical viscosity with kappa_m > 0.0.'
  endif

  write(stdout(),'(/a,f10.2)')'==>Note from ocean_vert_mix_coeff_mod: using forward time step for vert-frict of (secs)', &
                                  Time_steps%dtime_u 
  write(stdout(),'(/a,f10.2)')'==>Note from ocean_vert_mix_coeff_mod: using forward time step for vert-diff  of (secs)', &
                                  Time_steps%dtime_t 

  Dom => Domain
  Grd => Grid

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  do k=1,nk
    if ((Time_steps%dtime_t*kappa_h)/Grid%dzt(k)**2 >= 0.5) then
      write (stdout(),'(a,a,i3)')'==> Warning: vertical diffusive criteria exceeded for "kappa_h".',&
                         ' use a smaller "dtts" and/or "kappa_h" at level k=',k
    endif
  enddo
  dzmin  = 1.e10  ! meters 
  do k=1,nk
    dzmin = min(dzmin,Grid%dzt(k))
  enddo
  if ((Time_steps%dtime_u*kappa_m)/dzmin**2 >= 0.5) then
    write (stdout(),'(a,a)')'==> Warning: vertical diffusive criteria exceeded on "kappa_m".',&
                       ' use a smaller "dtuv" and/or "kappa_m"              '
  endif

  num_prog_tracers = size(T_prog)
  num_diag_tracers = size(T_diag)

  index_temp = -1
  index_salt = -1
  
  do n = 1, num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
     if (trim(T_prog(n)%name) == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,'==>Error in ocean_vert_mix_coeff_mod: temp or salt not present in tracer array')
  endif 

end subroutine ocean_vert_mix_coeff_init
! </SUBROUTINE>  NAME="ocean_vert_mix_coeff_init"


!#######################################################################
! <SUBROUTINE NAME="vertical_mix_coeff">
!
! <DESCRIPTION>
! This function computes the vertical diffusivity and viscosity.  
! These mixing coefficients are time independent but generally 
! arbitrary functions of space. 
! </DESCRIPTION>
!
  subroutine vertical_mix_coeff(aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, swflx, sw_frac, pme, &
                                river, visc_cbu, diff_cbt, hblt_depth)

  real, intent(in)                                     :: aidif
  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_thickness_type), intent(in)               :: Thickness
  type(ocean_velocity_type), intent(in)                :: Velocity
  type(ocean_prog_tracer_type), intent(in)             :: T_prog(num_prog_tracers)
  type(ocean_diag_tracer_type), intent(in)             :: T_diag(num_diag_tracers)
  type(ocean_density_type), intent(in)                 :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed)         :: swflx
  real, intent(in), dimension(isd:ied,jsd:jed,0:nk)    :: sw_frac
  real, intent(in), dimension(isd:ied,jsd:jed)         :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)         :: river
  real, intent(in), dimension(isd:ied,jsd:jed)         :: hblt_depth
  real, intent(inout), dimension(isd:ied,jsd:jed,nk)   :: visc_cbu
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt

  real, dimension(isd:ied,jsd:jed,nk) ::  delta_rho
  integer                             :: i, j, k
  integer                             :: tau

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_mix_coeff (vertical_mix_coeff): module must be initialized')
  endif 

  do k=1,nk
    visc_cbu(isc:iec,jsc:jec,k)   = kappa_m
    diff_cbt(isc:iec,jsc:jec,k,:) = kappa_h
  enddo

  tau = Time%tau

  ! set vertical diffusivity to diff_cbt_limit where gravitationally unstable
  if(aidif==1.0) then
      
    delta_rho(:,:,:) = density_delta_z (Dens%rho(:,:,:), &
         T_prog(index_salt)%field(:,:,:,tau), &
         T_prog(index_temp)%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:))
    
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
          if (delta_rho(i,j,k)*Grd%tmask(i,j,k) > 0.0) then
            diff_cbt(i,j,k,1) = diff_cbt_limit*Grd%tmask(i,j,k+1) 
          endif
        enddo
      enddo
    enddo
  endif

  do k=1,nk-1
     do i=isc,iec
        do j=jsc,jec
           diff_cbt(i,j,k,1) = max(diff_cbt(i,j,k,1),diff_cbt_back(i,j,k))
        enddo
     enddo
  enddo
  diff_cbt(isc:iec,jsc:jec,:,2) = diff_cbt(isc:iec,jsc:jec,:,1)

end subroutine vertical_mix_coeff
! </SUBROUTINE> NAME="vertical_mix_coeff"

end module ocean_vert_mix_coeff_mod
