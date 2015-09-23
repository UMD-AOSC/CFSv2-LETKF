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
module ocean_lap_friction_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes the thickness weighted time tendency for  
! horizontal velocity arising from horizontal Laplacian friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for
! horizontal velocity arising from horizontal Laplacian friction. 
! The viscosity used to determine the strength of the tendency 
! can be a general function of space and time as specified by 
! the Smagorinsky approach as well as a grid-scale dependent
! background viscosity.  The form of the friction operator 
! can be isotropic or anisotropic in the horizontal plane. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies and R.W. Hallberg, 
! Biharmonic friction with a Smagorinsky viscosity for use in large-scale
! eddy-permitting ocean models
! Monthly Weather Review vol 128 (2000) pages 2935--2946
! </REFERENCE>
!
! <REFERENCE>
! R. D. Smith and J. C. McWilliams, 
! Anisotropic viscosity for ocean models
! Ocean Modelling submitted 2002 
! </REFERENCE>
!
! <NOTE>
! The ocean model can generally run with both Laplacian and biharmonic friction
! enabled at the same time.  Such has been found useful for some eddying 
! ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lap_friction_general_nml">
!  <DATA NAME="lap_friction_on" TYPE="logical">
!  Must be true to use this module. 
!  </DATA> 
!
!  <DATA NAME="k_smag_iso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="k_smag_aniso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky anisotropic viscosity. 
!  </DATA> 
!  <DATA NAME="viscosity_ncar" TYPE="logical">
!  Anisotropic background viscosities used by NCAR. 
!  </DATA> 
!  <DATA NAME="vel_micom_iso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_aniso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity. 
!  </DATA> 
!
!  <DATA NAME="equatorial_zonal" TYPE="logical">
!  Orient the anisotropic friction within a latitudinal band according to zonal direction. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal_lat" TYPE="real">
!  Latitudinal band to use the zonal friction orientation. 
!  </DATA> 
!  <DATA NAME="ncar_isotropic_off_equator" TYPE="logical">
!  Polewards of equatorial_zonal_lat, revert NCAR scheme to isotropic 
!  </DATA> 
!  <DATA NAME="equatorial_no_smag" TYPE="logical">
!  Turn smag off within equatorial_zonal_lat region. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_iso" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity within
!  a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_aniso" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity within
!  a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity is set according 
!  to eq_vel_micom_iso and eq_vel_micom_aniso.
!  </DATA> 
!
!  <DATA NAME="restrict_polar_visc" TYPE="logical">
!  For restricting the background viscosity poleward of a 
!  latitude.  This method may be useful for coupling to an ice model
!  in which case the horizontal viscosity may need to be a bit 
!  smaller to maintain time step constraints.  This is because the 
!  effective friction is larger than that just within the ocean.  
!  </DATA> 
!  <DATA NAME="restrict_polar_visc_lat" TYPE="real">
!  Latitude poleward of which we restrict the viscosity.
!  </DATA> 
!  <DATA NAME="restrict_polar_visc_ratio" TYPE="real">
!  Ratio of the normal critical value that we limit the 
!  viscosity to be no greater than.  If restrict_polar_visc_ratio=1.0
!  then there is no special limitation of the viscosity beyond that 
!  of the one-dimensional stability constraint.  
!  </DATA> 
!
!  <DATA NAME="bottom_5point" TYPE="logical">
!  To alleviate problems with small partial cells, it is often necessary to reduce the 
!  operator to the traditional 5-point Laplacian at the ocean bottom.  This logical 
!  implements this mixing. 
!  </DATA> 
!
!  <DATA NAME="neptune" TYPE="logical">
!  Set to true for computing friction relative to Neptune barotropic velocity. 
!  </DATA> 
!  <DATA NAME="neptune_length_eq" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at equator.  
!  </DATA> 
!  <DATA NAME="neptune_length_pole" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at pole. 
!  </DATA> 
!
!  <DATA NAME="vconst_1" UNITS="cm^2/sec" TYPE="real">
!  Background viscosity for NCAR algorithm.
!  </DATA> 
!  <DATA NAME="vconst_2" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_3" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_4" UNITS="1/cm" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_5" TYPE="integer">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_6" UNITS="cm^2/sec" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_7" UNITS="cm/sec">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="debug_ncar_A" TYPE="logical">
!  Sets f_perp=f_para for debugging purposes with the NCAR scheme.
!  </DATA> 
!  <DATA NAME="debug_ncar_B" TYPE="logical">
!  Sets f_para=f_perp for debugging purposes with the NCAR scheme.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: pi, radius, epsln, omega, radian
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, write_version_number, close_file
use mpp_domains_mod,  only: mpp_update_domains, BGRID_NE, mpp_global_field
use mpp_mod,          only: mpp_sum, mpp_pe, mpp_error, mpp_max
use mpp_mod,          only: FATAL, NOTE, stdout, stdlog
use mpp_mod,          only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: BAY, BAX, BDX_EU, BDY_NU, FDX_PU, FDY_PU, FMX, FMY
use ocean_operators_mod, only: FAY, FAX, FDX_NT, FDY_ET 
use ocean_types_mod,     only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_velocity_type


implicit none

private

public ocean_lap_friction_init
public horz_lap_friction
public horz_lap_viscosity_check
public horz_lap_reynolds_check

private BDX_EU_smag
private BDY_NU_smag
private compute_neptune_velocity

! clock ids
integer :: id_lap_friction
 
! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_aiso=-1
integer :: id_aaniso=-1
integer :: id_along_back=-1
integer :: id_across_back=-1
integer :: id_neptune_u=-1
integer :: id_neptune_v=-1
integer :: id_neptune_psi=-1
integer :: id_neptune_ft=-1
integer :: id_along=-1
integer :: id_across=-1
integer :: id_crit=-1
integer :: id_lap_fric_u=-1
integer :: id_lap_fric_v=-1
logical :: used

real :: k_smag_iso          = 4.0        ! smag scaling coeff for isotropic viscosity (dimensionless) 
real :: k_smag_aniso        = 0.0        ! smag scaling coeff for anisotripic viscosity (dimensionless) 
real :: vel_micom_iso       = 0.0        ! background scaling velocity for isotropic viscosity (m/sec) 
real :: vel_micom_aniso     = 0.0        ! background scaling velocity for anisotropic viscosity (m/sec) 
real :: eq_vel_micom_iso    = 0.0        ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real :: eq_vel_micom_aniso  = 0.0        ! background scaling velocity (m/sec) within equatorial band for anisotropic visc
real :: eq_lat_micom          = 0.0      ! equatorial latitude band for micom (degrees)
real :: equatorial_zonal_lat  = 0.0      ! latitudinal band for orienting the friction along zonal direction
logical :: equatorial_zonal   = .false.  ! zonally orient anisotropic friction w/i equatorial_zonal_lat
logical :: equatorial_no_smag = .false.  ! remove smag within equatorial_zonal_lat
logical :: bottom_5point      = .true.   ! for bottom Laplacian 5point mixing to avoid problems with thin partial cells. 

! for restricting critical viscosity in high latitudes.
! this has been found of use for coupling to ice models, where 
! the effective friction is larger than just that within the ocean.   
logical :: restrict_polar_visc=.false.  
real    :: restrict_polar_visc_lat=60.0
real    :: restrict_polar_visc_ratio=0.35
                                 
! From NCAR viscosity algorithm.  Defaults are from CCSM2.0 ocean component gx1v3 
integer :: nx, ny                    ! number of global points in generalized x and y directions            
real :: vconst_1            = 1.e7   ! (cm^2/s)
real :: vconst_2            = 0.0    ! (dimensionless)
real :: vconst_3            = 0.16   ! (dimensionless)
real :: vconst_4            = 2.e-8  ! (1/cm) 
integer :: vconst_5         = 3      ! (number of grid points)
real :: vconst_6            = 1.e7   ! (cm^2/sec)
real :: vconst_7            = 100.0  ! (cm/sec)
logical :: viscosity_ncar   = .false.! to get the (x,y,z) dependent isotropic and anisotropic viscosities used by NCAR.  
logical :: debug_ncar_A     = .false.
logical :: debug_ncar_B     = .false.
logical :: ncar_isotropic_off_equator = .false.
logical :: ncar_only_equatorial = .false.  ! for use of the ncar scheme only within the tropical band 

! for Holloway's Neptune scheme  
real    :: neptune_length_eq   = 1.2e3   ! (metres)
real    :: neptune_length_pole = 3.0e3   ! (metres)
logical :: neptune             = .false. ! for computing friction relative to barotropic Neptune velocity

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
  
real, dimension(isd:ied,jsd:jed)    :: fsmag_iso   ! (m^2) combination of terms for computing isotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: fsmag_aniso ! (m^2) combination of terms for computing anisotropic smag visc
real, dimension(isd:ied,jsd:jed,nk) :: aiso_back   ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk) :: aaniso_back ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk) :: aiso        ! isotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: aaniso      ! anisotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed)    :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed)    :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed,nk) :: volqc       ! (1/4)*dxu*dyu*dhu (m^3) for time-independent quarter-cell volume
real, dimension(isd:ied,jsd:jed)    :: visc_crit   ! critical value of the viscosity (m^2/sec) for linear stability 

real, dimension(isd:ied,jsd:jed,0:1,0:1,2) :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdelx
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdely
real, dimension(isd:ied,jsd:jed)           :: tmp
real, dimension(isc:iec,jsc:jec)           :: cos2theta
real, dimension(isc:iec,jsc:jec)           :: sin2theta

! for the Neptune scheme 
real, dimension(isd:ied,jsd:jed,2)         :: neptune_velocity ! barotropic velocity (m/s) from Neptune 

#else

real, dimension(:,:), allocatable   :: fsmag_iso    ! (m^2) combination of terms for computing isotropic smag visc
real, dimension(:,:), allocatable   :: fsmag_aniso  ! (m^2) combination of terms for computing anisotropic smag visc
real, dimension(:,:,:), allocatable :: aiso_back    ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(:,:,:), allocatable :: aaniso_back  ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(:,:,:), allocatable :: aiso         ! isotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(:,:,:), allocatable :: aaniso       ! anisotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(:,:), allocatable   :: daur_dxur    ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(:,:), allocatable   :: daur_dyur    ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(:,:,:), allocatable :: volqc        ! (1/4)*dxu*dyu*dhu (m^3) for time-independent quarter-cell volume
real, dimension(:,:), allocatable   :: visc_crit    ! critical value of the viscosity (m^2/sec) for linear stability 

real, dimension(:,:,:,:,:), allocatable        :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(:,:,:), allocatable            :: tmpfdelx
real, dimension(:,:,:), allocatable            :: tmpfdely
real, dimension(:,:), allocatable              :: tmp
real, dimension(:,:), allocatable              :: cos2theta
real, dimension(:,:), allocatable              :: sin2theta

! for the Neptune scheme 
real, dimension(:,:,:), allocatable            ::  neptune_velocity   ! barotropic velocity (m/s) from Neptune 

#endif

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()


character(len=256) :: version=&
     '=>Using: /general/ocean_lap_friction.f90 ($Id$)'
character (len=128) :: tagname = &
     '$Name$'

logical :: lap_friction_on       = .false.
logical :: module_is_initialized = .FALSE.

namelist /ocean_lap_friction_general_nml/ lap_friction_on, bottom_5point, k_smag_iso, k_smag_aniso, &
                         vel_micom_iso, vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom, &
                         equatorial_zonal, equatorial_zonal_lat, equatorial_no_smag, &
                         viscosity_ncar, ncar_isotropic_off_equator, ncar_only_equatorial, &
                         vconst_1, vconst_2, vconst_3, vconst_4, vconst_5, vconst_6, vconst_7, &
                         debug_ncar_A, debug_ncar_B, &
                         neptune, neptune_length_eq, neptune_length_pole, &
                         restrict_polar_visc, restrict_polar_visc_lat, restrict_polar_visc_ratio 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_lap_friction_init">

! <DESCRIPTION>
! Initialize the horizontal Laplacian friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_lap_friction_init(Grid, Domain, Time, d_time)

  type(ocean_grid_type), target, intent(in)   :: Grid
  type(ocean_domain_type), target, intent(in) :: Domain
  type(ocean_time_type), intent(in)           :: Time
  real, intent(in)                            :: d_time

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real :: coeff_iso, coeff_aniso, dxdy

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod(ocean_lap_friction_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_lap_friction_general_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_lap_friction_general_nml)  
  write (stdlog(),ocean_lap_friction_general_nml)
  ierr = check_nml_error(io_status,'ocean_lap_friction_general_nml')
  call close_file(ioun)

#ifndef STATIC_MEMORY

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(lap_friction_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_lap_friction_mod from lap/general.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_lap_friction_mod.')
    return
  endif 

  dtime = d_time
  write(stdout(),'(/a,f10.2)')'==> Note from ocean_lap_friction_mod: using forward time step of (secs)', dtime 

  if(bottom_5point) then 
    write(stdout(),'(1x,a)') '==> Note: Will reduce horizontal friction &
         &to a 5point Laplacian on the bottom'
    write(stdout(),'(5x,a)') 'This helps to alleviate numerical problems&
         & with thin bottom partial cells.'
  endif 

  if(neptune) then 
    write(stdout(),'(1x,a)') '==> Note: Computing friction relative to Neptune'
    write(stdout(),'(5x,a,e10.5,a)') 'Equatorial length for computing Neptune velocity is ',neptune_length_eq,' meters' 
    write(stdout(),'(5x,a,e10.5,a)') 'Polar length for computing Neptune velocity is      ',neptune_length_pole,' meters' 
  endif 

  if(viscosity_ncar)  then 
    write( stdout(),'(1x,a)')'==> NOTE: USING background horz viscosities &
         &according to NCAR CCSM2.0 algorithm.' 
    write( stdout(),'(a,e15.4)')'  NCAR vconst_1 (cm^2/sec)  = ',vconst_1
    write( stdout(),*)' NCAR vconst_2             = ',vconst_2
    write( stdout(),*)' NCAR vconst_3             = ',vconst_3
    write( stdout(),*)' NCAR vconst_4 (1/cm)      = ',vconst_4
    write( stdout(),*)' NCAR vconst_5             = ',vconst_5
    write( stdout(),'(a,e15.4)')'  NCAR vconst_6 (cm^2/sec)  = ',vconst_6
    write( stdout(),*)' NCAR vconst_7 (cm/sec)    = ',vconst_7
    if(debug_ncar_A) then 
      write( stdout(),'(1x,a)')'==> NOTE: USING debug_ncar_A=.true., &
           &so will set f_perp=f_para'
    endif 
    if(debug_ncar_B) then 
      write( stdout(),'(1x,a)')'==> NOTE: USING debug_ncar_B=.true.,&
           & so will set f_para=f_perp'
    endif 
    if(ncar_isotropic_off_equator) then 
      write( stdout(),'(1x,a,f10.4)')'==> ncar_isotropic_off_equator=.true. =>f_para=f_perp poleward of lat',equatorial_zonal_lat
    endif 
    if(ncar_only_equatorial) then 
      write( stdout(),'(1x,a,f10.4)')'==> ncar_only_equatorial=.true. =>ncar scheme only in band +/- lat',equatorial_zonal_lat
    endif 

  endif  

  if(k_smag_iso > 0.0) then 
     write( stdout(),'(1x,a)')'==> NOTE: USING horz isotropic viscosity &
          &via Smagorinsky.'
  endif 
  if(k_smag_aniso > 0.0) then 
     write( stdout(),'(1x,a)')'==> NOTE: USING horz anisotropic viscosity &
          &via Smagorinsky.'
  endif      
  if(k_smag_iso==0.0)    write( stdout(),'(1x,a)')'==> NOTE: Setting horz &
       &isotropic Smagorinsky viscosity to zero.'
  if(k_smag_aniso==0.0)  write( stdout(),'(1x,a)')'==> NOTE: Setting horz &
       &anisotropic Smagorinsky viscosity to zero.'

  if(.not. viscosity_ncar) then 
    write( stdout(),'(1x,a)')'==> NOTE: NOT using NCAR scheme for computing viscosities. '
    if(vel_micom_iso > 0.0) then 
       write( stdout(),'(1x,a)')'==> NOTE: USING background horz isotropic viscosity via MICOM.'
    endif
    if(vel_micom_aniso > 0.0) then 
       write( stdout(),'(1x,a)')'==> NOTE: USING background horz anisotropic viscosity via MICOM.'
    endif 
    if(vel_micom_iso==0.0)  write( stdout(),'(1x,a)')  '==> NOTE: USING zero &
         &background horz isotropic viscosity.'
    if(vel_micom_aniso==0.0)  write( stdout(),'(1x,a)')'==> NOTE: USING zero &
         &background horz anisotropic viscosity.' 
  endif  
  if(eq_lat_micom > 0) write( stdout(),'(1x,a)')'==> NOTE: USING different &
       &background horz viscosity within equatorial zone.'

  if(k_smag_aniso > 2.0*k_smag_iso) then
    write( stdout(),'(1x,a)')'==> ERROR: Smag horz anisotropic visc too large. &
         &Must be less than twice Smagorinsky isotropic visc.'
  endif  
  if(vel_micom_aniso > 2.0*vel_micom_iso) then
    write( stdout(),'(1x,a)')'==> ERROR: Background anisotropic visc too large.&
         & Must be less than twice back isotropic visc.'
  endif  

  if(equatorial_zonal) then
    write( stdout(),'(/a)')'If using anisotropic friction, zonally orient the friction within a latitudinal band'
    write( stdout(),'(a,f12.5,a)')'of width ',equatorial_zonal_lat,' degrees.'
  endif 

  if(restrict_polar_visc) then
    write( stdout(),'(/a,f12.5)')'Using restrict_polar_visc to lower visc_crit poleward of (deg)',restrict_polar_visc_lat
    write( stdout(),'(a,f12.5)')'by an amount given by the fraction ',restrict_polar_visc_ratio
    write( stdout(),'(a/)')     'This approach useful when coupling to ice, where effective (ocn+ice) visc > ocn visc.'
  endif 

  Grd => Grid
  Dom => Domain
  nx  =  Grd%ni
  ny  =  Grd%nj

#ifndef STATIC_MEMORY
  
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso(isd:ied,jsd:jed,nk)) 
  allocate (aaniso(isd:ied,jsd:jed,nk))
  allocate (daur_dxur(isd:ied,jsd:jed))
  allocate (daur_dyur(isd:ied,jsd:jed))
  allocate (volqc(isd:ied,jsd:jed,nk))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (neptune_velocity(isd:ied,jsd:jed,2))
  allocate (stress(isd:ied,jsd:jed,0:1,0:1,2))  !stress tensor: last index 1=stress_xx, 2=stress_xy
  allocate (tmpfdelx(isd:ied,jsd:jed,2))    
  allocate (tmpfdely(isd:ied,jsd:jed,2))   
  allocate (tmp(isd:ied,jsd:jed))
  allocate (cos2theta(isc:iec,jsc:jec))
  allocate (sin2theta(isc:iec,jsc:jec))
  
#endif

  tmp(:,:) = 0.0
  
  ! Some commonly used grid arrays 
  daur_dxur(:,:) = 0.0
  daur_dyur(:,:) = 0.0
  daur_dxur(:,:) = Grd%daur(:,:)*Grd%dxur(:,:)
  daur_dyur(:,:) = Grd%daur(:,:)*Grd%dyur(:,:)

  ! critical value of viscosity, above which 2-dim linear stability is not satisfied
  visc_crit(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      dxdy = 1.0/( 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j)) )
      visc_crit(i,j) = 0.5*dxdy/(dtime + epsln)
    enddo
  enddo
  
  ! often need to use a more conservative restriction in high latitudes when coupling to ice models 
  if(restrict_polar_visc) then 
      do j=jsc,jec
         do i=isc,iec
            if(abs(Grd%yu(i,j)) > restrict_polar_visc_lat) then    
                visc_crit(i,j) = restrict_polar_visc_ratio*visc_crit(i,j)
            endif
         enddo
      enddo
  endif

  id_crit = register_static_field ('ocean_model', 'visc_crit', Grd%vel_axes_uv(1:2), 'critical viscosity', 'm^2/sec',&
                                 missing_value=-10.0, range=(/0.0,1.e20/))
  if (id_crit > 0) used = send_data (id_crit, visc_crit(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))


  ! compute some smag fields and parameters 

  ! f_smag = (delta s * k_smag/pi)^2 where (delta s) is for velocity cell
  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = (k_smag_iso/pi)**2
  coeff_aniso      = (k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = coeff_iso  *((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**2
  fsmag_aniso(:,:) = coeff_aniso*((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**2

  ! remove smag within equatorial_zonal_lat band 
  if(equatorial_no_smag) then
    do j=jsd,jed
      do i=isd,ied
        if(abs(Grd%yu(i,j)) < equatorial_zonal_lat) then
          fsmag_iso(i,j)   = 0.0
          fsmag_aniso(i,j) = 0.0
        endif 
      enddo
    enddo 
  endif 


  ! compute time independent background viscosities 
  aiso_back(:,:,:)   = 0.0
  aaniso_back(:,:,:) = 0.0

  ! grid scale dependent "MICOM" background viscosities
  if(.not. viscosity_ncar)  then
      do j=jsd,jed
         do i=isd,ied
            aiso_back(i,j,:)     = vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            aaniso_back(i,j,:)   = vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            if(abs(Grd%yu(i,j)) < eq_lat_micom) then
                aiso_back(i,j,:)   = eq_vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
                aaniso_back(i,j,:) = eq_vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            endif
         enddo
      enddo
  endif

  ! (x,y,z) viscosities according to NCAR CCSM2.0 algorithm
  ! if use ncar only for equatorial region, revert to micom for off-equatorial 
  if(viscosity_ncar)  then
     call anisotropic_ncar
     if(ncar_only_equatorial) then    
        do j=jsd,jed
           do i=isd,ied  
              if(abs(Grd%yu(i,j)) > equatorial_zonal_lat) then
                 aiso_back(i,j,:)   = vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
                 aaniso_back(i,j,:) = vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
              endif
           enddo
        enddo
     endif  
  endif 

  ! ensure background viscosities are not too large 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(aiso_back(i,j,k)   > visc_crit(i,j))  aiso_back(i,j,k)   = visc_crit(i,j)
           if(aaniso_back(i,j,k) > visc_crit(i,j))  aaniso_back(i,j,k) = visc_crit(i,j)
        enddo
     enddo
  enddo

  ! ensure viscosities maintain proper relative values so friction dissipates 
  do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            if(aaniso_back(i,j,k) >= 2.0*aiso_back(i,j,k) .and. aaniso_back(i,j,k) > 0.0)  then 
              write(stdout(),'(a,i4,a,i4,a,i3,a)')'Violating lap iso/aniso constraint at (',i,',',j,',',k,')'
              aaniso_back(i,j,k) = 1.9*aiso_back(i,j,k)
            endif 
         enddo
      enddo
  enddo

  call mpp_update_domains (aiso_back(:,:,:),  Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:),  Dom%domain2d)    

  neptune_velocity=0.0
  if(neptune) then 
     call compute_neptune_velocity(Time)
  endif 

  ! for clocks 
  id_lap_friction = mpp_clock_id('(Ocean general lap frict) ' ,grain=CLOCK_MODULE)

  ! send static fields to diag_manager 
  id_along_back   = register_static_field('ocean_model', 'along_lap_back', Grd%vel_axes_uv(1:3), &
              'U-cell background along-stream visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_across_back = register_static_field('ocean_model', 'across_lap_back', Grd%vel_axes_uv(1:3), &
              'U-cell background cross-stream visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))

  if (id_along_back > 0) then 
    used = send_data (id_along_back, aiso_back(isc:iec,jsc:jec,:)+0.5*aaniso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 
  if (id_across_back > 0) then 
    used = send_data (id_across_back, aiso_back(isc:iec,jsc:jec,:)-0.5*aaniso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  ! other viscosities for diagnostic output 

  aiso(:,:,:)   = 0.0
  aaniso(:,:,:) = 0.0
  id_aiso   = register_diag_field ('ocean_model', 'aiso_lap', Grd%vel_axes_uv(1:3), Time%model_time, &
              'U-cell isotropic visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_aaniso = register_diag_field ('ocean_model', 'aaniso_lap', Grd%vel_axes_uv(1:3), Time%model_time, &
              'U-cell anisotropic visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_along  = register_diag_field ('ocean_model', 'along', Grd%vel_axes_uv(1:3),  Time%model_time, &
              'U-cell along-stream visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_across = register_diag_field ('ocean_model', 'across', Grd%vel_axes_uv(1:3),  Time%model_time, &
              'U-cell cross-stream visc', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_lap_fric_u = register_diag_field ('ocean_model', 'lap_fric_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz lap frict on u-zonal', 'm^2/s^2', missing_value=-10.0, range=(/-1.e20,1.e20/))
  id_lap_fric_v = register_diag_field ('ocean_model', 'lap_fric_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz lap frict on v-merid','m^2/s^2', missing_value=-10.0, range=(/-1.e20,1.e20/))


end subroutine ocean_lap_friction_init
! </SUBROUTINE>  NAME="ocean_lap_friction_init"


!#######################################################################
! <FUNCTION NAME="horz_lap_friction">
!
! <DESCRIPTION>
! This function computes thickness weighted time tendency for horizontal 
! velocity (i.e., thickness weighted acceleration) from horizontal 
! Laplacian friction.  
!
! The algorithm is derived from a functional approach that ensures kinetic 
! energy is consistenty dissipated for all flow configurations. 
! The triad do-loops are expanded in order to enhance the 
! ability of cache-based machines to keep most of the variables 
! on-cache.  
! 
! Fundamental to the scheme are the rates of horizontal deformation  <BR/> 
! horizontal tension = DT = (dy)(u/dy)_x - (dx)(v/dx)_y              <BR/>    
! horizontal strain  = DS = (dx)(u/dx)_y + (dy)(v/dy)_x              <BR/> 
! Units of the tension and strain are sec^-1.
!
! Four tensions and four strains are computed for each velocity point, <BR/>
! corresponding to the four triads surrounding the point.              <BR/>  
! The following notation is used to distinguish the triads:            <BR/>  
! (0,1)=northwest triad  (1,1)=northeast triad,                        <BR/> 
! (0,0)=southwest triad, (1,0)=southeast triad
!
! A triad contributes when at least one of its velocities is            
! not a land point.  In order to obtain the correct tension          
! and strain next to boundaries, tension and strain should not be   
! masked with umask. 
!
! </DESCRIPTION>
!
function horz_lap_friction(Time, Thickness, Velocity, diag_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  logical, intent(in), optional          :: diag_flag 
  logical                                :: send_diagnostics

  real, dimension(0:1,0:1)              :: aiso_smag
  real, dimension(0:1,0:1)              :: aaniso_smag
  real, dimension(isc:iec,jsc:jec,nk,2) :: horz_lap_friction

  integer :: i, j, k, n, ip, jq, taum1, tau

  real :: u1_m10, u2_m10
  real :: u1_10,  u2_10
  real :: u1_00,  u2_00
  real :: u1_01,  u2_01
  real :: u1_0m1, u2_0m1 
  real :: dxuer_ip0, dxuer_ip1
  real :: dyunr_jq0, dyunr_jq1
  real :: usqrd, vsqrd, umagr
  real :: tension, strain, deform, delta 
  real :: tension_metric, strain_metric

  if(.not. lap_friction_on) then 
    horz_lap_friction=0.0
    return 
  endif 

  call mpp_clock_begin(id_lap_friction)

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_friction): module must be initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  stress(:,:,:,:,:) = 0.0

  taum1 = Time%taum1
  tau   = Time%tau

  volqc = 0.0
  do k=1,nk
    volqc(:,:,k) = 0.25*Grd%dau(:,:)*Thickness%dhu(:,:,k,tau)
  enddo 
  
  do k=1,nk

     do j=jsc,jec
        do i=isc,iec

           tension_metric = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh2dx(i,j) &
                            +(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh1dy(i,j)  

           strain_metric  = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh1dy(i,j) &
                            -(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh2dx(i,j)  

           if(equatorial_zonal .and. abs(Grd%yu(i,j)) <= equatorial_zonal_lat) then  
             sin2theta(i,j) = 0.0
             cos2theta(i,j) = 1.0
           else 
             usqrd = (Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))**2
             vsqrd = (Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))**2
!
! XW bug fixed based on email from Stephen Griffies
!
! modules where bug appears:
! ocean_param/mixing/horz/velocity/lap/general/ocean_lap_friction.F90
! ocean_param/mixing/horz/velocity/bih/general/ocean_bih_friction.F90
!
!bug (same in both modules):
!  umagr = 1.0/(epsln + sqrt(usqrd + vsqrd))
!
!bugfix (same in both modules):
! umagr = 1.0/(epsln + usqrd + vsqrd)
!            umagr = 1.0/(epsln + sqrt(usqrd + vsqrd))
             umagr = 1.0/(epsln + usqrd + vsqrd)
             sin2theta(i,j) = 2.0*Velocity%u(i,j,k,1,taum1)*Velocity%u(i,j,k,2,taum1)*umagr
             cos2theta(i,j) = (usqrd-vsqrd)*umagr
           endif  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are quicker. 

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
                           u2_m10 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 
           ip=1  ; jq=0  ; u1_10  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
                           u2_10  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)
           ip=0  ; jq=0  ; u1_00  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
                           u2_00  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)
           ip=0  ; jq=1  ; u1_01  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)  
                           u2_01  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 
           ip=0  ; jq=-1 ; u1_0m1 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
                           u2_0m1 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)

           ! viscosities for diagnostics 
           aiso(i,j,k)   = 0.25*(aiso_smag(0,0)   + aiso_smag(1,0)   + aiso_smag(0,1)   + aiso_smag(1,1))
           aaniso(i,j,k) = 0.25*(aaniso_smag(0,0) + aaniso_smag(1,0) + aaniso_smag(0,1) + aaniso_smag(1,1))

        enddo
     enddo

     call mpp_update_domains (stress(:,:,:,:,:), Dom%domain2d)    

     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains

     tmpfdelx(:,:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx(i,j,1) = (stress(i,j,1,0,1)   + stress(i,j,1,1,1)  )*volqc(i,j,k)    &
                +(stress(i+1,j,0,0,1) + stress(i+1,j,0,1,1))*volqc(i+1,j,k)
           tmpfdelx(i,j,2) = (stress(i,j,1,0,2)   + stress(i,j,1,1,2)  )*volqc(i,j,k)    &  
                +(stress(i+1,j,0,0,2) + stress(i+1,j,0,1,2))*volqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely(:,:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely(i,j,1) = (stress(i,j,0,1,2)   + stress(i,j,1,1,2)  )*volqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,2) + stress(i,j+1,1,0,2))*volqc(i,j+1,k) 
           tmpfdely(i,j,2) = (stress(i,j,0,1,1)   + stress(i,j,1,1,1)  )*volqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,1) + stress(i,j+1,1,0,1))*volqc(i,j+1,k)
        enddo
     enddo


     ! compute acceleration (m/sec^2) arising from horizontal friction
     do n=1,2

        tmp(:,:) =  BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n))    
        horz_lap_friction(isc:iec,jsc:jec,k,n) = tmp(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,k)

        ! reduce to 5-point laplacian at bottom to avoid problems with thin bottom partial cells
        if (bottom_5point) then
            tmp(:,:) = BDX_EU(aiso_back(:,:,k)*FMX(Thickness%dhu(:,:,k,tau))&
                      *FDX_PU(Velocity%u(:,:,k,n,taum1)-neptune_velocity(:,:,n))) &
                      +BDY_NU(aiso_back(:,:,k)*FMY(Thickness%dhu(:,:,k,tau)) &
                      *FDY_PU(Velocity%u(:,:,k,n,taum1)-neptune_velocity(:,:,n)))
            do j=jsc,jec
               do i=isc,iec
                  if(k==Grd%kmu(i,j)) horz_lap_friction(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)
               enddo
            enddo
        endif

     enddo

  enddo  !end of k-loop

  ! send viscosity and friction for diagnostic output
  if(send_diagnostics) then

      if (id_aiso > 0)    used = send_data (id_aiso, aiso(:,:,:),   &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      if (id_aaniso > 0)  used = send_data (id_aaniso, aaniso(:,:,:), &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      if (id_along > 0)   used = send_data (id_along,  aiso(:,:,:)+0.5*aaniso(:,:,:),  &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      if (id_across > 0)  used = send_data (id_across, aiso(:,:,:)-0.5*aaniso(:,:,:),  &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_lap_fric_u > 0) used = send_data(id_lap_fric_u, horz_lap_friction(isc:iec,jsc:jec,:,1), &
                                    Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))

      if (id_lap_fric_v > 0) used = send_data(id_lap_fric_v, horz_lap_friction(isc:iec,jsc:jec,:,2), &
                                    Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif

  call mpp_clock_end(id_lap_friction)

end function horz_lap_friction
! </FUNCTION> NAME="horz_lap_friction"


!#######################################################################
! <FUNCTION NAME="BDX_EU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in X of a quantity defined on the east 
! face of a U-cell. Slightly modified version of BDX_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i-1/2,j) BDX_EU_smag changes dimensions by m^-3 
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the east face of a U-cell
! </IN>
!
function BDX_EU_smag(a)

  real, dimension(isd:ied,jsd:jed) :: a, BDX_EU_smag
  integer :: i, j

  do j=jsd,jed
    do i=isd+1,ied
      BDX_EU_smag(i,j) = (Grd%dyue_dxuer(i,j)*a(i,j) - Grd%dyue_dxuer(i-1,j)*a(i-1,j))*daur_dyur(i,j)
    enddo
    BDX_EU_smag(isd,j) = 0.0
  enddo

end function BDX_EU_smag
! </FUNCTION> NAME="BDX_EU_smag"


!#######################################################################
! <FUNCTION NAME="BDY_NU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in Y of a quantity defined on the north
! face of a U-cell. Slightly modified version of BDY_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i,j-1/2) BDY_EU_smag changes dimensions by m^-3 
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the north face of a U-cell
! </IN>
!
function BDY_NU_smag(a)

  real, dimension(isd:ied,jsd:jed) :: a, BDY_NU_smag
  integer :: i, j

  do j=jsd+1,jed
    do i=isd,ied            
      BDY_NU_smag(i,j) = (Grd%dxun_dyunr(i,j)*a(i,j) - Grd%dxun_dyunr(i,j-1)*a(i,j-1))*daur_dxur(i,j)
    enddo
  enddo
  BDY_NU_smag(:,jsd) = 0.0

end function BDY_NU_smag
! </FUNCTION> NAME="BDY_EU_smag"


!#######################################################################
! <SUBROUTINE NAME="anisotropic_ncar">
!
! <DESCRIPTION>
!
!     Spatially-varying anisotropic viscosity initialization
!
!     This routine defines NCOM-like spatial distributions of
!     viscosity coefficients F_PARA and F_PERP.
!     Uses NCAR CCSM2.0 algorithm with cm^2/sec --> m^2/sec.  
!
!     written by:     Stephen Yeager 3/2000                           <BR/> 
!     modified by:    Gokhan Danabasoglu (08/2001)                    <BR/>
!     port to mom4:   Stephen.Griffies@noaa.gov (9/2002)  
!
!   "A_viscosity" = F_PARA = Along = viscosity parallel to flow 
!                  = max{0.5*visc_vel_scale(z)*A*max[dx,dy],vconst_6}
!
!   where                                                                <BR/>    
!          A = 0.425 * cos(pi*y*radian/30) + 0.575   for |y*radian| < 30 <BR/>
!          A = 0.15                                  otherwise 
!
!   Here, A provides a horizontal variation for visc_vel_scale.
!
!   "B_viscosity" = F_PERP = Across = viscosity perpendicular to flow = max( bu, bv)
!
!   and                                                                  <BR/>    
!        F_PARA = min(F_PARA, AMAX_CFL),                                 <BR/> 
!        F_PERP = min(F_PERP, AMAX_CFL),                                 <BR/>   
!        F_PARA = max(F_PARA, F_PERP)                                    <BR/> 
!   are enforced 
!
!   In the above equations, 
!
!        bu  = vconst_1 * ( 1 + vconst_2  * ( 1 + cos( 2*y + pi ) ) )        <BR/>
!        bv  = vconst_3 * beta_f * dx^3   * exp( - (vconst_4 * distance)^2 ) <BR/>
!
!   with                                                                     <BR/>        
!        beta_f         (x,y)   = 2 * omega * cos(ULAT(i,j)) / radius        <BR/>  
!        distance       (x,y,z) = actual distance to "vconst_5" points       <BR/>  
!                                 west of the nearest western boundary       <BR/>  
!        dx             (x,y)   = DXU(i,j)                                   <BR/>  
!        dy             (x,y)   = DYU(i,j)                                   <BR/> 
!        visc_vel_scale (z)     = vconst_7 * exp(-zt(k)/visc_vel_scale_length)  <BR/> 
!        visc_vel_scale_length  = e-folding scale ( = 1500.0e2 cm)           <BR/> 
!        y              (x,y)   = ULAT(i,j), latitude of "u/v" grid pts in radians   <BR/> 
!        In mom4, ULAT(radians) = xu*pi/180 with xu(i,j) the longitude of U grid points in degrees
!
!   "vconst_#" are input parameters defined in namelist ocean_lap_friction_general_nml. 
!   "vconst_1", "vconst_6", and "vconst_4" have dimensions of cm^2/s,
!   cm^2/s, and 1/cm, respectively. "vconst_5" is an INTEGER.
!
!   NOTE: The nearest western boundary computations are done along the
!         model i-grid lines. Therefore, viscosity based on these are 
!         only approximate in the high Northern Hemisphere when using 
!         generalized coordinates with coordinate pole(s) shifted onto land. 
!
! </DESCRIPTION>
!
subroutine anisotropic_ncar

  integer :: i, j, k
  integer :: n, ncount
  integer :: is, ie, ip1, ii
  integer :: index, indexo
  integer, dimension(nx) :: iwp
  integer, dimension(nx) :: nwbp_global

  real, dimension(nx,ny)            :: dxtn_global
  real, dimension(nx,ny)            :: kmu_global
  real, dimension(isd:ied,jsd:jed)  :: kmu_tmp
  real, dimension(nx)               :: dist_g 

  real :: dist_max, vvsl  
  real :: a_para, bu, bv 
  real :: f_para, f_perp
  real :: visc_vel_scale

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (anisotropic_ncar): module must be initialized')
  endif 

  dist_max         = 1.e10    ! distance for ACC region (cm)
  vvsl             = 1500.e2  ! exponential decay scale (cm) for decay away from western boundaries

  dxtn_global(:,:) = 0.0
  kmu_global(:,:)  = 0.0
  kmu_tmp(:,:)     = Grd%kmu(:,:)

  call mpp_global_field(Dom%domain2d,kmu_tmp,kmu_global)  
  call mpp_global_field(Dom%domain2d,Grd%dxtn,dxtn_global)

  do k=1,nk

! 100.0 factor in exponential is to convert mom4 zt from (m) to (cm)
! visc_vel_scale is in cm/sec.  It is used for along-stream viscosity. 
     visc_vel_scale = vconst_7 * exp(-100.0*Grd%zt(k)/vvsl)   

     do j=jsc,jec


!determine nearest western boundary
        ncount         = 0
        nwbp_global(:) = 0
        iwp(:)         = 0
        dist_g(:)      = 0.0

        do i=1,Grd%ni
           ip1 = i+1
           if(i==Grd%ni) then 
               if(Grd%cyclic) then 
                   ip1=1
               else
                   ip1=i  
               endif
           endif
           if(kmu_global(i,j+Dom%joff)<k .and. kmu_global(ip1,j+Dom%joff) >= k) then 
               ncount      = ncount+1
               iwp(ncount) = i
           endif
        enddo

        if (ncount > 0) then
            do n=1,ncount-1
               is = iwp(n)
               ie = iwp(n+1)-1
               do i=is,ie
                  nwbp_global(i)=is
               enddo
            enddo
            do i=1,Grd%ni
               if(nwbp_global(i)==0) nwbp_global(i) = iwp(ncount)
            enddo
        endif


!determine distance (cm) to nearest western boundary
        do i=1,Grd%ni

           index  = nwbp_global(i)
           indexo = index + vconst_5
           if ( index .eq. 0) then
               dist_g(i) = dist_max
           elseif ( i .ge. index  .and.  i .le. indexo ) then
               dist_g(i) = 0.0
           elseif ( (i .gt. indexo) ) then
               dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
           elseif ( i .lt. index ) then
               if (indexo .le. Grd%ni) then
                   if (i .eq. 1) then
                       dist_g(i) = 0.0
                       do ii=indexo+1,Grd%ni
                          dist_g(i)=dxtn_global(ii,j+Dom%joff) + dist_g(i)
                       enddo
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i)
                   else
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
                   endif
               else
                   if (i .le. (indexo - Grd%ni)) then
                       dist_g(i) = 0.0
                   else
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
                   endif
               endif
           endif

        enddo

! convert distance in mom4 (m) to POP1.4 (cm)
        dist_g(:) = 100.0*dist_g(:)


!calculate viscosity over comp domain 
        do i=isc,iec

! f_para is viscosity parallel to flow 
! 100.0 factor in f_para converts (m) to (cm)
           a_para = 0.15
           if ( abs(Grd%yu(i,j)) < 30.0) then  ! yu is in degrees, not radians 
               a_para = 0.425*cos(pi*Grd%yu(i,j)/30.0) + 0.575
           endif
           f_para = max(0.5 * visc_vel_scale * a_para * 100.0*max(Grd%dxu(i,j),Grd%dyu(i,j)), vconst_6 )

! f_perp is viscosity perpendicular to flow 
! 1e4 factor in bv converts (m) to (cm) for dxu and beta_f 
           bu = vconst_1*(1.0 + vconst_2*(1.0 + cos((2.0*Grd%yu(i,j)/radian)+pi)))
           bv = 1e4*vconst_3*Grd%beta(i,j)*Grd%dxu(i,j)**3
           bv = bv*exp(-(vconst_4*dist_g(i+Dom%ioff))**2)
           f_perp = min(max(bu,bv),0.99*f_para)

           if(debug_ncar_A) then 
             f_perp=f_para
           endif 
           if(debug_ncar_B) then 
             f_para=f_perp
           endif 

! revert to isotropic background viscosity outside equatorial_zonal_lat
           if(ncar_isotropic_off_equator .and. abs(Grd%yu(i,j)) >= equatorial_zonal_lat) then 
             f_perp=f_para
           endif   

! 1.e-4 converts from cm^2/sec to m^2/sec 
           aiso_back(i,j,k)   = 1.e-4*0.5*(f_para + f_perp) 
           aaniso_back(i,j,k) = 1.e-4*(f_para - f_perp)

        enddo

     enddo ! j
  enddo    ! k

end subroutine anisotropic_ncar
! </SUBROUTINE>  NAME="anisotropic_ncar"


!#######################################################################
! <SUBROUTINE NAME="horz_lap_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the Laplacian 
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
subroutine horz_lap_viscosity_check

  integer       :: num, i, j, k
  integer, save :: unit=6, max_num=3

  if(.not. lap_friction_on) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_viscosity_check): module must be initialized')
  endif 

  write (stdout(),'(/60x,a/)') ' Excessive horizontal Laplacian friction summary:'
  write (stdout(),'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then 
          if (aiso(i,j,k) > visc_crit(i,j) .and. num < max_num) then
            num = num + 1
            write (unit,9600) 'aiso(',aiso(i,j,k), visc_crit(i,j), i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum(num)
  if (num > max_num) write (stdout(),*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es10.3,' m^2/s) exceeds max value (',es10.3,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')
9601  format(/' Warning: ',a,es10.3,' m^2/s) exceeds max value (',es10.3,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')

end subroutine horz_lap_viscosity_check
! </SUBROUTINE> NAME="horz_lap_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="horz_lap_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the LLaplacian grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
! <IN NAME="u" TYPE="real" DIM="(isd:ied,jsd:jed,nk,2)">
! Horizontal velocity field at time tau
! </IN>
subroutine horz_lap_reynolds_check(Time, Velocity)

  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  if(.not. lap_friction_on) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_reynolds_check): module must be initialized')
  endif 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso(i,j,k) + epsln)
        ramn = 1.0/(aiso(i,j,k) + epsln)
        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j))*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j))*ramn
        if (reyy > reyny) then
          ireyny = i
          jreyny = j
          kreyny = k
          reyny  = reyy
          reynv  = Velocity%u(i,j,k,2,tau)
          reynmv = 1.0/ramn
        endif
      enddo
    enddo
  enddo
  write (stdout(),'(/60x,a/)') ' Horizontal Laplacian Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0 .and. reynx > 1.e-6) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, Grd%xu(ireynx,jreynx), Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0 .and. reyny > 1.e-6) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, Grd%xu(ireyny,jreyny), Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)

end subroutine horz_lap_reynolds_check
! </SUBROUTINE> NAME="horz_lap_reynolds_check"


!#######################################################################
! <SUBROUTINE NAME="compute_neptune_velocity">
!
! <DESCRIPTION>
! Compute Neptune velocity.  
!
! Method follows that used in MOM2/3 as implemented by 
! Greg Holloway (zounds@ios.bc.ca) and Michael Eby (eby@uvic.ca) 
!
! Neptune is calculated as an equilibrium streamfunction given by 
! pnep = -f*snep*snep*ht and is applied through friction 
!
! ht    = depth of tracer cells 
! snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
!
! Neptune length scale snep has a value of senep at the
! equator and smoothly changes to spnep at the poles
!
! Reference:
! Holloway, G., 1992: Representing topographic stress for large
! scale ocean models, J. Phys. Oceanogr., 22, 1033-1046
!
! </DESCRIPTION>
!
subroutine compute_neptune_velocity(Time)

  type(ocean_time_type), intent(in) :: Time
  real, dimension(isd:ied,jsd:jed)  :: neptune_ft, neptune_psi
  real                              :: neptune_length
  integer                           :: i,j

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (compute_neptune_velocity): module must be initialized')
  endif 

  do j=jsd,jed
    do i=isd,ied
      neptune_ft(i,j)  = 2.0*omega*sin(Grd%phit(i,j))
      neptune_length   = neptune_length_pole + 0.5*(neptune_length_eq-neptune_length_pole)*(1.0 + cos(2.0*Grd%phiu(i,j)))
      neptune_psi(i,j) = -neptune_ft(i,j)*neptune_length*Grd%ht(i,j)
    enddo
  enddo
  neptune_velocity(:,:,1) = -FDY_ET(FAX(neptune_psi(:,:)))  
  neptune_velocity(:,:,2) =  FDX_NT(FAY(neptune_psi(:,:)))  
  call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                          Dom%domain2d,gridtype=BGRID_NE)


  id_neptune_u   = register_static_field('ocean_model', 'neptune_u', Grd%vel_axes_uv(1:2), &
              'Zonal velocity from Neptune', 'm/sec', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_neptune_v   = register_static_field('ocean_model', 'neptune_v', Grd%vel_axes_uv(1:2), &
              'Meridional velocity from Neptune', 'm/sec', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_neptune_psi   = register_static_field('ocean_model', 'neptune_psi', Grd%tracer_axes(1:2), &
              'Streamfunction for Neptune', 'm^2/sec', missing_value=-1.e12, range=(/-1.e12,1.e12/))
  id_neptune_ft   = register_static_field('ocean_model', 'neptune_ft', Grd%tracer_axes(1:2), &
              'Coriolis parameter used for Neptune Psi', '1/sec', missing_value=-1.e3, range=(/-1.e3,1.e3/))

  if (id_neptune_ft  > 0) then 
    used = send_data (id_neptune_ft, neptune_ft(isc:iec,jsc:jec), &
                      Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))
  endif 
  if (id_neptune_psi  > 0) then 
    used = send_data (id_neptune_psi, neptune_psi(isc:iec,jsc:jec), &
                      Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))
  endif 
  if (id_neptune_u  > 0) then 
    used = send_data (id_neptune_u, neptune_velocity(isc:iec,jsc:jec,1), &
                      Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))
  endif 
  if (id_neptune_v  > 0) then 
    used = send_data (id_neptune_v, neptune_velocity(isc:iec,jsc:jec,2), &
                      Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))
  endif 


end subroutine compute_neptune_velocity
! </SUBROUTINE> NAME="compute_neptune_velocity"



end module ocean_lap_friction_mod
      
      




