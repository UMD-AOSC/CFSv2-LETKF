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
module ocean_bih_friction_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes thickness weighted time tendency for horizontal 
! velocity arising from horizontal biharmonic friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted time tendency for horizontal 
! velocity arising from horizontal biharmonic friction. 
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
!<NAMELIST NAME="ocean_bih_friction_general_nml">
!  <DATA NAME="bih_friction_on" TYPE="logical">
!  Must be true to use this module. 
!  </DATA> 
!  <DATA NAME="k_smag_iso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="k_smag_aniso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky anisotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_iso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_aniso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal" TYPE="real">
!  Orient the anisotropic friction within a latitudinal band according to zonal direction. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal_lat" TYPE="real">
!  Latitudinal band to use the zonal friction orientation. 
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
!  <DATA NAME="bottom_5point" TYPE="logical">
!  To alleviate problems with small partial cells, it is often necessary to reduce the 
!  operator to the traditional 5-point Laplacian at the ocean bottom.  This logical 
!  implements this mixing. 
!  <DATA NAME="vel_micom_bottom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity for 5point Laplacian
!  at the bottom. 
!  </DATA> 

!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use mpp_mod,             only: mpp_error, mpp_max, mpp_sum, mpp_pe

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: BAY, BAX, BDX_EU, BDY_NU, FDX_PU, FDY_PU, FMX, FMY
use ocean_types_mod,     only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_velocity_type

implicit none

private

public ocean_bih_friction_init
public horz_bih_friction
public horz_bih_viscosity_check
public horz_bih_reynolds_check

private BDX_EU_smag
private BDY_NU_smag

real :: k_smag_iso          = 0.0    ! smag scaling coeff for isotropic viscosity (dimensionless) 
real :: k_smag_aniso        = 0.0    ! smag scaling coeff for anisotripic viscosity (dimensionless) 
real :: vel_micom_iso       = 0.0    ! background scaling velocity for isotropic viscosity (m/sec) 
real :: vel_micom_aniso     = 0.0    ! background scaling velocity for anisotropic viscosity (m/sec) 
real :: eq_vel_micom_iso    = 0.2    ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real :: eq_vel_micom_aniso  = 0.0    ! background scaling velocity (m/sec) within equatorial band for anisotropic visc
real :: eq_lat_micom        = 0.0    ! equatorial latitude band for micom (degrees)
real :: vel_micom_bottom    = 0.01   ! velocity scale for determining viscosity at the bottom 
real :: equatorial_zonal_lat= 0.0    ! latitudinal band for orienting the friction along zonal direction
logical :: equatorial_zonal = .false.! zonally orient anisotropic friction w/i equatorial_zonal_lat
logical :: bottom_5point    =.true.  ! for bottom Laplacian 5point mixing to avoid problems with thin partial cells. 

! for diagnostics 
integer :: id_aiso=-1
integer :: id_aaniso=-1
integer :: id_along=-1
integer :: id_across=-1
integer :: id_along_back=-1
integer :: id_across_back=-1
integer :: id_crit=-1
integer :: id_bih_fric_u=-1
integer :: id_bih_fric_v=-1
logical :: used

! clock ids
integer :: id_bih_friction

! time step 
real ::  dtime = 0.0  ! time step used for friction tendency (2*dtuv if threelevel, dtuv if twolevel) 

#include <ocean_memory.h>


#ifdef STATIC_MEMORY

real, dimension(isd:ied,jsd:jed)    :: fsmag_iso   ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: fsmag_aniso ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: aiso_back   ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed)    :: aaniso_back ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk) :: aiso        ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: aaniso      ! anisotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed)    :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed)    :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed,nk) :: volqc       ! (1/4)*dxu*dyu*dhu (m^3) for time-independent quarter-cell volume
real, dimension(isd:ied,jsd:jed)    :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(isd:ied,jsd:jed)    :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(isd:ied,jsd:jed,0:1,0:1,2) :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(isd:ied,jsd:jed,2)         :: tmplap, tmpfdelx, tmpfdely
real, dimension(isd:ied,jsd:jed)           :: tmp
real, dimension(isc:iec,jsc:jec)           :: cos2theta
real, dimension(isc:iec,jsc:jec)           :: sin2theta

#else

real, dimension(:,:), allocatable   :: fsmag_iso   ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(:,:), allocatable   :: fsmag_aniso ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(:,:), allocatable   :: aiso_back   ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(:,:), allocatable   :: aaniso_back ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(:,:,:), allocatable :: aiso        ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(:,:,:), allocatable :: aaniso      ! anisotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(:,:), allocatable   :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(:,:), allocatable   :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(:,:,:), allocatable :: volqc       ! (1/4)*dxu*dyu*dhu (m^3) for time-independent quarter-cell volume
real, dimension(:,:), allocatable   :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(:,:), allocatable   :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(:,:,:,:,:), allocatable        :: stress             !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(:,:,:), allocatable            :: tmplap, tmpfdelx, tmpfdely
real, dimension(:,:), allocatable              :: tmp
real, dimension(:,:), allocatable              :: cos2theta
real, dimension(:,:), allocatable              :: sin2theta
#endif

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: /general/ocean_bih_friction.f90 ($Id$)'
character (len=128) :: tagname = &
     '$Name$'

logical :: bih_friction_on       = .false.
logical :: module_is_initialized = .FALSE.

namelist /ocean_bih_friction_general_nml/ bih_friction_on, k_smag_iso, k_smag_aniso, vel_micom_iso, &
                                          vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom, &
                                          vel_micom_bottom, bottom_5point, equatorial_zonal, equatorial_zonal_lat 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_friction_init">

! <DESCRIPTION>
! Initialize the horizontal biharmonic friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bih_friction_init(Grid, Domain, Time, d_time)

  type(ocean_grid_type), target, intent(in)   :: Grid
  type(ocean_domain_type), target, intent(in) :: Domain
  type(ocean_time_type), intent(in)           :: Time
  real, intent(in)                            :: d_time

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real    :: coeff_iso, coeff_aniso, dsmin

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bih_friction_mod(ocean_bih_friction_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_bih_friction_general_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_bih_friction_general_nml)  
  write (stdlog(),ocean_bih_friction_general_nml)
  ierr = check_nml_error(io_status,'ocean_bih_friction_general_nml')
  call close_file(ioun)

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(bih_friction_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_bih_friction_mod from bih/general.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_bih_friction_mod.')
    return
  endif 

  Grd => Grid
  Dom => Domain

  dtime = d_time

  write(stdout(),*) ' ' 

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_bih_friction_mod: using forward time step of (secs)', dtime 

  if (Dom%xhalo > 1 .or. Dom%yhalo > 1) then
    write (stdout(),'(/1x,a)') ' ==> NOTE: General biharmonic friction allows xhalo=yhalo=1.'
  endif

  if(bottom_5point) then
    write(stdout(),'(/)')  
    write(stdout(),'(/1x,a)') ' ==> NOTE: Will make horizontal friction to a 5point Laplacian on the bottom'
    write(stdout(),'(a)') '        This helps alleviate numerical problems with thin bottom partial cells.'
  endif 

  if(k_smag_iso > 0.0) then 
     write( stdout(),*)' Computing horz isotropic bih viscosity via Smagorinsky.'
  endif
  if(k_smag_iso==0.0)  write( stdout(),*)' Setting horz isotropic bih Smagorinsky viscosity to zero.'

  if(k_smag_aniso > 0.0) then 
     write( stdout(),*)' Computing horz anisotropic bih viscosity via Smagorinsky.'
     write( stdout(),'(1x,a)')' ==>WARNING: anisotropic biharmonic Smagorinsky not well tested.'
  endif
  if(k_smag_aniso==0.0)  write( stdout(),*)' Setting horz anisotropic bih Smagorinsky viscosity to zero.'

  if(vel_micom_iso > 0.0) then 
     write( stdout(),*)' Computing background horz bih isotropic viscosity via MICOM.'
  endif
  if(vel_micom_iso==0.0)  write( stdout(),*)' Setting the background horz bih isotropic viscosity to zero.'

  if(vel_micom_aniso > 0.0) then
     write( stdout(),*)' Computing background horz anisotropic bih viscosity via MICOM.'
     write( stdout(),'(1x,a)')' ==> WARNING: anisotropic biharmonic background viscosity not well tested.'
  endif 
  if(vel_micom_aniso==0.0)  write( stdout(),*)' Setting the background horz bih anisotropic viscosity to zero.' 

  if(k_smag_aniso > 2.0*k_smag_iso) then
     call mpp_error(FATAL,'==>Error: Smag horz bih anisotropic visc too large. Must be < 2*isotropic-visc.')
  endif  
  if(vel_micom_aniso > 2.0*vel_micom_iso) then
    call mpp_error(FATAL,'==>Error: Background bih horz anisotropic visc too large. Must be < 2.0*backgound isotropic visc.')
  endif  
  if(eq_lat_micom > 0) write( stdout(),*)'Setting different background horz bih viscosity within equatorial zone.'

  if(equatorial_zonal) then
    write( stdout(),*)'If using anisotropic biharmonic friction, zonally orient friction within latitudinal'
    write( stdout(),*)'band of width ',equatorial_zonal_lat,' degrees.'
  endif 


#ifndef STATIC_MEMORY
  
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (aiso_back(isd:ied,jsd:jed))
  allocate (aaniso_back(isd:ied,jsd:jed))
  allocate (aiso(isd:ied,jsd:jed,nk)) 
  allocate (aaniso(isd:ied,jsd:jed,nk))
  allocate (daur_dxur(isd:ied,jsd:jed))
  allocate (daur_dyur(isd:ied,jsd:jed))
  allocate (volqc(isd:ied,jsd:jed,nk))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (visc_bottom(isd:ied,jsd:jed))
  allocate (stress(isd:ied,jsd:jed,0:1,0:1,2))  !stress tensor: last index 1=stress_xx, 2=stress_xy
  allocate (tmplap(isd:ied,jsd:jed,2))    
  allocate (tmpfdelx(isd:ied,jsd:jed,2))    
  allocate (tmpfdely(isd:ied,jsd:jed,2))   
  allocate (tmp(isd:ied,jsd:jed))
  allocate (cos2theta(isc:iec,jsc:jec))
  allocate (sin2theta(isc:iec,jsc:jec))
  
#endif

  tmp(:,:) = 0.0
  tmplap(:,:,:) = 0.0
  
  ! Some commonly used grid arrays 
  daur_dxur(:,:) = 0.0
  daur_dyur(:,:) = 0.0
  daur_dxur(:,:) = Grd%daur(:,:)*Grd%dxur(:,:)
  daur_dyur(:,:) = Grd%daur(:,:)*Grd%dyur(:,:)

  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = 0.125*(k_smag_iso/pi)**2
  coeff_aniso      = 0.125*(k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = Grd%umask(:,:,1)*coeff_iso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4
  fsmag_aniso(:,:) = Grd%umask(:,:,1)*coeff_aniso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4

  aiso_back(:,:)   = 0.0
  aaniso_back(:,:) = 0.0
  visc_bottom(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied 
      visc_bottom(i,j)   = Grd%umask(i,j,1)*vel_micom_bottom* &
                            ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))
      aiso_back(i,j)     = Grd%umask(i,j,1)*vel_micom_iso* &
                            ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      aaniso_back(i,j)   = Grd%umask(i,j,1)*vel_micom_aniso* &
                            ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      if(abs(Grd%yu(i,j)) < eq_lat_micom) then
        aiso_back(i,j)   = Grd%umask(i,j,1)*eq_vel_micom_iso* &
                            ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
        aaniso_back(i,j) = Grd%umask(i,j,1)*eq_vel_micom_aniso* &
                            ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      endif
    enddo  
  enddo

  ! critical value of viscosity, above which 2-dim linear stability is not satisfied.
  visc_crit(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      dsmin = min(Grd%dxu(i,j),Grd%dxu(i+1,j),Grd%dyu(i,j),Grd%dyu(i,j+1))
      visc_crit(i,j) =  0.0625*dsmin**4/(dtime+epsln)
    enddo
  enddo

  id_crit = register_static_field ('ocean_model', 'visc_crit', Grd%vel_axes_uv(1:2), 'critical viscosity', 'm^4/sec',&
                                 missing_value=-10.0, range=(/0.0,1.e20/))
  if (id_crit > 0) used = send_data (id_crit, visc_crit(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))


  ! ensure that background viscosities are not too large 
  do j=jsc,jec
     do i=isc,iec
        if(aiso_back(i,j)   > visc_crit(i,j))  aiso_back(i,j)   = visc_crit(i,j)
        if(aaniso_back(i,j) > visc_crit(i,j))  aaniso_back(i,j) = visc_crit(i,j)
     enddo
  enddo

  ! ensure viscosities maintain proper relative values so that friction dissipates 
  do j=jsc,jec
     do i=isc,iec
        if(aaniso_back(i,j) >= 2.0*aiso_back(i,j) .and. aaniso_back(i,j) > 0.0)  then 
          write(stdout(),'(a,i4,a,i4,a,i3,a)')'Violating bih iso/aniso constraint at (',i,',',j,')'
          aaniso_back(i,j) = 1.9*aiso_back(i,j)
        endif 
     enddo
  enddo

  call mpp_update_domains (aiso_back(:,:),  Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:),  Dom%domain2d)    

  ! for clocks 
  id_bih_friction  = mpp_clock_id('(Ocean general bih friction) ' ,grain=CLOCK_MODULE)

  ! viscosity for diagnostic output 
  aiso(:,:,:)    = 0.0
  aaniso(:,:,:)  = 0.0
  id_aiso    = register_diag_field ('ocean_model', 'aiso_bih', Grd%vel_axes_uv(1:3), Time%model_time,&
     'U-cell isotropic bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e20/))
  id_aaniso  = register_diag_field ('ocean_model', 'aaniso_bih', Grd%vel_axes_uv(1:3), Time%model_time,&
     'U-cell  anisotropic bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e20/))
  id_along  = register_diag_field ('ocean_model', 'along_bih', Grd%vel_axes_uv(1:3), Time%model_time,&
     'U-cell along-stream bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e20/))
  id_across  = register_diag_field ('ocean_model', 'across_bih', Grd%vel_axes_uv(1:3), Time%model_time, &
     'U-cell cross-stream bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e20/))
  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz bih frict on u-zonal', 'm^2/s^2', missing_value=-10.0, range=(/-1.e20,1.e20/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz bih frict on v-merid','m^2/s^2', missing_value=-10.0, range=(/-1.e20,1.e20/))

  id_along_back   = register_static_field('ocean_model', 'along_bih_back', Grd%vel_axes_uv(1:2), &
              'U-cell background along-stream bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_across_back = register_static_field('ocean_model', 'across_bih_back', Grd%vel_axes_uv(1:2), &
              'U-cell background cross-stream bih visc', 'm^4/sec', missing_value=-10.0, range=(/-10.0,1.e10/))

  if (id_along_back > 0) then 
    used = send_data (id_along_back, aiso_back(isc:iec,jsc:jec)+0.5*aaniso_back(isc:iec,jsc:jec), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))
  endif 

  if (id_across_back > 0) then 
    used = send_data (id_across_back, aiso_back(isc:iec,jsc:jec)-0.5*aaniso_back(isc:iec,jsc:jec), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))
  endif 

end subroutine ocean_bih_friction_init
! </SUBROUTINE>  NAME="ocean_bih_friction_init"


!#######################################################################
! <FUNCTION NAME="horz_bih_friction">
!
! <DESCRIPTION>
! This function computes the time tendency for horizontal 
! velocity (i.e., the acceleration) from horizontal biharmonic friction.  
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
! As shown in Griffies and Hallberg (2000), 
! a biharmonic operator with a nonconstant viscosity is guaranteed to 
! dissipate kinetic energy *only* when using the sqrt of the biharmonic
! viscosity at each of the two stages of the algorithm. 
! The sqrt approach is employed here.  
!
! </DESCRIPTION>
!
function horz_bih_friction(Time, Thickness, Velocity, diag_flag)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_thickness_type), intent(in)   :: Thickness
  type(ocean_velocity_type), intent(in)    :: Velocity
  logical, intent(in), optional            :: diag_flag 
  logical                                  :: send_diagnostics
  real, dimension(isc:iec,jsc:jec,0:1,0:1) :: aiso_smag
  real, dimension(isc:iec,jsc:jec,0:1,0:1) :: aaniso_smag
  real, dimension(isc:iec,jsc:jec,nk,2)    :: horz_bih_friction

  integer :: i, j, k, n, ip, jq, tau, taum1

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

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bih_friction_mod (horz_bih_friction): module needs to be initialized')
  endif 

  if(.not. bih_friction_on) then 
    horz_bih_friction=0.0
    return 
  endif 

  call mpp_clock_begin(id_bih_friction)

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


     ! Laplacian part of algorithm
     do j=jsc,jec
        do i=isc,iec

           tension_metric = -Velocity%u(i,j,k,1,taum1)*Grd%dh2dx(i,j) + Velocity%u(i,j,k,2,taum1)*Grd%dh1dy(i,j)  
           strain_metric  = -Velocity%u(i,j,k,1,taum1)*Grd%dh1dy(i,j) - Velocity%u(i,j,k,2,taum1)*Grd%dh2dx(i,j)  

           if(equatorial_zonal .and. abs(Grd%yu(i,j)) <= equatorial_zonal_lat) then  
             sin2theta(i,j) = 0.0
             cos2theta(i,j) = 1.0
           else 
             usqrd = Velocity%u(i,j,k,1,taum1)*Velocity%u(i,j,k,1,taum1)
             vsqrd = Velocity%u(i,j,k,2,taum1)*Velocity%u(i,j,k,2,taum1)
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
          ! dimensions[aiso_smag]=m^2/s^0.5
          ! dimensions[tension and strain]=1/s
          ! dimensions[stress]=m^2/s^1.5

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_m10 = Velocity%u(i+ip,j+jq,k,2,taum1)
           ip=1  ; jq=0  ; u1_10  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_10  = Velocity%u(i+ip,j+jq,k,2,taum1) 
           ip=0  ; jq=0  ; u1_00  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_00  = Velocity%u(i+ip,j+jq,k,2,taum1) 
           ip=0  ; jq=1  ; u1_01  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_01  = Velocity%u(i+ip,j+jq,k,2,taum1)
           ip=0  ; jq=-1 ; u1_0m1 = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_0m1 = Velocity%u(i+ip,j+jq,k,2,taum1)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ! viscosities for diagnostics 
           aiso(i,j,k)   = 0.25*(  aiso_smag(i,j,0,0)**2   + aiso_smag(i,j,1,0)**2   &
                                 + aiso_smag(i,j,0,1)**2   + aiso_smag(i,j,1,1)**2)
           aaniso(i,j,k) = 0.25*(  aaniso_smag(i,j,0,0)**2 + aaniso_smag(i,j,1,0)**2 &
                                 + aaniso_smag(i,j,0,1)**2 + aaniso_smag(i,j,1,1)**2)

        enddo
     enddo

     call mpp_update_domains (stress(:,:,:,:,:), Dom%domain2d)    

     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=m^5/s^1.5

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

     ! compute laplacian operator
     ! dimensions[tmplap]=m/s^1.5
     do n=1,2
        tmp(:,:)                  = BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n))    
        tmplap(isc:iec,jsc:jec,n) = tmp(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,k)/Thickness%dhu(isc:iec,jsc:jec,k,tau)
     enddo
     call mpp_update_domains (tmplap(:,:,:), Dom%domain2d)    

     ! Second part of the iteration
     ! tmplap[m/s^1.5] replaces velocity[m/s]
     ! dimensions[aiso_smag]=m^2/s^0.5
     ! dimensions[tension and strain]=1/s^1.5
     ! dimensions[stress]=m^2/s^2

     stress(:,:,:,:,:)   = 0.0

     do j=jsc,jec
        do i=isc,iec

           tension_metric = -tmplap(i,j,1)*Grd%dh2dx(i,j) + tmplap(i,j,2)*Grd%dh1dy(i,j)  
           strain_metric  = -tmplap(i,j,1)*Grd%dh1dy(i,j) - tmplap(i,j,2)*Grd%dh2dx(i,j)  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are quicker. 

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = tmplap(i+ip,j+jq,1) ; u2_m10 = tmplap(i+ip,j+jq,2)
           ip=1  ; jq=0  ; u1_10  = tmplap(i+ip,j+jq,1) ; u2_10  = tmplap(i+ip,j+jq,2) 
           ip=0  ; jq=0  ; u1_00  = tmplap(i+ip,j+jq,1) ; u2_00  = tmplap(i+ip,j+jq,2) 
           ip=0  ; jq=1  ; u1_01  = tmplap(i+ip,j+jq,1) ; u2_01  = tmplap(i+ip,j+jq,2)
           ip=0  ; jq=-1 ; u1_0m1 = tmplap(i+ip,j+jq,1) ; u2_0m1 = tmplap(i+ip,j+jq,2)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

        enddo
     enddo

     call mpp_update_domains (stress(:,:,:,:,:), Dom%domain2d)    

     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=m^5/s^2

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

    ! compute acceleration from horizontal biharmonic friction
    ! dimensions[tmp]=m^2/s^2
    ! dimensions[horz_bih_friction]=m/s^2
     do n=1,2

        tmp(:,:)                               =  BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n)) 
        horz_bih_friction(isc:iec,jsc:jec,k,n) = -tmp(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,k)

        ! reduce to 5-point laplacian at bottom to avoid problems with thin bottom partial cells
        if (bottom_5point) then
            tmp(:,:) = BDX_EU(visc_bottom(:,:)*FMX(Thickness%dhu(:,:,k,tau))*FDX_PU(Velocity%u(:,:,k,n,taum1))) &
                 +BDY_NU(visc_bottom(:,:)*FMY(Thickness%dhu(:,:,k,tau))*FDY_PU(Velocity%u(:,:,k,n,taum1)))
            do j=jsc,jec
               do i=isc,iec
                  if(k==Grd%kmu(i,j)) horz_bih_friction(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)
               enddo
            enddo
        endif

     enddo !end of n-loop

  enddo  !end of k-loop

  ! send viscosity and friction for diagnostic output
  if(send_diagnostics) then

      if (id_aiso > 0)    used = send_data (id_aiso, aiso(:,:,:), &
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

      if (id_bih_fric_u > 0) used = send_data(id_bih_fric_u, horz_bih_friction(isc:iec,jsc:jec,:,1), &
                                    Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
      if (id_bih_fric_v > 0) used = send_data(id_bih_fric_v, horz_bih_friction(isc:iec,jsc:jec,:,2), &
                                    Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif

  call mpp_clock_end(id_bih_friction)

end function horz_bih_friction
! </FUNCTION> NAME="horz_bih_friction"


!#######################################################################
! <FUNCTION NAME="BDX_EU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in X of a quantity defined on the east 
! face of a U-cell. Slightly modified version of BDX_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i-1/2,j).
!
! BDX_EU_smag changes dimensions by m^-3 
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
! at (i,j-1/2).
!
! BDY_EU_smag changes dimensions by m^-3 
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
! <SUBROUTINE NAME="horz_bih_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the biharmonic
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine horz_bih_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: dsmin, crit

  if(.not. bih_friction_on) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bih_friction_mod (horz_bih_viscosity_check): module needs to be initialized')
  endif 

  write (stdout(),'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdout(),'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then            
          if (aiso(i,j,k) > visc_crit(i,j) .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'aiso(',aiso(i,j,k), visc_crit(i,j), i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdout(),*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es14.6,' m^4/s) exceeds max value (',es14.6,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')

end subroutine horz_bih_viscosity_check
! </SUBROUTINE> NAME="horz_bih_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="horz_bih_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the biharmonic grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine horz_bih_reynolds_check(Time, Velocity)

  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bih_friction_mod (horz_bih_reynolds_check): module needs to be initialized')
  endif 

  if(.not. bih_friction_on) return 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso(i,j,k) + epsln)
        ramn = 1.0/(aiso(i,j,k) + epsln)

        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j)**3)*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j)**3)*ramn
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
  write (stdout(),'(/60x,a/)') ' Horizontal biharmonic Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, Grd%xu(ireynx,jreynx), Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, Grd%xu(ireyny,jreyny), Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine horz_bih_reynolds_check
! </SUBROUTINE> NAME="horz_bih_reynolds_check"


end module ocean_bih_friction_mod
      
      




