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
module ocean_sigma_diffuse_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted time tendency for tracer from sigma diffusion 
! model grid cells.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for 
! tracer arising from Laplacian diffusion within the bottom-most 
! model grid cells. The diffusivity used to determine the strength 
! of the tendency is generally set to be a function of the local 
! horizontal grid spacing.  Diffusivity is the sum of an a priori 
! background plus a velocity dependent diffusivity.  It is large 
! if there is a heavier parcel living adjacent within the 
! "sigma layer" above a lighter parcel.  It is small otherwise. 
! 
! The thickness of the "sigma layer" is time independent and 
! equated to the thickness of the bottom partial cell thickness dht.
!
! When the logical sigma_tracer_diffusion=.true., other diffusive 
! fluxes passing across bottom cell boundaries are masked to zero.
! Generally, sigma_tracer_diffusion=.false. has been found useful. 
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R. Doscher and A. Beckmann
! Effects of a bottom boundary layer parameterization 
! in a coarse-resolution model of the North Atlantic Ocean
! Journal of Atmospheric and Oceanic Technology (1999), vol 17 pages 698--707
! </REFERENCE>
!
! <REFERENCE>
! A. Beckmann and R. Doscher",
! A method for improved representation of dense water spreading over 
! topography in geopotential--coordinate models
! Journal of Physical Oceanography (1997) vol 27, pages 581--59
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, and A. Rosati,
! A Guide to MOM4 for Users and Developers (2002)
! </REFERENCE>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_sigma_diffuse_nml">
!  <DATA NAME="sigma_diffuse_on" TYPE="logical">
!  Must be true to use this module. 
!  </DATA> 
!  <DATA NAME="tmask_sigma_on" TYPE="logical">
!  IF .true. then masks out fluxes passing into the sigma layer, except those 
!  associated with sigma diffusion. Typically set to .false.  
!  </DATA> 
!  <DATA NAME="sigma_diffusivity" UNITS="m^2/sec" TYPE="real">
!  Sigma tracer diffusivity for use if not using micom diffusivity.   
!  </DATA> 
!  <DATA NAME="sigma_diffusivity_ratio" UNITS="dimensionless" TYPE="real">
!  When flow along sigma surface is stable (i.e., heavy parcels are below lighter parcels)
!  then sigma diffusivity is reduced by sigma_diffusivity_ratio from the case where 
!  heavy parcels are above lighter parcels.  
!  </DATA>
!  <DATA NAME="sigma_thickness_max" UNITS="meter" TYPE="real">
!  The maximum thickness of the bottom sigma diffusion layer.  
!  If dht > sigma_thickness_max, then bbl_thickness is set to sigma_thickness_max.
!  Otherwise, bbl_thickness=dht.  
!  </DATA>   
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the sigma diffusivity is set according to a velocity scale 
!  times the grid spacing. 
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!</NAMELIST>

use constants_mod,    only: rho0, grav, c2dbars, pi
use diag_manager_mod, only: register_static_field, register_diag_field, send_data
use fms_mod,          only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,          only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE
use mpp_mod,          only: mpp_error

use ocean_domains_mod,   only: get_local_indices, set_ocean_domain
use ocean_density_mod,   only: density
use ocean_operators_mod, only: FMX, FMY, FDX_PT_flat, FDY_PT_flat, BDX_ET, BDY_NT
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_adv_vel_type, ocean_thickness_type

implicit none

public ocean_sigma_diffuse_init
public sigma_diffusion

private

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_sigma_diff   ! thickness weighted time tendency for tracer from sigma diffusion
integer, dimension(:), allocatable  :: id_xflux_sigma  ! for i-directed heat flux from sigma diffusion
integer, dimension(:), allocatable  :: id_yflux_sigma  ! for j-directed heat flux from sigma diffusion
integer :: id_gradX_ztp=-1
integer :: id_gradY_ztp=-1 
integer :: id_tmask_sigma=-1
integer :: id_bbl_thickness=-1
integer :: id_diff_cet=-1
integer :: id_diff_cnt=-1

#include <ocean_memory.h>

#ifdef STATIC_MEMORY

real, public, dimension(isd:ied,jsd:jed) :: tmask_sigma  ! Mask defining which bottom cells are active sigma cells 
real, dimension(isd:ied,jsd:jed)  :: bbl_thickness       ! Spatially dependent bbl thickness (m)
real, dimension(isd:ied,jsd:jed)  :: gradX_ztp           ! X-derivative of the T-cell bottom topography (dimensionless)
real, dimension(isd:ied,jsd:jed)  :: gradY_ztp           ! Y-derivative of the T-cell bottom topography (dimensionless)
real, dimension(isd:ied,jsd:jed)  :: diff_cet            ! Sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(isd:ied,jsd:jed)  :: diff_cnt            ! Sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(isd:ied,jsd:jed)  :: diff_max            ! Maximum a priori diffusivity (m2/sec)
real, dimension(isd:ied,jsd:jed)  :: tracer_sigma        ! Tracer concentration in the sigma layer 

#else

real, public, dimension(:,:),   allocatable :: tmask_sigma    ! Mask defining which bottom cells are active sigma cells 
real, dimension(:,:),           allocatable :: bbl_thickness  ! Spatially dependent bbl thickness (m)
real, dimension(:,:),           allocatable :: gradX_ztp      ! X-derivative of the T-cell bottom topography (dimensionless)
real, dimension(:,:),           allocatable :: gradY_ztp      ! Y-derivative of the T-cell bottom topography (dimensionless)
real, dimension(:,:),           allocatable :: diff_cet       ! Sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:),           allocatable :: diff_cnt       ! Sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:),           allocatable :: diff_max       ! Maximum a priori diffusivity (m2/sec)
real, dimension(:,:),           allocatable :: tracer_sigma   ! Tracer concentration in the sigma layer 

#endif


type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers = -1
integer :: index_temp=-1
integer :: index_salt=-1
real    :: press_constant    ! A useful constant

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

integer :: unit=6  
logical :: module_is_initialized = .FALSE.


! nml settings 
logical, public  :: tmask_sigma_on  = .false. ! will mask out non-sigma diffusive fluxes passing into bottom level
logical :: sigma_diffuse_on         = .false. ! must be true to have any sigma diffusion occur 
logical :: verbose_init             = .true.  ! for verbose initialization printout
real     :: sigma_diffusivity       = 0.0e3   ! sigma tracer diffusivity (m^2/sec)
real     :: sigma_diffusivity_ratio = 1.0e-6  ! ratio of min to max sigma diffusivities
real     :: sigma_thickness_max     = 50.0    ! max thickness (m) of bottom sigma diffusion layer
real     :: vel_micom               = 0.50    ! constant velocity scale (m/s) for setting micom diffusivity
logical  :: tracer_mix_micom        =.false.  ! if true, diffusivity made a function of horz grid spacing

namelist /ocean_sigma_diffuse_nml/ sigma_diffuse_on, tmask_sigma_on, sigma_diffusivity,&
                                   sigma_diffusivity_ratio, tracer_mix_micom, vel_micom,&
                                   sigma_thickness_max, verbose_init

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_sigma_diffuse_init">
!
! <DESCRIPTION>
! Initialize the sigma diffusion module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_sigma_diffuse_init(Grid, Domain, Time, Thickness, T_prog, dtime)

  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, intent(in)                            :: dtime

  integer :: ioun, io_status, ierr
  integer :: i, j, k, n, kip1, kjp1, num
  real    :: dxdymn, asigma_crit
  
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_sigma_diffuse_mod (ocean_sigma_diffuse_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_sigma_diffuse_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_sigma_diffuse_nml)  
  write (stdlog(),ocean_sigma_diffuse_nml)
  ierr = check_nml_error(io_status,'ocean_sigma_diffuse_nml')
  call close_file(ioun)

  write (stdout(),'(/a,e12.5/)') '==> Note from ocean_sigma_diffuse_mod: max thickness of sigma layer is ', sigma_thickness_max
  write(stdout(),'(/a,f10.2)')'==> Note from ocean_sigma_diffuse_mod: using forward time step of (secs)', dtime 
  
  num_prog_tracers = size(T_prog(:))

  press_constant = rho0*grav*c2dbars

#ifndef STATIC_MEMORY  
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(sigma_diffuse_on) then 
    call mpp_error(NOTE, '==>Note from ocean_sigma_diffuse_mod: USING ocean_sigma_diffuse_mod.')
  else 
    call mpp_error(NOTE, '==>Note from ocean_sigma_diffuse_mod:s NOT using ocean_sigma_diffuse_mod.')
    return
  endif 

  Dom => Domain
  Grd => Grid

  call set_ocean_domain(Dom_flux, Grid, xhalo=Dom%xhalo, yhalo=Dom%yhalo, name='flux dom sigma')

  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,'==>Error in ocean_sigma_diffuse_mod: temp and/or salt not identified in tracer array')
  endif 

#ifndef STATIC_MEMORY
  allocate (tracer_sigma(isd:ied,jsd:jed))
  allocate (diff_cet(isd:ied,jsd:jed))
  allocate (diff_cnt(isd:ied,jsd:jed))
  allocate (diff_max(isd:ied,jsd:jed))
  allocate (gradX_ztp(isd:ied,jsd:jed))
  allocate (gradY_ztp(isd:ied,jsd:jed))
  allocate (tmask_sigma(isd:ied,jsd:jed))
  allocate (bbl_thickness(isd:ied,jsd:jed))
#endif

  ! Tracer concentration in the sigma layer 

  tracer_sigma(:,:) = 0.0

  ! Sigma diffusivities 

  diff_cet(:,:) = 0.0
  diff_cnt(:,:) = 0.0
  diff_max(:,:) = sigma_diffusivity 

  ! Micom diffusivity
  if(tracer_mix_micom .and. sigma_diffuse_on) then
    do j=jsd,jed
      do i=isd,ied
        diff_max(i,j) = vel_micom*(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)/(Grd%dxt(i,j)+Grd%dyt(i,j)))
      enddo
    enddo
    if(verbose_init) then 
      do j=jsc,jec
        write (stdout(),'(a,i4,a,e14.7,a)') ' Laplacian diffusivity in sigma layer at (isc,',j,') = ',diff_max(isc,j),' m^2/s'
      enddo 
    endif 
  endif

  ! Gradient of bottom topography

  gradX_ztp=0.0 ;  gradY_ztp=0.0
  do j=jsd,jed-1
      do i=isd,ied-1 
        kip1 = max(1,Grd%kmt(i+1,j))
        kjp1 = max(1,Grd%kmt(i,j+1))
        k    = max(1,Grd%kmt(i,j))
        gradX_ztp(i,j) = Grd%tmask(i+1,j,1)*(Thickness%ztp(i+1,j,kip1)-Thickness%ztp(i,j,k))*Grd%dxter(i,j)
        gradY_ztp(i,j) = Grd%tmask(i,j+1,1)*(Thickness%ztp(i,j+1,kjp1)-Thickness%ztp(i,j,k))*Grd%dytnr(i,j)
      enddo
  enddo

  ! Set sigma mask.  If wish to mask out regions, do so here.
  ! Otherwise, all bottom cells in regions of more than 2 vertical cells  
  ! are treated as sigma in the tracer diffusion equation.  

  tmask_sigma(:,:) = Grd%tmask(:,:,1)
  do j=jsd,jed
    do i=isd,ied
      if(Grd%kmt(i,j) < 3) tmask_sigma(i,j) = 0.0
    enddo
  enddo

  ! Set sigma layer thickness equal bottom z-level thickness

  bbl_thickness(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied
      k = Grd%kmt(i,j)
      if(tmask_sigma(i,j)==1.0) bbl_thickness(i,j) = min(sigma_thickness_max,Thickness%dht(i,j,k,1))
    enddo
  enddo

  ! checks

  if (dtime /= 0.0 .and. sigma_diffuse_on) then
    num = 0
    do j=jsc,jec
      do i=isc,iec
        dxdymn = 2.0/(1.0/(Grd%dxt(i,j))**2 + 1.0/Grd%dyt(i,j)**2)
        asigma_crit = 0.5*dxdymn/dtime
        if (diff_max(i,j)*Grd%tmask(i,j,1)  > asigma_crit .and. num <= 10) then
          num = num + 1
          if (num == 1) write (unit,'(/,(1x,a))')&
          '==> Warning: Diffusive criteria exceeded for "asigma". use a smaller "dtts", or "asigma".  Show first 10 violations:'
          write (stdout(),'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') &
          ' at (i,j)= (',i,',',j,'), (lon,lat)= (', Grd%xt(i,j),',',Grd%yt(i,j),&
          '),  "ah" = ', diff_max(i,j),' m^2/s. the critical value =',asigma_crit,' m^2/s'
        endif
      enddo
    enddo
  endif

  ! static diagnostics 

  id_bbl_thickness = register_static_field ('ocean_model', 'bbl_thick', Grd%tracer_axes(1:2), &
                    'Thickness of BBL', 'm', missing_value=-1.e3, range=(/-1.e3,1.e3/))
  id_gradX_ztp  = register_static_field ('ocean_model', 'gradX_ztp', Grd%tracer_axes(1:2), &
                  'X-derivative of bottom depth', 'm/m', missing_value=-1.e3, range=(/-1.e3,1.e3/))
  id_gradY_ztp  = register_static_field ('ocean_model', 'gradY_ztp', Grd%tracer_axes(1:2), &
                  'Y-derivative of bottom depth', 'm/m', missing_value=-1.e3, range=(/-1.e3,1.e3/))
  id_tmask_sigma = register_static_field ('ocean_model', 'tmask_sigma', Grd%tracer_axes(1:2), &
                  'Mask for bottom sigma layer', 'm', missing_value=-1.e3, range=(/-1.e3,1.e5/))

  if (id_bbl_thickness  > 0) then 
    used = send_data (id_bbl_thickness, bbl_thickness(isc:iec,jsc:jec), Time%model_time)
  endif
  if (id_gradX_ztp  > 0) then 
    used = send_data (id_gradX_ztp, gradX_ztp(isc:iec,jsc:jec), Time%model_time)
  endif
  if (id_gradY_ztp  > 0) then 
    used = send_data (id_gradY_ztp, gradY_ztp(isc:iec,jsc:jec), Time%model_time)
  endif 
  if (id_tmask_sigma  > 0) then 
    used = send_data (id_tmask_sigma, tmask_sigma(isc:iec,jsc:jec), Time%model_time)
  endif 


  ! dynamic diagnostics 

  id_diff_cet = register_diag_field ('ocean_model', 'diff_cet_sigma', Grd%tracer_axes(1:2), &
                Time%model_time, 'X-sigma diffusivity', 'm^2/s', missing_value=-10.0, &
                range=(/-10.0,1.e8/))
  id_diff_cnt = register_diag_field ('ocean_model', 'diff_cnt_sigma', Grd%tracer_axes(1:2), &
                 Time%model_time, 'Y-sigma diffusivity', 'm^2/s', missing_value=-10.0, &
                 range=(/-10.0,1.e8/))

  allocate (id_sigma_diff(num_prog_tracers))
  id_sigma_diff = -1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') then
        id_sigma_diff(n) = register_diag_field ('ocean_model', 'sigma_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:2), Time%model_time, &
                               'thk wghtd sigma-diffusion heating ', 'Watts/m^2', &
                               missing_value=-1.e10, range=(/-1.e10,1.e10/))
     else 
        id_sigma_diff(n) = register_diag_field ('ocean_model', 'sigma_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:2), Time%model_time, &
                              'thk wghtd sigma-physics on '//trim(T_prog(n)%name), 'm*kg/sec', &
                               missing_value=-1.e10, range=(/-1.e10,1.e10/))
     endif 
  enddo

  ! Fluxes are nonzero only for bottom cell. 
  ! So the thickness weighted fluxes are equivalent to the vertically integrated fluxes
  allocate (id_xflux_sigma(num_prog_tracers))
  allocate (id_yflux_sigma(num_prog_tracers))
  id_xflux_sigma = -1
  id_yflux_sigma = -1

  do n=1,num_prog_tracers
    if(n == index_temp) then 
      id_xflux_sigma(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_sigma', &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time, 'rho_cp*sigma_xflux*dyt*dht*temp', &
                   'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_yflux_sigma(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_sigma', &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time, 'rho_cp*sigma_yflux*dxt*dht*temp', &
                   'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
  else
      id_xflux_sigma(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_sigma', &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time, &
                   'rho0*sigma_xflux*dyt*dht*tracer for'//trim(T_prog(n)%name),&
                   trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_yflux_sigma(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_sigma', &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time, &
                   'rho0*sigma_yflux*dxt*dht*tracer for'//trim(T_prog(n)%name),&
                   trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
    endif 
  enddo


end subroutine ocean_sigma_diffuse_init
! </SUBROUTINE> NAME="ocean_sigma_diffuse_init"


!#######################################################################
! <SUBROUTINE NAME="sigma_diffusion">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted time tendency for
! tracer from sigma diffusion and stores result in trace th_tendency. 
! </DESCRIPTION>
!
subroutine sigma_diffusion (Time, Thickness, T_prog, Adv_vel)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_adv_vel_type), intent(in)        :: Adv_vel

  real, dimension(isd:ied,jsd:jed)  :: fx, fy, tchg, tmp_flux
  real, dimension(2)                :: density_test
  real :: conversion_factor, check, press, temp0, salt0, tempip1, saltip1, tempjp1, saltjp1
  integer :: i, j, k, kip1, kjp1, taum1, tau, n

  if (.not. sigma_diffuse_on) return

  taum1 = Time%taum1
  tau   = Time%tau

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_sigma_diffuse_mod (sigma_diffusion): size mismatch for tracer array')
  endif 
  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_sigma_diffuse_mod (sigma_diffusion): module must be initialized')
  endif 

  ! Determine sigma diffusivities.
  ! Allow enhanced diffusion even if topography is flat (check=0).
  ! Use deepest depth to approximate hydrostatic pressure for computing density. 

  diff_cet(:,:)     = 0.0
  diff_cnt(:,:)     = 0.0

  do j=jsd,jed-1
     do i=isd,ied-1

        if(tmask_sigma(i,j)==1.0) then 

           k       = Grd%kmt(i,j)
           temp0   = T_prog(index_temp)%field(i,j,k,taum1)
           salt0   = T_prog(index_salt)%field(i,j,k,taum1)
           kip1    = max(1,Grd%kmt(i+1,j))
           tempip1 = T_prog(index_temp)%field(i+1,j,kip1,taum1)
           saltip1 = T_prog(index_salt)%field(i+1,j,kip1,taum1)
           kjp1    = max(1,Grd%kmt(i,j+1))
           tempjp1 = T_prog(index_temp)%field(i,j+1,kjp1,taum1)
           saltjp1 = T_prog(index_salt)%field(i,j+1,kjp1,taum1)
           press = press_constant*max(Thickness%ztp(i,j,k),Thickness%ztp(i+1,j,kip1))
           density_test(1) = density(salt0,temp0,press)
           density_test(2) = density(saltip1,tempip1,press)
           check = (density_test(2)-density_test(1))*gradX_ztp(i,j)        
           if(check <= 0.0) then 
              diff_cet(i,j)     = diff_max(i,j)              
           else 
              diff_cet(i,j)     = diff_max(i,j)*sigma_diffusivity_ratio              
           endif
           if(check <= 0.0 .and. Adv_vel%uh_et(i,j,k)*gradX_ztp(i,j) >= 0.0) then 
              diff_cet(i,j)     =   min(diff_cet(i,j)  &
                                  + abs(Adv_vel%uh_et(i,j,k))*Grd%dxt(i,j)/Thickness%dht(i,j,k,tau), 2.0*diff_max(i,j))     
           endif

           press = press_constant*max(Thickness%ztp(i,j,k),Thickness%ztp(i,j+1,kjp1))
           density_test(1) = density(salt0,temp0,press)     
           density_test(2) = density(saltjp1,tempjp1,press)
           check = (density_test(2)-density_test(1))*gradY_ztp(i,j)        
           if(check <= 0.0) then 
              diff_cnt(i,j)     = diff_max(i,j)              
           else 
              diff_cnt(i,j)     = diff_max(i,j)*sigma_diffusivity_ratio               
           endif
           if(check <= 0.0 .and. Adv_vel%vh_nt(i,j,k)*gradY_ztp(i,j) >= 0.0) then 
              diff_cnt(i,j)     = min(diff_cnt(i,j)  &
                                  + abs(Adv_vel%vh_nt(i,j,k))*Grd%dyt(i,j)/Thickness%dht(i,j,k,tau), 2.0*diff_max(i,j))     
           endif
           
        endif
        
     enddo
  enddo

  ! set tracer concentration in the sigma layer
  
  do n = 1, num_prog_tracers
     tracer_sigma(:,:) = 0.0 
     do j=jsd,jed
        do i=isd,ied
           if(tmask_sigma(i,j)==1.0) then
              k=Grd%kmt(i,j)
              tracer_sigma(i,j) = T_prog(n)%field(i,j,k,taum1)     
           endif
        enddo
     enddo

     ! Compute diffusion time tendency, weighted by sigma-layer thickness
     fx(:,:) = diff_cet(:,:)*FDX_PT_flat(tracer_sigma(:,:))*FMX(bbl_thickness(:,:)*tmask_sigma(:,:))
     fy(:,:) = diff_cnt(:,:)*FDY_PT_flat(tracer_sigma(:,:))*FMY(bbl_thickness(:,:)*tmask_sigma(:,:))

     !for redundancies at Arctic fold 
     if (Grd%tripolar) call mpp_update_domains(fx(:,:), fy(:,:), Dom_flux%domain2d, gridtype=CGRID_NE)   

     tchg(:,:) = tmask_sigma(:,:)*(BDX_ET(fx(:,:)) + BDY_NT(fy(:,:)))

     do j=jsc,jec
        do i=isc,iec
           k=max(1,Grd%kmt(i,j))
           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + tchg(i,j)
        enddo
     enddo

     if (id_sigma_diff(n) > 0) used = send_data (id_sigma_diff(n), T_prog(n)%conversion*tchg(:,:), &
                                      Time%model_time, rmask=tmask_sigma(:,:), &
                                      is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


     ! Fluxes are nonzero only for bottom cell. 
     ! So the thickness weighted fluxes are equivalent to the vertically integrated fluxes
     ! minus sign accounts for mom4 sign convention for fluxes
     if(id_xflux_sigma(n) > 0) then 
        used = send_data (id_xflux_sigma(n), -1.0*T_prog(n)%conversion*Grd%dyte(:,:)*fx(:,:), &
               Time%model_time, rmask=tmask_sigma(:,:), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

     if(id_yflux_sigma(n) > 0) then 
        used = send_data (id_yflux_sigma(n), -1.0*T_prog(n)%conversion*Grd%dxtn(:,:)*fy(:,:), &
               Time%model_time, rmask=tmask_sigma(:,:), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

     if(n==index_temp) then 
       if(id_diff_cet > 0) used = send_data (id_diff_cet, diff_cet(:,:), &
                           Time%model_time, rmask=tmask_sigma(:,:), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if(id_diff_cnt > 0) used = send_data (id_diff_cnt, diff_cnt(:,:), &
                           Time%model_time, rmask=tmask_sigma(:,:), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

  enddo

end subroutine sigma_diffusion
! </SUBROUTINE> 


end module ocean_sigma_diffuse_mod
      
      





