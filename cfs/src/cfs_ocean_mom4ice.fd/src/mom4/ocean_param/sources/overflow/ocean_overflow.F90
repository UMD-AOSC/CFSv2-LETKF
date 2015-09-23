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
module ocean_overflow_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Tracer source from discharging dense shallow water into the abyss
! at the parcel's depth of neutral buoyancy.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Tracer source from discharging dense shallow water into the abyss
! at the parcel's depth of neutral buoyancy.  Follow the approach
! of Campin and Goosse (1999).
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Campin and Goosse (1999): Parameterization of density-driven downsloping flow 
! for a coarse-resolution model in z-coordinate", Tellus 51A, pages 412-430
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_overflow_nml">
!  <DATA NAME="overflow_mu" TYPE="real" UNITS="inverse seconds">
!  Dissipation rate for the bottom friction.  Campin and Goosse 
!  suggest overflow_mu=10^-4
!  </DATA> 
!  <DATA NAME="overflow_delta" TYPE="real" UNITS="dimensionless">
!  Fraction of a grid cell participating in the overflow process. 
!  Campin and Goosse suggest overflow_delta=1/3. 
!  </DATA> 
!  <DATA NAME="overflow_umax" TYPE="real" UNITS="m/s">
!  Maximum downslope speed. 
!  </DATA> 
!  <DATA NAME="debug_overflow" TYPE="logical">
!  For debugging 
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!</NAMELIST>
!

use axis_utils_mod,      only: nearest_index, frac_index
use constants_mod,       only: epsln, grav, rho0r
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use fms_mod,             only: write_version_number, error_mesg, FATAL, NOTE, mpp_error
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,             only: mpp_error, mpp_chksum

use ocean_density_mod,   only: density
use ocean_domains_mod,   only: get_local_indices, set_ocean_domain
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_density_type, ocean_thickness_type
use ocean_util_mod,      only: write_timestamp
use ocean_workspace_mod, only: wrk1, wrk2, wrk3, wrk4

implicit none

private 

public ocean_overflow_init
public overflow

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

#include <ocean_memory.h>

! for the algorithm
real    :: overflow_factor ! grav*overflow_delta/(rho0*overflow_mu)  (m^4/kg/sec)
integer :: ip(4),jq(4)     ! to label four surrounding cells 
integer :: ijsc(4),ijec(4) ! for shifting do-loop indices to minimize calls to mpp_update_domains

#ifdef STATIC_MEMORY
real, dimension(isd:ied,jsd:jed)        :: slope_x         ! topographic slope (dimensionless) in the i-direction
real, dimension(isd:ied,jsd:jed)        :: slope_y         ! topographic slope (dimensionless) in the j-direction 
real, dimension(isd:ied,jsd:jed,4)      :: topog_slope     ! topographic slope (dimensionless) for surrounding directions 
real, dimension(isd:ied,jsd:jed,4)      :: topog_step      ! designates where downslope flow is topographically possible 
real, dimension(isd:ied,jsd:jed,4)      :: flux_int_z      ! for diagnostics 
real, dimension(isd:ied,jsd:jed,nk)     :: source_overflow ! thickness weighted tendency from overflow transport 
#else 
real, dimension(:,:),   allocatable   :: slope_x         ! topographic slope (dimensionless) in the i-direction
real, dimension(:,:),   allocatable   :: slope_y         ! topographic slope (dimensionless) in the j-direction 
real, dimension(:,:,:), allocatable   :: topog_slope     ! topographic slope (dimensionless) for surrounding directions 
real, dimension(:,:,:), allocatable   :: topog_step      ! designates where downslope flow is topgraphically possible 
real, dimension(:,:,:), allocatable   :: flux_int_z      ! for diagnostics 
real, dimension(:,:,:), allocatable   :: source_overflow ! thickness weighted tendency from overflow transport 
#endif

! for diagnostic manager 
logical :: used
integer :: id_slope_x=-1
integer :: id_slope_y=-1
integer :: id_topog_step_1=-1
integer :: id_topog_step_2=-1
integer :: id_topog_step_3=-1
integer :: id_topog_step_4=-1
integer :: id_overflow_xflux=-1 
integer :: id_overflow_yflux=-1
integer, dimension(:), allocatable :: id_overflow
integer, dimension(:), allocatable :: id_overflow_xflux_int_z
integer, dimension(:), allocatable :: id_overflow_yflux_int_z


character(len=128) :: version=&
       '=>Using: ocean_overflow.f90 ($Id$)'
character (len=128) :: tagname=&
     '$Name$'

integer :: num_prog_tracers=0
integer :: unit=6       ! processor zero writes to unit 6
logical :: module_is_initialized = .FALSE.
integer :: global_sum_flag      ! flag for mpp_global_sum

! parameters set from nml
logical :: use_overflow   = .false. ! must be set .true. in nml to enable this scheme
logical :: debug_overflow = .false. ! for debugging
real    :: overflow_mu    = 1.e-4   ! frictional dissipation rate (sec^-1) at bottom  
real    :: overflow_delta = 0.3333  ! fraction of a grid cell participating in overflow 
real    :: overflow_umax  = 0.01    ! maximum downslope speed (m/s)
logical :: do_bitwise_exact_sum = .false. 

namelist /ocean_overflow_nml/ use_overflow, debug_overflow, overflow_mu, overflow_delta, &
                              overflow_umax, do_bitwise_exact_sum

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_overflow_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers into the abyss next to topography.
! </DESCRIPTION>
!
  subroutine ocean_overflow_init(Grid, Domain, Time, T_prog, debug)

    type(ocean_grid_type), target               :: Grid
    type(ocean_domain_type), target             :: Domain
    type(ocean_time_type), target               :: Time
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    logical, intent(in), optional               :: debug

    integer :: chksum
    integer :: io_status, ioun, ierr
    integer :: i,j,m,n,kbot

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error from ocean_overflow_mod (ocean_overflow_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

#ifndef STATIC_MEMORY
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    Dom => Domain
    Grd => Grid

    num_prog_tracers = size(T_prog(:))

    if (PRESENT(debug)) debug_overflow = debug

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_overflow_nml,IOSTAT=io_status)
    write (stdlog(),ocean_overflow_nml)
    write (stdout(),'(/)') 
    write (stdout(),ocean_overflow_nml)
    ierr = check_nml_error(io_status,'ocean_overflow_nml')
    call close_file (ioun)

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    if(.not. use_overflow) then 
      call mpp_error(NOTE,'==>From ocean_overflow_mod: NOT using Campin and Goosse overflow scheme.')
      return
    endif 
    if(use_overflow) then 
      call mpp_error(NOTE,'==>From ocean_overflow_mod: USING Campin and Goosse overflow scheme.')
    endif
    if(debug_overflow) then 
      call mpp_error(NOTE,'==>From ocean_overflow_mod: USING debug_overflow')
    endif 

#ifndef STATIC_MEMORY
    allocate (source_overflow(isd:ied,jsd:jed,nk))
#endif 
    source_overflow(:,:,:) = 0.0    


    ! compute a common factor that is time independent 
    overflow_factor = grav*overflow_delta*rho0r/(epsln+overflow_mu)

    ! compute topographic slope arrays for the i-slope and j-slope
    ! slopes are centered on the i-face and j-face of tracer cells 
#ifndef STATIC_MEMORY
    allocate (slope_x(isd:ied,jsd:jed))
    allocate (slope_y(isd:ied,jsd:jed))
#endif 
    slope_x = 0.0 
    slope_y = 0.0
    do j=jsc,jec
      do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)     
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
      enddo
    enddo 
    call mpp_update_domains(slope_x(:,:),Dom%domain2d)
    call mpp_update_domains(slope_y(:,:),Dom%domain2d)

    ! topographic slope for the four surrounding directions 
#ifndef STATIC_MEMORY
    allocate (topog_slope(isd:ied,jsd:jed,4))
#endif 
    topog_slope(:,:,:) = 0.0
    do j=jsc,jec
       do i=isc,iec
          m=1 ; topog_slope(i,j,m) = slope_x(i,j)  
          m=2 ; topog_slope(i,j,m) = slope_y(i,j)  
          m=3 ; topog_slope(i,j,m) = slope_x(i-1,j)
          m=4 ; topog_slope(i,j,m) = slope_y(i,j-1)
       enddo
    enddo

    ! compute directions from an (i,j) point where topography deepens.
    ! these directions may potentially have downslope flow. 
    ! insist that downslope flow occurs only when there are more kmt cells
    ! in the adjacent column. 
    ! also insist that downslope flow does not involve k=1 cells. 
#ifndef STATIC_MEMORY
    allocate (topog_step(isd:ied,jsd:jed,4))
#endif
    topog_step(:,:,:) = 0.0
    do j=jsc,jec
      do i=isc,iec
         kbot = Grd%kmt(i,j)
         if(kbot > 1) then 
           if(Grd%kmt(i+1,j) > kbot) topog_step(i,j,1)=1.0
           if(Grd%kmt(i,j+1) > kbot) topog_step(i,j,2)=1.0
           if(Grd%kmt(i-1,j) > kbot) topog_step(i,j,3)=1.0
           if(Grd%kmt(i,j-1) > kbot) topog_step(i,j,4)=1.0
         endif 
      enddo
    enddo   
    ! Block out the bipolar fold in order to ensure tracer conservation.
    ! The reason we do so is related to how the algorithm reaches between  
    ! two adjacent columns of tracer points.  When the column straddles 
    ! the bipolar fold, the code logic is not general and so it actually 
    ! attempts to reach to a non-adjacent column.  Shy of generalizing 
    ! the logic, we simply do not consider overflow for points along fold. 
    if(jec+Dom%joff==Dom%jeg) then 
       topog_step(:,jec,:) = 0.0
    endif 

    ! labels for the four tracer cells surrounding a 
    ! central tracer cell (moving counter-clockwise)
    m=1 ; ip(m)=1  ; jq(m)=0  ; ijsc(m) = 1 ; ijec(m) = 0
    m=2 ; ip(m)=0  ; jq(m)=1  ; ijsc(m) = 1 ; ijec(m) = 0
    m=3 ; ip(m)=-1 ; jq(m)=0  ; ijsc(m) = 0 ; ijec(m) = 1
    m=4 ; ip(m)=0  ; jq(m)=-1 ; ijsc(m) = 0 ; ijec(m) = 1 

    ! for diagnosing vertically integrated tracer flux
#ifndef STATIC_MEMORY
    allocate (flux_int_z(isd:ied,jsd:jed,4))
#endif
    flux_int_z = 0.0


    ! register/send diagnostics 
    id_slope_x = register_static_field ('ocean_model', 'slope_x', Grd%tracer_axes_flux_x(1:2), &
                 '|d(ht)/dx| on T-cell face', 'm/m', missing_value=-1.e10, range=(/-1.e10,1.e10/))
    if (id_slope_x > 0) used = send_data (id_slope_x, slope_x(isc:iec,jsc:jec), &
                               Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

    id_slope_y = register_static_field ('ocean_model', 'slope_y', Grd%tracer_axes_flux_y(1:2), &
                 '|d(ht)/dy| on T-cell face', 'm/m', missing_value=-1.e10, range=(/-1.e10,1.e10/))
    if (id_slope_y > 0) used = send_data (id_slope_y, slope_y(isc:iec,jsc:jec), &
                               Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

    id_topog_step_1 = register_static_field ('ocean_model', 'topog_step_1', Grd%tracer_axes(1:2), &
                 'topog_step_1', 'dimensionless', missing_value=-1.0, range=(/-1.0,1.0/))
    if (id_topog_step_1 > 0) used = send_data (id_topog_step_1, topog_step(isc:iec,jsc:jec,1), &
                                    Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

    id_topog_step_2 = register_static_field ('ocean_model', 'topog_step_2', Grd%tracer_axes(1:2), &
                 'topog_step_2', 'dimensionless', missing_value=-1.0, range=(/-1.0,1.0/))
    if (id_topog_step_2 > 0) used = send_data (id_topog_step_2, topog_step(isc:iec,jsc:jec,2), &
                                    Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

    id_topog_step_3 = register_static_field ('ocean_model', 'topog_step_3', Grd%tracer_axes(1:2), &
                 'topog_step_3', 'dimensionless', missing_value=-1.0, range=(/-1.0,1.0/))
    if (id_topog_step_3 > 0) used = send_data (id_topog_step_3, topog_step(isc:iec,jsc:jec,3), &
                                    Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

    id_topog_step_4 = register_static_field ('ocean_model', 'topog_step_4', Grd%tracer_axes(1:2), &
                 'topog_step_4', 'dimensionless', missing_value=-1.0, range=(/-1.0,1.0/))
    if (id_topog_step_4 > 0) used = send_data (id_topog_step_4, topog_step(isc:iec,jsc:jec,4), &
                                    Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))


    id_overflow_xflux = register_diag_field ('ocean_model', 'overflow_xflux', Grd%tracer_axes_flux_x(1:2), &
         Time%model_time, 'off-shelf i-transport from overflow', &
        'Sv', missing_value=-1.e10, range=(/-1.e10,1.e10/))

    id_overflow_yflux = register_diag_field ('ocean_model', 'overflow_yflux', Grd%tracer_axes_flux_y(1:2), &
         Time%model_time, 'off-shelf j-transport from overflow', &
        'Sv', missing_value=-1.e10, range=(/-1.e10,1.e10/))


    allocate (id_overflow(num_prog_tracers))
    allocate (id_overflow_xflux_int_z(num_prog_tracers))
    allocate (id_overflow_yflux_int_z(num_prog_tracers))
    id_overflow = -1
    id_overflow_xflux_int_z=-1
    id_overflow_yflux_int_z=-1

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_overflow(n) = register_diag_field ('ocean_model', &
                'overflow_'//trim(T_prog(n)%name), &
                Grd%tracer_axes(1:3), Time%model_time, &
                'rho_cp*overflow*dht*temp', &
                'Watt/m^2', missing_value=-1.e10, range=(/-1.e10,1.e10/))
           id_overflow_xflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_overflow_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time, &
              'vert-integ rho_cp*overflow_xflux*dyt*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
           id_overflow_yflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_overflow_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time, &
              'vert-integ rho_cp*overflow_yflux*dxt*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
       else
           id_overflow(n) = register_diag_field ('ocean_model', 'overflow_'//trim(T_prog(n)%name), &
                Grd%tracer_axes(1:3), Time%model_time, &
                'rho0*overflow*dht*tracer for '//trim(T_prog(n)%name),&
                trim(T_prog(n)%units)//' m*kg/sec', missing_value=-1.e10, range=(/-1.e10,1.e10/))
           id_overflow_xflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_overflow_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time, &
              'vert-integ rho0*overflow_xflux*dyt*dht*tracer for', &
              'kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
           id_overflow_yflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_overflow_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time, &
              'vert-integ rho0*overflow_yflux*dxt*dht*tracer for', &
              'kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/)) 
       endif
    enddo

    if (debug_overflow) then
       write(stdout(),*) ' '
       write(stdout(),*) '==Global sums from ocean_overflow_mod== '
       call write_timestamp(Time%model_time)
       write(stdout(),'(a,es24.17)') 'slope_x        = ',mpp_global_sum(Dom%domain2d,slope_x(:,:), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'slope_y        = ',mpp_global_sum(Dom%domain2d,slope_y(:,:), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_slope(1) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,1), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_slope(2) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,2), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_slope(3) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,3), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_slope(4) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,4), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_step(1)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,1), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_step(2)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,2), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_step(3)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,3), global_sum_flag)
       write(stdout(),'(a,es24.17)') 'topog_step(4)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,4), global_sum_flag)
    endif

  end subroutine ocean_overflow_init
! </SUBROUTINE> NAME="ocean_overflow_init"


!#######################################################################
! <SUBROUTINE NAME="overflow">
!
! <DESCRIPTION>
! Compute thickness weighted tracer source [tracer*m/s]
! due to upstream tracer advection in regions where 
! density-driven overflows are favorable. 
!
! The MOM4 implementation of the Campin and Goosse (1999)
! algorithm is detailed in the MOM4 Technical Guide.
!
! </DESCRIPTION>
!
subroutine overflow (Time, Thickness, T_prog, Dens, index_temp, index_salt, dtime)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(inout)  :: T_prog(:)
  type(ocean_density_type), intent(in)         :: Dens
  integer, intent(in)                          :: index_temp
  integer, intent(in)                          :: index_salt
  real,    intent(in)                          :: dtime

  integer :: tau, taum1
  integer :: i, j, k, m, n
  integer :: kup(isd:ied,jsd:jed,4),kdw(isd:ied,jsd:jed,4) 

  real    :: temp_so, salt_so, press, density_check
  real    :: overflow_thickness, overflow_speed
  real    :: source_check(0:nk)
  real    :: source_column(nk)
  real    :: flux(0:nk)
  real    :: overflow_xflux(isd:ied,jsd:jed)
  real    :: overflow_yflux(isd:ied,jsd:jed)
  real    :: overflow_flux(isd:ied,jsd:jed,4) 
  real    :: tmp(isd:ied,jsd:jed) 

  if(.not. use_overflow) return 

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_overflow_mod (overflow): module must be initialized')
  endif 

  tau     = Time%tau
  taum1   = Time%taum1

  ! compute grid details for cells participating in downslope flow 
  ! note:  "so" = "shallow ocean" cell 
  ! note:  "do" = cells within the "deep ocean" column 

  kup(:,:,:) = 0
  kdw(:,:,:) = 0
  overflow_flux(:,:,:) = 0.0
  overflow_xflux(:,:)  = 0.0
  overflow_yflux(:,:)  = 0.0

  do m=1,4
    do j=jsc,jec
       do i=isc,iec

           ! check to see if downslope flow is topographically possible
           if(topog_step(i,j,m) == 1.0) then  

               ! check to see if density of so-cell > density of adjacent do-cell  
               k = Grd%kmt(i,j) 
               if(Dens%rho(i,j,k) > Dens%rho(i+ip(m),j+jq(m),k)) then   

                   ! k-level of shallow-ocn (so) bottom cell
                   kup(i,j,m) = Grd%kmt(i,j) 

                   ! speed (m/sec) of downslope flow 
                   k=kup(i,j,m) 
                   temp_so     = T_prog(index_temp)%field(i,j,k,tau)         
                   salt_so     = T_prog(index_salt)%field(i,j,k,tau)           
                   overflow_speed = overflow_factor*topog_slope(i,j,m)*(Dens%rho(i,j,k)-Dens%rho(i+ip(m),j+jq(m),k))
                   overflow_speed = min(overflow_speed,overflow_umax)

                   ! convert overflow_speed to volume flux (m^3/sec) 
                   k=kup(i,j,m)
                   overflow_thickness = min(Thickness%dht(i,j,k,tau),Thickness%dht(i+ip(m),j+jq(m),k,tau)) 
                   if(m==1 .or. m==3) overflow_flux(i,j,m) = overflow_speed*Grd%dyt(i,j)*overflow_thickness
                   if(m==2 .or. m==4) overflow_flux(i,j,m) = overflow_speed*Grd%dxt(i,j)*overflow_thickness

                   ! kdw = the k-level within the do-column where so-cell is neutrally buoyant 
                   kdw(i,j,m) = kup(i,j,m)
                   do k=kup(i,j,m)+1,Grd%kmt(i+ip(m),j+jq(m))
                      press         = Dens%pressure_at_depth(i+ip(m),j+jq(m),k)
                      density_check = density(salt_so,temp_so,press)
                      if(density_check > Dens%rho(i+ip(m),j+jq(m),k)) then 
                          kdw(i,j,m)=k
                      else 
                          exit
                      endif
                   enddo

                   ! for diagnostics: volume transport (m^3/sec) off the shelf
                   if(m==1) overflow_xflux(i,j) = overflow_xflux(i,j) + overflow_flux(i,j,m)
                   if(m==2) overflow_yflux(i,j) = overflow_yflux(i,j) + overflow_flux(i,j,m)
                   if(m==3) overflow_xflux(i,j) = overflow_xflux(i,j) - overflow_flux(i,j,m)
                   if(m==4) overflow_yflux(i,j) = overflow_yflux(i,j) - overflow_flux(i,j,m)

               endif    ! if-test for density difference
           endif      ! if-test for topog_step(i,j,m)==1.0
        enddo      ! i-end  
     enddo      ! j-end
  enddo      ! m-end for the four surrounding cells 

  call mpp_update_domains(overflow_flux(:,:,:),Dom%domain2d)
  call mpp_update_domains(kup(:,:,:),Dom%domain2d)
  call mpp_update_domains(kdw(:,:,:),Dom%domain2d)

  if (debug_overflow) then
     write(stdout(),*) ' '
     write(stdout(),*) '==Global sums from ocean_overflow_mod== '
     call write_timestamp(Time%model_time)
     write(stdout(),'(a,es24.17)')  'overflow_flux(1)(m^3/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,1), global_sum_flag)
     write(stdout(),'(a,es24.17)')  'overflow_flux(2)(m^3/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,2), global_sum_flag)
     write(stdout(),'(a,es24.17)')  'overflow_flux(3)(m^3/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,3), global_sum_flag)
     write(stdout(),'(a,es24.17)')  'overflow_flux(4)(m^3/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,4), global_sum_flag)
     write(stdout(),'(a,i22)')      'kup(1)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,1), global_sum_flag)
     write(stdout(),'(a,i22)')      'kup(2)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,2), global_sum_flag)
     write(stdout(),'(a,i22)')      'kup(3)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,3), global_sum_flag)
     write(stdout(),'(a,i22)')      'kup(4)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,4), global_sum_flag)
     write(stdout(),'(a,i22)')      'kdw(1)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,1), global_sum_flag)
     write(stdout(),'(a,i22)')      'kdw(2)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,2), global_sum_flag)
     write(stdout(),'(a,i22)')      'kdw(3)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,3), global_sum_flag)
     write(stdout(),'(a,i22)')      'kdw(4)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,4), global_sum_flag)
  endif

  ! compute fluxes and the source for those 
  ! cells participating in downslope transport 
  ! tracer source = convergence of upstream advective fluxes 

  do n=1,num_prog_tracers 

     source_overflow(:,:,:) = 0.0
     flux_int_z(:,:,:)      = 0.0 
     wrk1=0.0 ; wrk2=0.0 ; wrk3=0.0 ; wrk4=0.0 

     ! extended do-loop limits are to reduce the need to call 
     ! mpp_update_domains by performing extra computations.
     do m=1,4
        do j=jsc-ijsc(m),jec+ijec(m)
           do i=isc-ijsc(m),iec+ijec(m)

              ! see if (i,j) is a so-cell participating in downslope flow in the m-direction
              if(kup(i,j,m) > 0) then 

                  ! initialize upstream tracer flux leaving tracer cells   
                  flux(:) = 0.0

                  ! initialize thickness weighted tendency in column 
                  source_column(:) = 0.0

                  ! flux horizontally leaving so-cell at k=kup and entering do-cell at k=kdw 
                  k=kup(i,j,m) 
                  flux(0) = overflow_flux(i,j,m)*T_prog(n)%field(i,j,k,taum1)

                  ! flux leaving do-cells
                  do k=kup(i,j,m),kdw(i,j,m)
                      flux(k) = overflow_flux(i,j,m)*T_prog(n)%field(i+ip(m),j+jq(m),k,taum1)
                  enddo
 
                  ! for diagnosing vertically integrated horizontal volume flux   
                  flux_int_z(i,j,m) = flux_int_z(i,j,m) + flux(0)-flux(kup(i,j,m))                   

                  ! source for so-cell at k=kup   
                  k = kup(i,j,m)
                  source_overflow(i,j,k) = source_overflow(i,j,k) + Grd%datr(i,j)*(flux(k)-flux(0)) 

                  ! source for do-cell at k=kdw   
                  k = kdw(i,j,m)
                  source_column(k) = Grd%datr(i+ip(m),j+jq(m))*(flux(0)-flux(k))

                  ! source for do-cell with kup <= k <= kdw-1
                  if(kdw(i,j,m) > kup(i,j,m)) then 
                      do k = kup(i,j,m),kdw(i,j,m)-1
                         source_column(k) = Grd%datr(i+ip(m),j+jq(m))*(flux(k+1)-flux(k)) 
                      enddo
                  endif

                  if(m==1) wrk1(i,j,:) = source_column(:) 
                  if(m==2) wrk2(i,j,:) = source_column(:) 
                  if(m==3) wrk3(i,j,:) = source_column(:) 
                  if(m==4) wrk4(i,j,:) = source_column(:)


              endif  ! kup(i,j,m) > 0
           enddo ! i-end 
        enddo  ! j-end
     enddo  ! m-end for the four directions 

     ! each cell can have a source five different ways:
     ! 1:   when cell is the so-cell
     ! 2-5: when cell is one of the do-cells for adjacent so-cells.
     do k=1,nk
        do j=jsc,jec   
           do i=isc,iec
              source_overflow(i,j,k)  = source_overflow(i,j,k)   &
                                        + wrk1(i-1,j,k) & 
                                        + wrk2(i,j-1,k) &
                                        + wrk3(i+1,j,k) &
                                        + wrk4(i,j+1,k)
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + source_overflow(i,j,k)
           enddo
        enddo
     enddo

     if (debug_overflow) then
        source_check(:) = 0.0
        do k=1,nk
           do j=jsc,jec  
             do i=isc,iec
               tmp(i,j) =  T_prog(n)%conversion*dtime*Grd%dat(i,j)*source_overflow(i,j,k)
             enddo
           enddo    
           source_check(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
        enddo
        do k=1,nk
           source_check(0) = source_check(0) + source_check(k)
        enddo
        if(T_prog(n)%name == 'temp' ) then 
           write(stdout(),'(a,i3,a,es24.17)') 'For tracer(',n,'), (temp   ) total overflow source is (J) ', source_check(0)
        elseif(T_prog(n)%name == 'salt' ) then 
           write(stdout(),'(a,i3,a,es24.17)') 'For tracer(',n,'), (salt   ) total overflow source is (kg)', source_check(0)
        else 
           write(stdout(),'(a,i3,a,es24.17)') 'For tracer(',n,'), (passive) total overflow source is (kg)', source_check(0)
        endif 
     endif

     if(id_overflow(n) > 0) then 
         used = send_data (id_overflow(n), T_prog(n)%conversion*source_overflow(:,:,:), &
                Time%model_time, rmask=Grd%tmask(:,:,:), &
                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_overflow_xflux_int_z(n) > 0) then 
         used = send_data (id_overflow_xflux_int_z(n), &
                T_prog(n)%conversion*(flux_int_z(:,:,1)-flux_int_z(:,:,3)), &
                Time%model_time, rmask=Grd%tmask(:,:,1), &
                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 
     if(id_overflow_yflux_int_z(n) > 0) then 
         used = send_data (id_overflow_yflux_int_z(n), &
                T_prog(n)%conversion*(flux_int_z(:,:,2)-flux_int_z(:,:,4)), &
                Time%model_time, rmask=Grd%tmask(:,:,1), &
                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif


  enddo  ! n-end for num_prog_tracers 


  ! 1e-6 converts (m^3/sec) to Sv
  if(id_overflow_xflux > 0) then 
    used = send_data (id_overflow_xflux, 1e-6*overflow_xflux(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1), &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if(id_overflow_yflux > 0) then 
    used = send_data (id_overflow_yflux, 1e-6*overflow_yflux(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1), &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


end subroutine overflow
! </SUBROUTINE> NAME="overflow"


end module ocean_overflow_mod
      
      




