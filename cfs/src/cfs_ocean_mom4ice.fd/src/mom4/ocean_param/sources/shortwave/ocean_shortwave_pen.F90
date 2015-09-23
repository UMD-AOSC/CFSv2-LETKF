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
module ocean_shortwave_pen_mod
!
!<CONTACT EMAIL="Tony.Rosati@noaa.gov"> A. Rosati
!</CONTACT>
!
!<CONTACT EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Colm.Sweeney@noaa.gov"> Colm Sweeney
!</CONTACT>
!
!<REVIEWER EMAIL="russell.fiedler@csiro.au"> Russell Fiedler 
!</REVIEWER>
!
!<OVERVIEW>
! This module returns thickness weighted temperature 
! tendency [deg C *m/sec] from penetrative shortwave heating.
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted tendency [deg C *m/sec]
! of tracer associated with penetrative shortwave heating in the upper
! ocean. Generally penetration is taken as a function of monthly optical 
! properties of the upper ocean, where optical properties are read 
! in from a file of climatological data. Presently there is account taken 
! only of chlorophyll-a on the optical properties of ocean water.  Other 
! particulates can be added so to have a more complete picture of the ocean 
! optical properties.  Also, this module provide a framework for 
! incorporating the effects from a prognostic biology model on ocean optics.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Jerlov (1968)
! Optical Oceanography
! Elsevier Press
! </REFERENCE>
!
! <REFERENCE>
! Morel and Antoine (1994)
! Heating rate in the upper ocean in relation to its bio-optical state 
! Journal of Physical Oceanography vol 24 pages 1652-1664
! </REFERENCE>
!
! <REFERENCE>
! Paulson and Simpson (1977)
! Irradiance measurements in the upper ocean
! Journal of Physical Oceanography vol 7 pages 952-956
! </REFERENCE>
!
! <REFERENCE>
! Rosati and Miyakoda (1988)
! A General Circulation Model for Upper Ocean Simulation
! Journal of Physical Oceanography vol 18 pages 1601-1626.
! </REFERENCE>
!
! <NOTE>
! Optimized for vector peformance by R. Fiedler (russell.fiedler@csiro.au)
! June 2003 on the Australian NEC computer. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_shortwave_pen_nml">
!  <DATA NAME="shortwave_pen_on=" TYPE="logical">
!  Must be .true. to run with module. 
!  </DATA> 
!  <DATA NAME="read_chl" TYPE="logical">
!  If .true. then read in climatological data of chlorophyll-a.
!  </DATA> 
!  <DATA NAME="sw_frac_top" TYPE="real">
!  The fraction of shortwave radiation that should be incorporated into 
!  the sw_source array at k=1.  The generic treatment in mom4 is to assume
!  that shortwave radiation is already contained inside the 
!  T_prog(index_temp)%stf field. Hence, to avoid   
!  double counting, sw_frac(k=0)=sw_frac_top should=0.0.
!  If one removes shortwave from stf, then set sw_frac_top=1.0.
!  </DATA> 
!  <DATA NAME="zmax_pen" UNITS="meter" TYPE="real">
!   Maximum depth of penetration of shortwave radiation. 
!   Below this depth, shortwave penetration is exponentially 
!   small and so is ignored.
!  </DATA>
!  <DATA NAME="chl_default" UNITS="mg/m^3" TYPE="real">
!   Default concentration chl_default=0.08 roughly yields Jerlov Type 1A optics.
!  </DATA>
!  <DATA NAME="enforce_sw_frac" TYPE="logical">
!  To ensure the shortwave fraction is monotonically decreasing with depth. 
!  </DATA> 
!  <DATA NAME="debug_sw_pen" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>

use axis_utils_mod,           only: frac_index
use constants_mod,            only: rho_cp
use diag_manager_mod,         only: register_diag_field, send_data
use field_manager_mod,        only: fm_get_index
use fms_mod,                  only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,                  only: stdout, stdlog, FATAL, NOTE
use mpp_mod,                  only: mpp_error, mpp_max, mpp_min
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_ROUTINE
use time_interp_external_mod, only: time_interp_external, init_external_field

use ocean_domains_mod,        only: get_local_indices
use ocean_tpm_util_mod,       only: T_diag
use ocean_types_mod,          only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,          only: ocean_prog_tracer_type

implicit none

private

! for diagnostics 
integer :: id_sw_frac=-1
integer :: id_sat_chl=-1
integer :: id_sw_heat=-1
logical :: used

! clock ids
integer :: id_sw_pen

integer :: index_chl
integer :: index_irr
integer :: sbc_chl
logical :: verbose_flag=.false.

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
  real, dimension(isd:ied,jsd:jed)  :: sw_fk_zt   ! sw (radiation) fractional decay on t grid
  real, dimension(isd:ied,jsd:jed)  :: sw_fk_zw   ! sw (radiation) fractional decay on w grid
  real, dimension(isd:ied,jsd:jed)  :: sat_chl    ! chlorophyll concentration (mg/m^3) (a proxy for color)   
#else  
  real, allocatable, dimension(:,:) :: sw_fk_zt   ! sw (radiation) fractional decay on t grid
  real, allocatable, dimension(:,:) :: sw_fk_zw   ! sw (radiation) fractional decay on w grid
  real, allocatable, dimension(:,:) :: sat_chl    ! chlorophyll concentration (mg/m^3) (a proxy for color)   
#endif
  
! defining coeficients for optical model (these coeficients represent a non uniform distribution
! of chlorophyll-a through the water column) taken from Morel and Antoine (1994):
! These may be more appropriate when using an interactive ecosystem model that predicts
! three-dimensional chl-a values.  
!
real, dimension(6), parameter::V1_coef=(/0.321,   0.008,   0.132,   0.038, -0.017,  -0.007/)
real, dimension(6), parameter::V2_coef=(/0.679,  -0.008,  -0.132,  -0.038,  0.017,   0.007/)
real, dimension(6), parameter::Z1_coef=(/1.540,  -0.197,   0.166,  -0.252, -0.055,   0.042/)
real, dimension(6), parameter::Z2_coef=(/7.925,  -6.644,   3.662,  -1.815, -0.218,   0.502/)
!
! defining coeficients for optical model (these coeficients represent a uniform distribution
! of chlorophyll-a through the water column) taken from Morel and Antoine (1994):
! These may be more appropriate when using chl-a from surface-viewing satellite data. 
!
!!$real, dimension(6), parameter::V1_coef=(/0.353,  -0.047,  0.083,  0.047, -0.011, -0.009/)
!!$real, dimension(6), parameter::V2_coef=(/0.647,   0.047, -0.083, -0.047,  0.011,  0.009/)
!!$real, dimension(6), parameter::Z1_coef=(/1.662,  -0.605,  0.128, -0.033, -0.051, -0.004/)
!!$real, dimension(6), parameter::Z2_coef=(/8.541,  -8.924,  4.020, -0.077, -0.536,  0.055/)


type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$Id$'
character (len=128) :: tagname = '$Name$'

public  ocean_shortwave_pen_init
private sw_pen
public  sw_source


logical :: shortwave_pen_on       = .false.
logical :: read_chl               = .false.
logical :: module_is_initialized  = .FALSE.
logical :: debug_sw_pen           = .false. 
logical :: enforce_sw_frac        = .true. 

real :: chl_default = 0.08  ! (mg/m^3) default concentration 0.08 roughly yields Jerlov Type 1A optics
real :: zmax_pen    = 120.0 ! maximum depth (m) of solar penetration. 
                            ! below, penetration is exponentially small and so is ignored
real :: sw_frac_top = 0.0   ! set to 1.0 if do not have shortwave radiation inside of T_prog(index_temp)%stf.

namelist /ocean_shortwave_nml/ shortwave_pen_on, read_chl, chl_default, &
                               zmax_pen, sw_frac_top, debug_sw_pen,     &
                               enforce_sw_frac

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_shortwave_pen_init">
!
! <DESCRIPTION>
! Initialization for the shorwave module
! </DESCRIPTION>
  subroutine ocean_shortwave_pen_init(Grid, Domain, Time)

    type(ocean_grid_type), intent(in), target    :: Grid
    type(ocean_domain_type), intent(in), target  :: Domain
    type(ocean_time_type), intent(in)            :: Time

    integer :: unit, io_status, ierr, i, j
#ifdef STATIC_MEMORY    
    real, dimension(isd:ied,jsd:jed)            :: chl_data
#else
    real, allocatable, dimension(:,:)           :: chl_data
#endif

    if ( module_is_initialized ) return
    
    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    unit = open_namelist_file()
    read(unit, ocean_shortwave_nml,iostat=io_status)
    write (stdout(),'(/)')
    write(stdout(),ocean_shortwave_nml)    
    write(stdlog(),ocean_shortwave_nml)
    ierr = check_nml_error(io_status, 'ocean_shortwave_nml')
    call close_file(unit)

    Dom => Domain
    Grd => Grid

#ifndef STATIC_MEMORY    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif

    if(shortwave_pen_on) then 
      call mpp_error(NOTE, '==>Note: USING shortwave_pen_mod.')
    else
      call mpp_error(NOTE, '==>Note: NOT using shortwave_pen_mod.')
      return
    endif

#ifndef STATIC_MEMORY   
    allocate( chl_data(isd:ied,jsd:jed))
    allocate( sat_chl(isd:ied,jsd:jed))
    allocate( sw_fk_zt(isd:ied,jsd:jed))
    allocate( sw_fk_zw(isd:ied,jsd:jed))
#endif

    ! set clock ids     
    id_sw_pen = mpp_clock_id('(Ocean shortwave penetrate) ' ,grain=CLOCK_ROUTINE)

    sat_chl(:,:)   = chl_default*Grd%tmask(:,:,1)
    sw_fk_zt          = 0.0
    sw_fk_zw          = 0.0

    index_chl = fm_get_index('/ocean_mod/diag_tracers/chl')
    index_irr = fm_get_index('/ocean_mod/diag_tracers/irr')

    if(read_chl) then 

       call mpp_error(NOTE, '==>Note: Reading in chlorophyll-a from data file for shortwave penetration.')

       ! get the unit number "sbc_chl" for reading chl data 
       sbc_chl = init_external_field('INPUT/chl','chl',domain=Domain%domain2d)
       if (sbc_chl == -1) call mpp_error(FATAL,&
            '==>Error in ocean_shortwave_pen_mod: failure to find sbc_chl data file')

       ! Need to update chl in case of restart
       chl_data = 0.0
       call time_interp_external(sbc_chl, Time%model_time, chl_data, verbose=debug_sw_pen)
       do j=jsc,jec
          do i=isc,iec
             sat_chl(i,j)    = Grd%tmask(i,j,1)*chl_data(i,j)
          enddo
       enddo

       id_sat_chl = register_diag_field ('ocean_model', 'sat_chl', &
            Grid%tracer_axes(1:2), Time%model_time, 'Chlorophyll', &
            'mg/m^3',missing_value=-10.0, range=(/-10.0,10.0/))

    elseif (index_chl .gt. 0) then
       call mpp_error(NOTE, '==>Note: Using chl from biology module for shortwave penetration.')
    else
       call mpp_error(NOTE, '==>Note: Setting chl=chl_default in shortwave_pen_init.')
    endif

    if(sw_frac_top==0.0) then 
       write(stdout(),*)'=>Note: computing effects of solar shortwave penetration.  Model assumes stf has sw-radiation'
       write(stdout(),*)'  field included.  Hence, solar shortwave penetration effects placed in sw_source will '
       write(stdout(),*)'  subtract out the effects of shortwave at k=1 to avoid double-counting.'
    elseif(sw_frac_top==1.0) then 
       write(stdout(),*)'=>Note: computing effects of solar shortwave penetration. We assume stf does not have sw-radiation'
       write(stdout(),*)' field included.  Hence, solar shortwave penetration effects are placed completely in sw_source.'
       write(stdout(),*)' This is not the standard approach used in mom4.'
    elseif(sw_frac_top/=1.0 .and. sw_frac_top/=0.0) then 
       write(stdout(),*)'=>Note: computing effects of solar shortwave penetration. We assume a portion of the sw-effects are'
       write(stdout(),*)'  included in stf and a portion in sw_source.  Are you sure you wish to partition shortwave'
       write(stdout(),*)'  in such a manner?'
    endif

    if(enforce_sw_frac) then  
       write(stdout(),*)'=>Note: enforce_sw_frac=.true. so will enforce monotonic decrease of sw_frac with depth.'
    else 
       write(stdout(),*)'=>Note: enforce_sw_frac=.false. so can have non-monotonic sw_frac with certain penetration profiles.'
    endif

    id_sw_frac = register_diag_field ('ocean_model', 'sw_frac', &
         Grid%tracer_axes(1:3), Time%model_time, 'sw_frac from shortwave pen ', &
        ' frac', missing_value=-1e10, range=(/-1.e10,1.e10/))

    id_sw_heat = register_diag_field ('ocean_model', 'sw_heat', &
         Grid%tracer_axes(1:3), Time%model_time, 'thickness wghtd shortwave heating', &
         'W/m^2', missing_value=-1e10, range=(/-1.e10,1.e10/))

#ifndef STATIC_MEMORY
    deallocate(chl_data)
#endif


end subroutine ocean_shortwave_pen_init
! </SUBROUTINE> NAME="ocean_shortwave_pen_init"



!#######################################################################
! <SUBROUTINE NAME="sw_source">
!
! <DESCRIPTION>
! Incorporate short wave penetration via the "source" 
! term. note that the divergence of shortwave for the first
! level "div_sw(0)" is compensating for the effect of having
! the shortwave component already included in the total
! surface tracer flux "stf(i,j,temp)"
!
! output:
!
! tracer_source = thickness weighted source from  penetrative short wave heating
!
! </DESCRIPTION>
subroutine sw_source (Time, swflx, Temp, sw_frac_zt)

  type(ocean_time_type), intent(in)                    :: Time
  real, intent(in), dimension(isd:ied,jsd:jed)         :: swflx
  type(ocean_prog_tracer_type), intent(inout)          :: Temp
  real, intent(inout), dimension(isd:ied,jsd:jed,1:nk) :: sw_frac_zt
  real, dimension(isd:ied,jsd:jed,0:nk)                :: sw_frac_zw
  real, dimension(isd:ied,jsd:jed)                     :: zt_sw
  real, dimension(isd:ied,jsd:jed)                     :: zw_sw
  real, dimension(isd:ied,jsd:jed)                     :: chl_data
  real    :: div_sw 
  integer :: i, j, k, irr_index

  if(.not. shortwave_pen_on) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error in ocean_shortwave_pen_mod (sw_source): module must be initialized ')
  endif 

  if (Temp%name /= 'temp') then 
    call mpp_error(FATAL,'==>Error in ocean_shortwave_pen_mod (sw_source): invalid tracer for sw_source')
  endif 

  if (read_chl) then
    chl_data=0.0
    call time_interp_external(sbc_chl, Time%model_time, chl_data, verbose=debug_sw_pen)
    do j=jsc,jec
      do i=isc,iec
        sat_chl(i,j) = Grd%tmask(i,j,1)*chl_data(i,j)
      enddo
    enddo
    if (id_sat_chl > 0) used = send_data (id_sat_chl, sat_chl(:,:), &
                               Time%model_time,rmask=Grd%tmask(:,:,1), &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! zero out the wrk1 array used to diagnose heating from shortwave 
  Temp%wrk1(:,:,:) = 0

  ! only compute 3-D sw_fract for ocean regions shallower than zmax_pen 
  zw_sw=0.0
  zt_sw=0.0
  sw_frac_zw(:,:,0) = sw_frac_top
  do k=1,nk-1
     if(Grd%zw(k) <= zmax_pen) then 
        zw_sw(isc:iec,jsc:jec) = Grd%zw(k)
        call sw_pen(zw_sw, sw_fk_zw)
        zt_sw(isc:iec,jsc:jec) = Grd%zt(k)
        call sw_pen(zt_sw, sw_fk_zt)
     else
       sw_fk_zt(:,:) = 0.0
       sw_fk_zw(:,:) = 0.0
     endif 
     sw_frac_zt(:,:,k) = sw_fk_zt(:,:)
     sw_frac_zw(:,:,k) = sw_fk_zw(:,:)
  enddo

  if(enforce_sw_frac) then   
      do k=2,nk-1
         do j=jsc,jec
            do i=isc,iec
               sw_frac_zt(i,j,k) = min(sw_frac_zt(i,j,k),sw_frac_zt(i,j,k-1))
            enddo
         enddo
      enddo
  endif
  
  ! when chlorophyll is being read through T_diag (not from climatology)
  if (index_chl .gt. 0) then  !{!
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
          T_diag(index_irr)%field(i,j,k) = swflx(i,j) * sw_frac_zt(i,j,k) * rho_cp
        enddo
       enddo
     enddo
  endif

  ! divergence of solar shortwave fraction 
  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        div_sw                  = (sw_frac_zw(i,j,k-1) - sw_frac_zw(i,j,k)*Grd%tmask(i,j,k+1))/Grd%dzt(k)
        Temp%wrk1(i,j,k)        = Grd%tmask(i,j,k)*swflx(i,j)*div_sw*Grd%dzt(k)
        Temp%th_tendency(i,j,k) = Temp%th_tendency(i,j,k) + Temp%wrk1(i,j,k)
      enddo
    enddo
  enddo

  if (id_sw_frac > 0) used = send_data (id_sw_frac, sw_frac_zt(:,:,1:nk), &
                             Time%model_time, rmask=Grd%tmask(:,:,:), &
                             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_sw_heat > 0) used = send_data (id_sw_heat, Temp%wrk1(:,:,:)*Temp%conversion, &
                             Time%model_time, rmask=Grd%tmask(:,:,:), &
                             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

end subroutine sw_source
! </SUBROUTINE> NAME="sw_source"


!#######################################################################
! <SUBROUTINE NAME="sw_pen">
!
! <DESCRIPTION>
!  Solar shortwave energy penetrates below the ocean surface and is aborbed 
!  by water and organic matter (both particulate and dissolved). This
!  routine estimates fraction of shortwave penetration using chlorophyll a.
!  Absorbtion of shortwave radiation in the water assumes energy partitions
!  between three exponentials:
!
!  The first exponential is for wavelength > 0.75 um (microns) and assumes a
!  single attenuation of 0.267 m if the "zenith_angle" is 0.  Presently the 
!  code assumes a zero zenith angle, but this could be modified easily. 
!
!  The second and third exponentials represent a parameterization of the
!  attenuation coeficient for light between 300 um and 750 um in the following
!  form:
!
!	E(z) = E(0) * [V1 *  exp(z/efold1) + V2 * exp(z/efold2)]
!       with z < 0 the ocean depth 
!
!  Here, V1+V2=1 represent the partitioning between long (V1) and short (V2)
!  wavelengths between 300 um and 750 um. Thoughout most of the ocean V1<0.5
!  and V2>0.5. The "efold1" and "efold2" are the efolding depth of the long and short
!  visable and ultra violet light. Throughout most of the ocean efold1 should not exceed 3 m
!  while the efold2 will vary between 30 m in oligotrophic waters and 4 m in coastal
!  regions. All of these constants are based on satellite estimates of chlorophyll a and
!  taken from Morel and Antoine (JPO 1994, (24) 1652-1664).
!
!  If the thickness of the first ocean level "dzt(1)" is 50 meters,
!  then shortwave penetration does not do much. However, for higher
!  vertical resolution, such as dzt(1) = 10 meters commonly used
!  in ocean climate models, the effect of shortwave heating can
!  be significant. This can be particularly noticable in the summer
!  hemisphere.
!
! </DESCRIPTION>
!
! <INFO>
!
! <NOTE>
!  z_sw is a 2d array since KPP calls sw_pen and passes hbl(i,j) into z_sw. 
! </NOTE> 
!
! <NOTE>
!  The terms contributing to sw_fk(i,j) are depth independent
!  when chl is depth independent.  However, we anticipate implementing 
!  a biological model, whereby chl will be depth dependent.  
! </NOTE> 
!
! <NOTE>
!  Simpson and Dickey (1981) and others have argued between one and 
!  two exponentials for light between 300 um and 750 um.  
!  With vertical grid resolution of 5 meters or finer
!  for the upper 20 meters, the second exponential will make a difference.
!  We anticipate using such resolutions, and so have implemented both 
!  exponentials. 
! </NOTE> 
!
! </INFO>
!
subroutine sw_pen (z_sw, sw_fk)

  real, intent(in),    dimension(isd:ied,jsd:jed) :: z_sw     ! vertical depth
  real, intent(inout), dimension(isd:ied,jsd:jed) :: sw_fk    ! sw fractional decay

! flag for setting k of bottom of boundary layer           
  logical :: keep_going(isd:iec)                            
! introducing kb as 1d vector improves vector performance  
  integer, dimension(isc:iec) :: kb   
  integer, dimension(isc:iec) :: kb_old   

! parameters for optical model using C=log10(chl)
  real    :: V1, V2, efold1, efold2 
  real    :: C, C2, C3, C4, C5
  real    :: F_vis, zenith_angle
  real    :: swmax, swmin, chmax, chmin
  integer :: i, j, k

  call mpp_clock_begin(id_sw_pen)  

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_shortwave_pen_mod (sw_pen): module must be initialized')
  endif 

! F_vis is the amount of light in the shortwave verses the long wave. F_vis=0.54  
! on sunny days and F_vis=0.60 on cloudy days. 
  F_vis        = 0.57     
  zenith_angle = 0.0

  ! split the i,j loops to allow vectorization in j-loop 
  do j=jsc,jec

     ! compute kb
     k = 0
     kb(:)=nk                                   
     keep_going(:) = .true.
     do while (k < nk .and. any(keep_going(:))) 
        k = k+1
        do  i=isc,iec
           if (z_sw(i,j) <= Grd%zw(k) .and.  keep_going(i)) then   
               kb(i) = k
               keep_going(i) = .false.
           endif
        enddo
     enddo
     do i=isc,iec
       kb(i)=min(kb(i)+1,nk)
     enddo

     ! check with older (non-vectorized) method for kb calculation
     if(debug_sw_pen) then 
         do i=isc,iec 
            kb_old(i) = ceiling(frac_index(z_sw(i,j), (/0.,Grd%zw(1:nk)/)))
            kb_old(i) = min(kb_old(i),nk)
         enddo
         do i=isc,iec 
            if(kb(i) - kb_old(i) /= 0) then 
                write(stdout(),*)'In sw_pen, kb computed two ways: kbnew(',i+Dom%ioff,',',j+Dom%joff,')= ', kb(i), &
                                  ' kbold(',i+Dom%ioff,',',j+Dom%joff,')= ',kb_old(i)
            endif
         enddo
     endif

     ! compute shortwave fraction based on triple exponential 
     ! note that 0.02 and 60.0 provide a floor and ceiling which 
     ! keep sw_fk between 0.0 and 1.0.  These values are dependent
     ! on details of the coefficients used in the exponential. 
     ! If the coefficients change, then the floor/ceiling needs 
     ! to be reevaluated. 
     do i=isc,iec

        if ( z_sw(i,j) > zmax_pen .or. Grd%tmask(i,j,kb(i)) == 0) then
            sw_fk(i,j) = 0.0
        else

            if (.not. read_chl .and. index_chl .gt. 0) then      ! chlorophyll-a from biology.
              C = log10(min(max(T_diag(index_chl)%field(i,j,kb(i)),0.02),60.0))
            else                                                 ! chlorophyll-a read from climatology
              C = log10(min(max(sat_chl(i,j),0.02),60.0))
            endif
            C2 = C  * C
            C3 = C2 * C
            C4 = C3 * C
            C5 = C4 * C
            V1 = V1_coef(1) + V1_coef(2) * C + V1_coef(3) * C2  &
                 + V1_coef(4) * C3 + V1_coef(5) * C4 + V1_coef(6) * C5
            V2 = V2_coef(1) + V2_coef(2) * C + V2_coef(3) * C2  &
                 + V2_coef(4) * C3 + V2_coef(5) * C4 + V2_coef(6) * C5
            efold1 = Z1_coef(1) + Z1_coef(2) * C + Z1_coef(3) * C**2  &
                 + Z1_coef(4) * C3 + Z1_coef(5) * C4 + Z1_coef(6) * C5
            efold2 = Z2_coef(1) + Z2_coef(2) * C + Z2_coef(3) * C**2  &
                 + Z2_coef(4) * C3 + Z2_coef(5) * C4 + Z2_coef(6) * C5

            sw_fk(i,j) = (1-F_vis) * exp( -z_sw(i,j)/(0.267 * cos(zenith_angle)) )   &
                 + F_vis  * ( V1 * exp( -z_sw(i,j)/efold1 )                 &
                 + V2 * exp( -z_sw(i,j)/efold2 ) )
        endif

     enddo  ! i-loop finish 

  enddo  ! j-loop finish

  if(debug_sw_pen) then 
      swmax=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_max(swmax);write(stdout(),*)'In ocean_shortwave_pen (sw_pen): max sw_fk=',swmax
      swmin=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_min(swmin);write(stdout(),*)'In ocean_shortwave_pen (sw_pen): min sw_fk=',swmin
      if (index_chl .le. 0) then
        chmax=maxval(sat_chl(isc:iec,jsc:jec))
        call mpp_max(chmax);write(stdout(),*)'In ocean_shortwave_pen (sw_pen): max chl=',chmax
        chmin=maxval(sat_chl(isc:iec,jsc:jec)) 
        call mpp_min(chmin);write(stdout(),*)'In ocean_shortwave_pen (sw_pen): min chl=',chmin
      endif
  endif
  call mpp_clock_end(id_sw_pen)

end subroutine sw_pen
! </SUBROUTINE> NAME="sw_pen"


end module ocean_shortwave_pen_mod
