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
module ocean_lap_friction_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! This module computes the thickness weighted acceleration for 
! horizontal velocity arising from horizontal Laplacian friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for 
! horizontal velocity arising from horizontal Laplacian friction. 
! The viscosity used to determine the strength of the tendency 
! can be a general function of space yet it is constant in time.  
! A namelist option exists that determines this viscosity 
! as a local function of the grid spacing. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
!  R. J. Murray and C. J. C. Reason,
!  A curvilinear version of the Bryan-Cox ocean model
!  Journal Computational Physics (2002), vol 171, pages 1--46
! </REFERENCE>
!
! <REFERENCE>
!  R.C. Pacanowski and A. Gnanadesikan,
!  Transient response in a z-level ocean model that resolves topography 
!  with partial-cells
!  Monthly Weather Review (1998), vol 126, pages 3248-3270
! </REFERENCE>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! <NOTE>
! This scheme has been found to be faster than the Smagorinsky viscosity 
! scheme.  However, the algorithm here is less robust since it 
! contains null modes in the terms associated with the sphericity 
! of the earth.  Hence, there may be flow configurations that are 
! not dissipated.
! </NOTE>
!
! <NOTE>
! The model can generally run with both Laplacian and biharmonic friction
! enabled at the same time.  Such has been found useful for some eddying 
! ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lap_friction_const_nml">
!  <DATA NAME="lap_friction_on" TYPE="logical">
!  Must be true to use this module. 
!  </DATA> 
!  <DATA NAME="alap" UNITS="m^2/sec" TYPE="real">
!  This is the value for the space-time constant Laplacian viscosity. 
!  </DATA> 
!  <DATA NAME="velocity_mix_micom" TYPE="logical">
!  If .true., then the viscosity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use mpp_mod,             only: mpp_sum, mpp_pe, mpp_max, mpp_error

use ocean_types_mod,     only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_velocity_type
use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: FAX, FAY, BAX, BAY, FMX, FMY
use ocean_operators_mod, only: BDX_EU, BDY_NU, FDX_PU, FDY_PU, FDX_NT, FDY_ET


implicit none

public ocean_lap_friction_init
public horz_lap_friction
public horz_lap_viscosity_check
public horz_lap_reynolds_check

private

! clock ids
integer :: id_lap_friction

! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_visc=-1
integer :: id_lap_fric_u=-1
integer :: id_lap_fric_v=-1
logical :: used

real :: alap                  = 0e5     ! constant horz viscosity for momentum (m^2/sec)
real :: vel_micom             = 0.0     ! background scaling velocity for micom viscosity (m/sec) 
real :: alap_taper_start_lat  = 60.     ! for tapering viscosity at high latitudes
real :: alap_taper_end_lat    = 90.     ! for tapering viscosity at high latitudes
real :: alap_taper_min_fac    = 30.     ! for tapering viscosity at high latitudes
logical :: alap_taper_hi_lats = .false. ! for tapering viscosity at high latitudes
logical :: velocity_mix_micom =.false.  ! if true, viscosity made a function of the grid spacing

real, dimension(:,:,:), allocatable :: visc_ceu ! viscosity on east face of U cells (m^2/s)
real, dimension(:,:,:), allocatable :: visc_cnu ! viscosity on north face of U cells (m^2/s)
real, dimension(:,:,:), allocatable :: visc_cu  ! variable viscosity averaged to U cell (m^2/s)

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: /const/ocean_lap_friction.f90 ($Id$)'
character (len=128) :: tagname = &
     '$Name$'

logical :: lap_friction_on       = .false.
logical :: module_is_initialized = .FALSE.

namelist /ocean_lap_friction_const_nml/ lap_friction_on, alap, velocity_mix_micom, vel_micom, &
                                        alap_taper_hi_lats, alap_taper_start_lat, alap_taper_end_lat, &
                                        alap_taper_min_fac

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_lap_friction_init">
!
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
  integer :: i, j, k, num

  integer, save :: unit=6     

  real    :: cbeta, alap_munk, dau_max
  real    :: alap_taper_min            
  real, dimension(Domain%isd:Domain%ied,Domain%jsd:Domain%jed) :: alap_const

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (ocean_lap_friction_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_lap_friction_const_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_lap_friction_const_nml)  
  write (stdlog(),ocean_lap_friction_const_nml)
  ierr = check_nml_error(io_status,'ocean_lap_friction_const_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(lap_friction_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_lap_friction_mod from lap/const.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_lap_friction_mod.')
    return
  endif 

  if(alap > 0 .or. (velocity_mix_micom .and. vel_micom > 0)) then 
    write (stdout(),'(/,a)') '==>USING constant horizontal Laplacian friction'
  else 
    write (stdout(),'(/,a)') '==>NOT using constant horizontal Laplacian friction'
  endif 

  Grd => Grid
  Dom => Domain
  dtime = d_time

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_lap_friction_mod: using forward time step of (secs)', dtime 

  alap_const = alap
  if (velocity_mix_micom) then
      alap_const(:,:) = vel_micom*(2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/&
           (epsln + Grd%dxu(:,:) + Grd%dyu(:,:))
  else
      ! use alternate velocity scaling based on cell area
      dau_max = 0.0
      do j=jsc,jec
         do i=isc,iec
            if (Grd%dau(i,j) > dau_max) dau_max = Grd%dau(i,j)
         enddo
      enddo
      call mpp_max(dau_max)
      alap_const(:,:) = alap*(Grd%dau(:,:)/dau_max)
  endif
  
  if (alap_taper_hi_lats) then
    alap_const = alap
    alap_taper_min = alap/alap_taper_min_fac
    do j=jsc,jec
      do i=isc,iec
        if (abs(Grid%yu(i,j)) > alap_taper_start_lat) then
           alap_const(i,j) = alap - &
            (alap-alap_taper_min)*(alap_taper_start_lat-abs(Grid%yu(i,j)))**2 &
            /(alap_taper_start_lat-alap_taper_end_lat)**2
        endif
        if (alap_const(i,j) < alap_taper_min) alap_const(i,j) = alap_taper_min
      enddo
    enddo
    call mpp_update_domains (alap_const, Dom%domain2d)
    do j=jsc,jec
      write (unit,'(a,i5,f10.3,f14.4)') ' lateral_friction_init: j, lat, alap_const =', &
                      j+Dom%joff, Grid%yu(isc,j), alap_const(isc,j)
    enddo
  endif

  allocate (visc_ceu(isd:ied,jsd:jed,nk)) ! east face of U cell
  allocate (visc_cnu(isd:ied,jsd:jed,nk)) ! north face of U cell
  allocate (visc_cu(isd:ied,jsd:jed,nk))  ! for metric term

  do k=1,nk
    visc_ceu(:,:,k) = BAY(alap_const(:,:))
    visc_cnu(:,:,k) = BAX(alap_const(:,:))
    visc_cu(:,:,k)  = BAY(BAX(alap_const(:,:)))
  enddo
  call mpp_update_domains (visc_ceu, Dom%domain2d)
  call mpp_update_domains (visc_cnu, Dom%domain2d)
  call mpp_update_domains (visc_cu,  Dom%domain2d)

  ! Munk boundary check
  if (dtime /= 0.0) then
    num = 0
    write (stdout(),'(/,(1x,a))')'==> Warning: locations (if any) where the Munk boundary layer is unresolved.'
    do j=jsc,jec
      do i=isc,iec
        cbeta = 2.28e-11
        alap_munk = cbeta*(Grd%h1u(i,j)/radius)*(Grd%dxu(i,j)*sqrt(3.0)/pi)**3
        if (visc_cu(i,j,1) <= alap_munk .and. num <= 10 .and. i == isc) then
          num = num + 1
          write (stdout(),'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') ' at (i,j)= (',i,',',j,'),(lon,lat)= (', &
                           Grd%xt(i,j),',',Grd%yt(i,j),&
                           '),  "alap" = ', visc_cu(i,j,1),&
                           ' m^2/s. the critical value =',alap_munk,' m^2/s'
        endif
      enddo
    enddo
    call mpp_sum (num)
    if (num > 10) write (stdout(),*) ' (a max of 10 violations per PE are listed)'
  endif

  ! for clocks 
  id_lap_friction  = mpp_clock_id('(Ocean constant lap friction) ' ,grain=CLOCK_MODULE)

  ! register diagnostics and send static viscosity 
  id_visc = register_static_field ('ocean_model', 'lap_visc', Grd%vel_axes_uv(1:2), 'lap viscosity', 'm^2/sec',&
                                 missing_value=-10.0, range=(/-10.0,1.e10/))
  if (id_visc > 0) used = send_data (id_visc, visc_cu(isc:iec,jsc:jec,1), &
                          Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))

  id_lap_fric_u = register_diag_field('ocean_model','lap_fric_u',Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz lap frict on u', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_lap_fric_v = register_diag_field('ocean_model','lap_fric_v',Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz lap frict on v', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))


end subroutine ocean_lap_friction_init
! </SUBROUTINE>  NAME="ocean_lap_friction_init"


!#######################################################################
! <FUNCTION NAME="horz_lap_friction">
!
! <DESCRIPTION>
! This function computes the thickness weighted time tendency for
! horizontal velocity from horizontal Laplacian friction.
! </DESCRIPTION>
!
function horz_lap_friction(Time, Thickness, Velocity, diag_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  logical, intent(in), optional          :: diag_flag 
  logical                                :: send_diagnostics

  real, dimension(isd:ied,jsd:jed)      :: metric
  real, dimension(isd:ied,jsd:jed,2)    :: uk
  real, dimension(isd:ied,jsd:jed)      :: pc_sink
  real, dimension(isd:ied,jsd:jed)      :: tmp, tmp1, tmp2, m1, m2, n1, n2
  real, dimension(isc:iec,jsc:jec,nk,2) :: horz_lap_friction
  integer :: i, j, k, n, taum1, tau

  if(.not. lap_friction_on) then 
    horz_lap_friction=0.0
    return 
  endif 

  call mpp_clock_begin(id_lap_friction)

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_friction): module needs to be initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  horz_lap_friction=0.0

  metric(:,:)  = 0.0
  pc_sink(:,:) = 0.0

  ! Murray and Reason Eqn 29

  do k=1,nk

    ! metric terms for generalized coordinate system
    m1(:,:) = 2.0*visc_cu(:,:,k)*Grd%dh1dy(:,:) + BAY(FDY_PU(visc_cu(:,:,k)))
    m2(:,:) = 2.0*visc_cu(:,:,k)*Grd%dh2dx(:,:) + BAX(FDX_PU(visc_cu(:,:,k)))
    tmp1(:,:) = visc_ceu(:,:,k)*Grd%dyue(:,:)
    tmp2(:,:) = visc_cnu(:,:,k)*Grd%dxun(:,:)
    n1(:,:) = - Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh2dx(:,:))) &
              - Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh1dy(:,:)))
    n2(:,:) = + Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh1dy(:,:))) &
              - Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh2dx(:,:)))

   ! Horz viscous bottom drag due to partial cells
    do j=jsc,jec
      do i=isc,iec              
        pc_sink(i,j) = - Grd%daur(i,j)*(  visc_ceu(i,j,k)*(&
                 (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i+1,j,k,tau),Thickness%dhu(i,j,k,tau)))*Grd%dyue_dxuer(i,j)  +    &
                 (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j,k,tau),Thickness%dhu(i-1,j,k,tau)))*Grd%dyue_dxuer(i-1,j)   ) &
               +                                     visc_cnu(i,j,k)*(&
                 (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j+1,k,tau),Thickness%dhu(i,j,k,tau)))*Grd%dxun_dyunr(i,j)  +    &
                 (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j,k,tau),Thickness%dhu(i,j-1,k,tau)))*Grd%dxun_dyunr(i,j-1) ) )
      enddo
    enddo

    do n=1,2
      uk(:,:,:) = Velocity%u(:,:,k,:,taum1)
      metric(:,:) = (3-2*n)*(m1(:,:)*FDX_NT(BAX(uk(:,:,3-n)))-&
                             m2(:,:)*FDY_ET(BAY(uk(:,:,3-n)))  )&
                          +  n1(:,:)*uk(:,:,n) + (3-2*n)*n2(:,:)*uk(:,:,3-n)
      tmp(:,:) = ( BDX_EU(visc_ceu(:,:,k)*FMX(Thickness%dhu(:,:,k,tau))*FDX_PU(uk(:,:,n)))&
                  +BDY_NU(visc_cnu(:,:,k)*FMY(Thickness%dhu(:,:,k,tau))*FDY_PU(uk(:,:,n))) ) &
                  +pc_sink(:,:)*uk(:,:,n) + metric(:,:)*Thickness%dhu(:,:,k,tau)
      horz_lap_friction(isc:iec,jsc:jec,k,n) = tmp(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,k)
    enddo

  enddo

  if(send_diagnostics) then
     if (id_lap_fric_u > 0) used = send_data(id_lap_fric_u, horz_lap_friction(isc:iec,jsc:jec,:,1), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_lap_fric_v > 0) used = send_data(id_lap_fric_v, horz_lap_friction(isc:iec,jsc:jec,:,2), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  call mpp_clock_end(id_lap_friction)

end function horz_lap_friction
! </FUNCTION> NAME="horz_lap_friction"



!#######################################################################
! <SUBROUTINE NAME="horz_lap_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the Laplacian 
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
subroutine horz_lap_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: crit, dsmin

  if(.not. lap_friction_on) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_viscosity_check): module needs to be initialized')
  endif 

  write (stdout(),'(/60x,a/)') ' Excessive horizontal Laplacian friction summary:'
  write (stdout(),'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then 
          dsmin = min(Grd%dxu(i,j),Grd%dxu(i+1,j),Grd%dyu(i,j),Grd%dyu(i,j+1))
          crit = 0.5*dsmin**2/max(dtime,epsln)
          if (visc_cu(i,j,k) > crit .and. num < max_num) then
            num = num + 1
            write (unit,9600) 'visc_cu(',visc_cu(i,j,k), crit, i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
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
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (horz_lap_reynolds_check): module needs to be initialized')
  endif 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(visc_cu(i,j,k) + epsln)
        ramn = 1.0/(visc_cu(i,j,k) + epsln)
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

end subroutine horz_lap_reynolds_check
! </SUBROUTINE> NAME="horz_lap_reynolds_check"




end module ocean_lap_friction_mod
      
      
