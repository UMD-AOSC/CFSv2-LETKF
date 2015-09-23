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
module ocean_bih_friction_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! Thickness weighted time tendency for velocity from horizontal biharmonic friction 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for 
! horizontal velocity arising from horizontal biharmonic friction. 
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
! The numerical implementation requires one call to mpp_update_domains if 
! running the model with halo=1.  
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
!<NAMELIST NAME="ocean_bih_friction_const_nml">
!  <DATA NAME="bih_friction_on" TYPE="logical">
!  Must be true to use this module.
!  </DATA> 
!  <DATA NAME="abih" UNITS="m^4/sec" TYPE="real">
!  This is the value for the space-time constant biharmonic viscosity. 
!  </DATA> 
!  <DATA NAME="velocity_mix_micom" TYPE="logical">
!  If .true., then the viscosity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity within a latitude
!  band surrounding the equator.  This is useful for some models of enhanced equatorial
!  resolution that can maintain numerical integrity in this region with less friction 
!  than outside the tropical band.  
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity is set according 
!  to eq_vel_micom.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: mpp_error, mpp_sum, mpp_pe, mpp_max, mpp_sum
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: FAX, FAY, BAX, BAY, FMX, FMY
use ocean_operators_mod, only: BDX_EU, BDY_NU, FDX_PU, FDY_PU, FDX_NT, FDY_ET
use ocean_types_mod,     only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_velocity_type

implicit none

public ocean_bih_friction_init
public horz_bih_friction
public horz_bih_viscosity_check
public horz_bih_reynolds_check
private delsq_velocity

private

! time step 
real ::  dtime = 0.0  ! forward time step for friction tendency (2*dtuv if threelevel, dtuv if twolevel) 

! clock ids
integer :: id_bih_friction

! for diagnostics 
integer :: id_aiso=-1
integer :: id_bih_fric_u=-1
integer :: id_bih_fric_v=-1
logical :: used

real, dimension(:,:,:,:), allocatable :: del2_vel ! del^2 of horizontal velocity
real, dimension(:,:,:),   allocatable :: visc_ceu ! viscosity on east face of U cells (m^4/s)
real, dimension(:,:,:),   allocatable :: visc_cnu ! viscosity on north face of U cells (m^4/s)
real, dimension(:,:,:),   allocatable :: visc_cu  ! variable viscosity averaged to U cell (m^4/s)

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
  '=>Using:const/ocean_bih_friction.f90 ($Id$)'

character (len=128) :: tagname = &
  '$Name$'

logical :: module_is_initialized = .FALSE.
logical :: bih_friction_on       = .false.
real ::    abih                  = 0.0    ! constant horz biharmonic viscosity (m^4/s)
real ::    vel_micom             = 0.0    ! background scaling velocity for micom viscosity (m/sec) 
real ::    eq_vel_micom          = 0.0    ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real ::    eq_lat_micom          = 0.0    ! equatorial latitude band for micom (degrees)
logical :: velocity_mix_micom    =.false. ! if true, diffusivity made a function of the grid spacing

namelist /ocean_bih_friction_const_nml/ bih_friction_on, abih, velocity_mix_micom, &
                                        vel_micom, eq_vel_micom, eq_lat_micom

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_friction_init">
!
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
  integer :: i, j, k, num
  real    :: cbeta, alap_munk
  real, dimension(Domain%isd:Domain%ied,Domain%jsd:Domain%jed) :: Biso

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_friction_mod (ocean_bih_friction_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_bih_friction_const_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_bih_friction_const_nml)  
  write (stdlog(),ocean_bih_friction_const_nml)
  ierr = check_nml_error(io_status,'ocean_bih_friction_const_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(bih_friction_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_bih_friction_mod from bih/const.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_bih_friction_mod.')
    return
  endif 

  Grd => Grid
  Dom => Domain

  dtime = d_time

  allocate (visc_ceu(isd:ied,jsd:jed,nk))   ! east face of U cell
  allocate (visc_cnu(isd:ied,jsd:jed,nk))   ! north face of U cell
  allocate (visc_cu(isd:ied,jsd:jed,nk))    ! for metric term
  allocate (del2_vel(isd:ied,jsd:jed,nk,2)) ! laplacian of velocity on U cell

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_bih_friction_mod: using forward time step of (secs)', dtime 

  Biso(:,:) = 0.0
  if(velocity_mix_micom) then 
     do j=jsd,jed
        do i=isd,ied
           Biso(i,j) = vel_micom*((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           if(abs(Grd%yu(i,j)) < eq_lat_micom) then
             Biso(i,j) = eq_vel_micom*((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           endif
       enddo  
     enddo
  else
    Biso(:,:) = abs(abih)
  endif 

  do k=1,nk
    visc_ceu(:,:,k) = BAY(Biso(:,:))  
    visc_cnu(:,:,k) = BAX(Biso(:,:))  
    visc_cu(:,:,k)  = BAY(BAX(Biso(:,:)))          
  enddo
  call mpp_update_domains (visc_ceu, Dom%domain2d)
  call mpp_update_domains (visc_cnu, Dom%domain2d)
  call mpp_update_domains (visc_cu,  Dom%domain2d)

  ! for clocks 
  id_bih_friction  = mpp_clock_id('(Ocean constant bih friction) ' ,grain=CLOCK_MODULE)

  ! diagnostics
  id_aiso = register_static_field ('ocean_model', 'aiso_bih', Grd%vel_axes_uv(1:2), 'bih viscosity', 'm^4/sec',&
                                 missing_value=-10.0, range=(/-10.0,1.e20/))
  if (id_aiso > 0) then 
    used = send_data (id_aiso, visc_cu(isc:iec,jsc:jec,1), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))
  endif 

  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz bih-frict on u', 'm^2/s^2', missing_value=-10.0, range=(/-10.0,1.e10/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz bih-frict on v','m^2/s^2', missing_value=-10.0, range=(/-10.0,1.e10/))


end subroutine ocean_bih_friction_init
! </SUBROUTINE>  NAME="ocean_bih_friction_init"


!#######################################################################
! <FUNCTION NAME="horz_bih_friction">
!
! <DESCRIPTION>
! This function computes the thickness weighted acceleration on
! horizontal velocity arising from horizontal biharmonic friction.
! </DESCRIPTION>
!
function horz_bih_friction(Time, Thickness, Velocity, diag_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  logical, intent(in), optional          :: diag_flag 
  logical                                :: send_diagnostics

  real, dimension(isd:ied,jsd:jed)      :: metric
  real, dimension(isd:ied,jsd:jed,2)    :: uk
  real, dimension(isd:ied,jsd:jed)      :: pc_sink
  real, dimension(isd:ied,jsd:jed)      :: tmp, tmp1, tmp2, m1, m2, n1, n2
  real, dimension(isc:iec,jsc:jec,nk,2) :: horz_bih_friction
  integer :: i, j, k, n, tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_friction_mod (horz_bih_friction): module needs to be initialized')
  endif 

  if(.not. bih_friction_on) then 
    horz_bih_friction=0.0
    return 
  endif 

  call mpp_clock_begin(id_bih_friction)

  tau = Time%tau 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif

  metric(:,:)  = 0.0
  pc_sink(:,:) = 0.0

  call delsq_velocity (Time, Thickness, Velocity)
  if(Dom%xhalo==1 .or. Dom%yhalo==1 ) call mpp_update_domains (del2_vel(:,:,:,:), Dom%domain2d)

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

   ! Bih viscous bottom drag due to partial cells
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
      uk(:,:,:) = -del2_vel(:,:,k,:)
      metric(:,:) = (3-2*n)*(m1(:,:)*FDX_NT(BAX(uk(:,:,3-n)))-&
                             m2(:,:)*FDY_ET(BAY(uk(:,:,3-n)))  )&
                          +  n1(:,:)*uk(:,:,n) + (3-2*n)*n2(:,:)*uk(:,:,3-n)
      tmp(:,:) = ( BDX_EU(visc_ceu(:,:,k)*FMX(Thickness%dhu(:,:,k,tau))*FDX_PU(uk(:,:,n)))&
                  +BDY_NU(visc_cnu(:,:,k)*FMY(Thickness%dhu(:,:,k,tau))*FDY_PU(uk(:,:,n))) ) &
                  +pc_sink(:,:)*uk(:,:,n) + metric(:,:)*Thickness%dhu(:,:,k,tau)
      horz_bih_friction(isc:iec,jsc:jec,k,n) = tmp(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,k)
    enddo
  enddo

  if(send_diagnostics) then
     if (id_bih_fric_u > 0) used = send_data(id_bih_fric_u, horz_bih_friction(isc:iec,jsc:jec,:,1), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_bih_fric_v > 0) used = send_data(id_bih_fric_v, horz_bih_friction(isc:iec,jsc:jec,:,2), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  call mpp_clock_end(id_bih_friction)

end function horz_bih_friction
! </FUNCTION> NAME="horz_bih_friction"


!#######################################################################
! <SUBROUTINE NAME="delsq_velocity">
!
! <DESCRIPTION>
! Subroutine computes the laplacian operator acting on velocity.
! Resulting array del2_vel has dimensions [m^-1 * sec^-1]
! </DESCRIPTION>
!
subroutine delsq_velocity (Time, Thickness, Velocity)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity

  real, dimension(isd:ied,jsd:jed)       :: fe, fn, sink_pc, metric, n1, n2
  integer                                :: i, j, k, n, taum1, tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_friction_mod (delsq_velocity): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  sink_pc(:,:) = 0.0
  metric(:,:)  = 0.0

  ! metric terms for generalized system
  n1(:,:) = - Grd%dyur(:,:)*BDX_EU(Grd%dyue(:,:)*FAX(Grd%dh2dx(:,:)))      &
            - Grd%dxur(:,:)*BDY_NU(Grd%dxun(:,:)*FAY(Grd%dh1dy(:,:)))
  n2(:,:) = + Grd%dyur(:,:)*BDX_EU(Grd%dyue(:,:)*FAX(Grd%dh1dy(:,:)))      &
            - Grd%dxur(:,:)*BDY_NU(Grd%dxun(:,:)*FAY(Grd%dh2dx(:,:)))
  do k=1,nk

    ! Sink due to partial cells
    do j=jsc,jec
      do i=isc,iec
        sink_pc(i,j) = -(Grd%daur(i,j)/Thickness%dhu(i,j,k,tau))*(&
            (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i+1,j,k,tau),Thickness%dhu(i,j,k,tau)))  *Grd%dyue_dxuer(i,j)  +&
            (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j,k,tau),  Thickness%dhu(i-1,j,k,tau)))*Grd%dyue_dxuer(i-1,j)+&
            (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j+1,k,tau),Thickness%dhu(i,j,k,tau)))  *Grd%dxun_dyunr(i,j)  +&
            (Thickness%dhu(i,j,k,tau)-min(Thickness%dhu(i,j,k,tau),  Thickness%dhu(i,j-1,k,tau)))*Grd%dxun_dyunr(i,j-1))
      enddo
    enddo

    do n=1,2

      ! Viscous flux across east and north face of U cells
      fe(:,:) = FDX_PU(Velocity%u(:,:,k,n,taum1))*FMX(Thickness%dhu(:,:,k,tau))
      fn(:,:) = FDY_PU(Velocity%u(:,:,k,n,taum1))*FMY(Thickness%dhu(:,:,k,tau))

      metric(:,:) = (3-2*n)*2.0*(Grd%dh1dy(:,:)*FDX_NT(BAX(Velocity%u(:,:,k,3-n,taum1)))-&
                                Grd%dh2dx(:,:)*FDY_ET(BAY(Velocity%u(:,:,k,3-n,taum1)))  )& 
       + n1(:,:)*Velocity%u(:,:,k,n,taum1) + (3-2*n)*n2(:,:)*Velocity%u(:,:,k,3-n,taum1)

      del2_vel(:,:,k,n) = ( (BDX_EU(fe(:,:)) + BDY_NU(fn(:,:)))/Thickness%dhu(:,:,k,tau)&
                        + metric(:,:) + sink_pc(:,:)*Velocity%u(:,:,k,n,taum1))*Grd%umask(:,:,k) 
    enddo
  enddo 

end subroutine delsq_velocity
! </SUBROUTINE> NAME="delsq_tracer"


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
    call mpp_error(FATAL, '==>Error in ocean_bih_friction_mod (horz_bih_viscosity_check): module needs to be initialized')
  endif 

  write (stdout(),'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdout(),'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then            
          dsmin = min(Grd%dxu(i,j),Grd%dxu(i+1,j),Grd%dyu(i,j),Grd%dyu(i,j+1))
          crit = 0.0625*dsmin**4/max(dtime,epsln)
          if (visc_cu(i,j,k) > crit .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'visc_cu(',visc_cu(i,j,k), crit, i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdout(),*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es10.3,' m^4/s) exceeds max value (',es10.3,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
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

  if(.not. bih_friction_on) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_friction_mod (horz_bih_reynolds_check): module needs to be initialized')
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
      
      
