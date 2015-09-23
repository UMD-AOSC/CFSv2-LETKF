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
module ocean_adv_vel_diag_mod
!  
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Numerical diagnostics for advection velocity related quantities. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Numerical diagnostics for advection velocity related quantities.
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_adv_vel_diag_nml">
!  <DATA NAME="cflcrt" UNITS="dimensionless" TYPE="real">
!  Critical CFL number for CFL check.
!  </DATA> 
!  <DATA NAME="max_cfl_violations" UNITS="dimensionless" TYPE="integer">
!  Maximum number of CFL violations allowed before model is brought down.
!  </DATA> 
!  <DATA NAME="diag_freq" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!</NAMELIST>
!
use constants_mod,       only: epsln, rho0r 
use diag_manager_mod,    only: register_diag_field, send_data, need_data
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, stdout, stdlog
use mpp_mod,             only: mpp_error, mpp_max, mpp_pe
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,    only: time_type, increment_time

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: BAX, BAY, BDX_ET, BDY_NT, BDX_EU, BDY_NU
use ocean_operators_mod, only: REMAP_ET_TO_EU, REMAP_NT_TO_NU, REMAP_BT_TO_BU
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,     only: ocean_adv_vel_type, ocean_density_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_util_mod,      only: matrix
use ocean_workspace_mod, only: wrk1 

implicit none

private

#include <ocean_memory.h>

real    :: dtts
real    :: dtuv
real    :: dtime_t
real    :: dtime_u

integer :: index_temp
integer :: diag_freq=-1

! for output
integer :: unit=6

! for diagnostics clocks 
integer :: id_adv_vel_numerics
integer :: id_transport_on_z
integer :: id_transport_on_rho
integer :: id_transport_on_theta

! for CFL checks 
real    :: cflcrt=1.5           
integer :: numcfl=0             ! counter for number of times the "cflcrt" factor was exceeded.
integer :: max_cfl_violations=3 ! maximum number of times the "cflcrt" factor can be exceeded before terminating.

! for transport_on_z 
integer :: id_tx_trans=-1
integer :: id_ty_trans=-1
integer :: id_tz_trans=-1
integer :: id_ux_trans=-1
integer :: id_uy_trans=-1
integer :: id_uz_trans=-1

! for transport_on_rho
integer :: id_tx_trans_rho=-1
integer :: id_ty_trans_rho=-1

! for transport_on_theta
integer :: id_tx_trans_theta=-1
integer :: id_ty_trans_theta=-1

logical :: used

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

public ocean_adv_vel_diag_init
public ocean_adv_vel_diagnostics 

private remapping_check 
private cfl_check1
private cfl_check2
private maximum_bottom_w
private vertical_reynolds_check
public max_continuity_error

private transport_on_z
private transport_on_rho
private transport_on_theta

namelist /ocean_adv_vel_diag_nml/ cflcrt, max_cfl_violations, diag_freq


contains

!#######################################################################
! <SUBROUTINE NAME="ocean_adv_vel_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_adv_vel_diag module containing subroutines
! diagnosing advection velocity related properties of the simulation.
! </DESCRIPTION>
!
subroutine ocean_adv_vel_diag_init(Grid, Domain, Time, Time_steps, Adv_vel, T_prog, Dens)

type(ocean_grid_type), target, intent(in)   :: Grid
type(ocean_domain_type), target, intent(in) :: Domain
type(ocean_time_type), intent(in)           :: Time
type(ocean_time_steps_type), intent(in)     :: Time_steps
type(ocean_adv_vel_type), intent(in)        :: Adv_vel
type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
type(ocean_density_type), intent(in)        :: Dens

integer :: n, ioun, io_status, ierr

if (module_is_initialized) return

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

ioun = open_namelist_file()
read(ioun, ocean_adv_vel_diag_nml, iostat=io_status)
write (stdlog(), ocean_adv_vel_diag_nml)
write (stdout(),'(/)')
write (stdout(), ocean_adv_vel_diag_nml)
ierr = check_nml_error(io_status,'ocean_adv_vel_diag_nml')
call close_file(ioun)

if (diag_freq == 0) diag_freq = 1

#ifndef STATIC_MEMORY
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid

dtuv    = Time_steps%dtuv
dtts    = Time_steps%dtts
dtime_t = Time_steps%dtime_t
dtime_u = Time_steps%dtime_u

! set index for potential temperature 
do n=1,size(T_prog,1)
  if(T_prog(n)%name == 'temp' ) then 
    index_temp = n
  endif 
enddo 

! register fields for diagnostic output

id_tx_trans     = register_diag_field ('ocean_model','tx_trans', Grd%tracer_axes_flux_x(1:3), Time%model_time, &
  'T-cell x-volume transport','Sv', missing_value=-1e6, range=(/-1e6,1e6/))
id_ty_trans     = register_diag_field ('ocean_model','ty_trans', Grd%tracer_axes_flux_y(1:3), Time%model_time, &
  'T-cell y-volume transport','Sv', missing_value=-1e6, range=(/-1e6,1e6/))
id_tz_trans     = register_diag_field ('ocean_model','tz_trans', Grd%tracer_axes_wt(1:3), Time%model_time, &
  'T-cell vert volume transport','Sv', missing_value=-1e6, range=(/-1e6,1e6/))

id_ux_trans = register_diag_field ('ocean_model','ux_trans', Grd%vel_axes_flux_x(1:3), Time%model_time, &
  'U-cell x-volume transport','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))
id_uy_trans = register_diag_field ('ocean_model','uy_trans', Grd%vel_axes_flux_y(1:3), Time%model_time, &
  'U-cell y-volume transport','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))
id_uz_trans = register_diag_field ('ocean_model','uz_trans', Grd%vel_axes_wu(1:3), Time%model_time, &
  'U-cell vert volume transport','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))

id_tx_trans_rho = register_diag_field ('ocean_model','tx_trans_rho',  &
   Dens%potrho_axes(1:3), Time%model_time, &
  'T-cell x-volume transport on pot_rho','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))
id_ty_trans_rho = register_diag_field ('ocean_model','ty_trans_rho',  &
   Dens%potrho_axes(1:3), Time%model_time, &
  'T-cell y-volume transport on pot_rho','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))

id_tx_trans_theta = register_diag_field ('ocean_model','tx_trans_theta',  &
   Dens%theta_axes(1:3), Time%model_time, &
  'T-cell x-volume transport on theta','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))
id_ty_trans_theta = register_diag_field ('ocean_model','ty_trans_theta',  &
   Dens%theta_axes(1:3), Time%model_time, &
  'T-cell y-volume transport on theta','Sv', missing_value=-1.e6, range=(/-1e6,1e6/))

! set ids for clocks
id_adv_vel_numerics   = mpp_clock_id('(Ocean adv_vel_diag: numerics)'    ,grain=CLOCK_ROUTINE)
id_transport_on_z     = mpp_clock_id('(Ocean adv_vel_diag: z-trans)'     ,grain=CLOCK_ROUTINE)
id_transport_on_rho   = mpp_clock_id('(Ocean adv_vel_diag: rho-trans)'   ,grain=CLOCK_ROUTINE)
id_transport_on_theta = mpp_clock_id('(Ocean adv_vel_diag: theta-trans)' ,grain=CLOCK_ROUTINE)

call remapping_check

end subroutine ocean_adv_vel_diag_init
! </SUBROUTINE>  NAME="ocean_adv_vel_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_adv_vel_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the velocity. 
! </DESCRIPTION>

subroutine ocean_adv_vel_diagnostics(Time, Thickness, Adv_vel, T_prog, Dens, visc_cbu)

  type(ocean_time_type), intent(in)               :: Time
  type(ocean_thickness_type), intent(in)          :: Thickness
  type(ocean_adv_vel_type), intent(in)            :: Adv_vel
  type(ocean_prog_tracer_type), intent(in)        :: T_prog(:)
  type(ocean_density_type), intent(in)            :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: visc_cbu

  
  call mpp_clock_begin(id_adv_vel_numerics)
  if (diag_freq > 0) then
      if (mod(Time%itt,diag_freq) == 0) then
          call cfl_check1(Time, Thickness, Adv_vel)
          call cfl_check2(Time, Thickness, Adv_vel)
          call maximum_bottom_w(Adv_vel)
          call max_continuity_error(Adv_vel, Dens)
          call vertical_reynolds_check(Adv_vel, Thickness, visc_cbu)
      endif
  endif
  call mpp_clock_end(id_adv_vel_numerics)
  
  call mpp_clock_begin(id_transport_on_z)
    call transport_on_z(Time, Adv_vel)
  call mpp_clock_end(id_transport_on_z)

  call mpp_clock_begin(id_transport_on_rho)
    call transport_on_rho(Time, Dens, Adv_vel)
  call mpp_clock_end(id_transport_on_rho)
  
  call mpp_clock_begin(id_transport_on_theta)
    call transport_on_theta(Time, Dens, T_prog(index_temp), Adv_vel)
  call mpp_clock_end(id_transport_on_theta)


end subroutine ocean_adv_vel_diagnostics
! </SUBROUTINE>  NAME="ocean_adv_vel_diagnostics"



!#######################################################################
! <SUBROUTINE NAME="remapping_check">
!
! <DESCRIPTION>
! Compute remapping error.  This error will be roundoff only for model
! grids where the tracer and velocity grid cell distances are 
! linearly related.  The spherical version of mom4 satisfies the
! appropriate relation, and so should maintain roundoff for the 
! remapping error.  The tripolar version of mom4 does not have 
! tracer and velocity grids related linearly, and so the 
! "remapping error" is nontrivial.  The significance of this error
! is unclear.  No adverse effects have been identified.
! </DESCRIPTION>
!
subroutine remapping_check

  real, dimension(isd:ied,jsd:jed) :: uet, vnt, ueu, vnu, u, v, dw, err
  integer                          :: i, j, ierr, jerr
  real                             :: remap_err, remap_err0, fudge

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_adv_vel_diag_mod (remapping_check) : module needs initialization ')
  endif 

  u(:,:)   = 1.0*Grd%umask(:,:,1)
  v(:,:)   = 1.0*Grd%umask(:,:,1)
  uet(:,:) = BAY(u(:,:)*Grd%dyu(:,:))/Grd%dyte(:,:)
  vnt(:,:) = BAX(v(:,:)*Grd%dxu(:,:))/Grd%dxtn(:,:)
  dw(:,:)  = REMAP_BT_TO_BU(BDX_ET(uet(:,:)) + BDY_NT(vnt(:,:)))
  
  ueu(:,:) = REMAP_ET_TO_EU(uet(:,:))  
  vnu(:,:) = REMAP_NT_TO_NU(vnt(:,:))
  err(:,:) = abs(dw(:,:)-(BDX_EU(ueu(:,:))   + BDY_NU(vnu(:,:))))

  fudge = 1 + 1.e-12*mpp_pe()  ! used to distinguish processors when result is independent of processor
  remap_err=0.0; ierr=isc; jerr=jsc
  do j=jsc,jec
    do i=isc,iec
      if (abs(err(i,j)) > abs(remap_err)) then 
        remap_err  = err(i,j)
        ierr = i
        jerr = j
      endif
    enddo
  enddo
  remap_err  = remap_err*fudge
  remap_err0 = remap_err
  remap_err  = abs(remap_err)
  call mpp_max(remap_err)

  if (abs(remap_err0) == remap_err) then
    remap_err = remap_err0
    write (unit,9000)  remap_err/fudge, ierr+Dom%ioff, jerr+Dom%joff, Grd%xu(ierr,jerr), Grd%yu(ierr,jerr)
    if(Grd%tripolar) write (unit,'(a)') '==>Note: remapping error will be small (i.e., order 1e-20) only for spherical grids.'
  endif

9000  format(//'==>Maximum remapping error = ',es10.3,' m/s  at (i,j) = ','(',i4,',',i4,'),',' (lon,lat) = (',f7.2,',',f7.2,')')

end subroutine remapping_check
! </SUBROUTINE> NAME="remapping_check"


!#######################################################################
! <SUBROUTINE NAME="cfl_check1">
!
! <DESCRIPTION>
! Perform the first of two CFL checks.  
!
! Vectorized version from Russell.Fiedler@csiro.au computes cfl
! values at a single latitude. The location of the maximum at this
! latitude is calculated via the maxloc() intrinsic. The maximum 
! value for this processor is then updated if necessary.
!
! </DESCRIPTION>
subroutine cfl_check1(Time, Thickness, Adv_vel)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_adv_vel_type), intent(in)   :: Adv_vel

  real     :: cflup,cflvp,cflwtp    ! percent of cfl criteria reached by velocity component
  real     :: cflum,cflvm,cflwtm    ! velocity component which comes closest to its cfl criteria
  integer  :: icflu,jcflu,kcflu     ! "i","j","k" coordinate of "cflum"
  integer  :: icflv,jcflv,kcflv     ! "i","j","k" coordinate of "cflvm"
  integer  :: icflwt,jcflwt,kcflwt  ! "i","j","k" coordinate of "cflwtm"
  real     :: dtmax, umax, vmax, wmax
  real     :: cflup0, cflvp0, cflwtp0
  real     :: u, v, w, fudge
  integer  :: i, j, k, tau
 
  ! Work Array containing cfl values 
  real , dimension(isc:iec) :: tempcfl
  ! Array required for maxloc() intrinsic function
  integer, dimension(1)      :: itemp

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error in ocean_adv_vel_diag_mod (cfl_check1): module needs initialization ')
  endif 

  tau = Time%tau   

  write (stdout(),'(/60x,a/)') ' CFL summary I:  '
  write (stdout(),'(1x,a/)')'Location of nearest approach to local CFL velocity limit...'

  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  icflu=isc; jcflu=jsc; kcflu=1; cflup=epsln; cflum=0.0
  icflv=isc; jcflv=jsc; kcflv=1; cflvp=epsln; cflvm=0.0
  icflwt=isc; jcflwt=jsc; kcflwt=1; cflwtp=epsln; cflwtm=0.0
  dtmax = max(dtime_t, dtime_u)  

!!$  ! factors of 100.0 are to convert to percentages 
!!$
!!$  do k=1,nk
!!$    do j=jsc,jec
!!$      do i=isc,iec
!!$        wmax  = Thickness%dhwu(i,j,k)/dtmax            
!!$        vmax  = Grd%dytn(i,j)/dtmax
!!$        umax  = Grd%dxte(i,j)/dtmax
!!$        u     = Adv_vel%uh_et(i,j,k)/Thickness%dht(i,j,k,tau)
!!$        v     = Adv_vel%vh_nt(i,j,k)/Thickness%dht(i,j,k,tau)
!!$        w     = Adv_vel%w_bt(i,j,k)
!!$        if (abs(100.0*u/umax) > cflup) then
!!$          cflup = abs(100.0*u/umax)
!!$          cflum = u
!!$          icflu = i
!!$          jcflu = j
!!$          kcflu = k
!!$        endif
!!$        if (abs(100.0*v/vmax) > cflvp) then
!!$          cflvp = abs(100.0*v/vmax)
!!$          cflvm = v
!!$          icflv = i
!!$          jcflv = j
!!$          kcflv = k
!!$        endif
!!$        if (abs(100.0*Grd%tmask(i,j,k)*w/wmax) > cflwtp) then
!!$          cflwtp = abs(100.0*w/wmax)
!!$          cflwtm = w
!!$          icflwt = i
!!$          jcflwt = j
!!$          kcflwt = k
!!$        endif
!!$      enddo
!!$    enddo
!!$  enddo

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tempcfl(i)=abs(Adv_vel%uh_et(i,j,k)/(Thickness%dht(i,j,k,tau)*Grd%dxte(i,j)))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflup) then
            cflup=tempcfl(itemp(1))
            icflu=itemp(1)
            jcflu=j
            kcflu=k
        endif
        do i=isc,iec
           tempcfl(i)=abs(Adv_vel%vh_nt(i,j,k)/(Thickness%dht(i,j,k,tau)*Grd%dytn(i,j)))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(1)>cflvp) then
            cflvp=tempcfl(itemp(1))
            icflv=itemp(1)
            jcflv=j
            kcflv=k
        endif
        do i=isc,iec
           tempcfl(i)=abs(Grd%tmask(i,j,k)*Adv_vel%w_bt(i,j,k)/Thickness%dhwu(i,j,k))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflwtp) then
            cflwtp=tempcfl(itemp(1))
            icflwt=itemp(1)
            jcflwt=j
            kcflwt=k
        endif
     enddo
  enddo

  ! multiply by 100.0 to convert to percentages 
  ! multiply by dtmax to convert to dimensionless CFL number 
  cflup  = 100.0*cflup*dtmax
  cflvp  = 100.0*cflvp*dtmax
  cflwtp = 100.0*cflwtp*dtmax

  cflum  = Adv_vel%uh_et(icflu,jcflu,kcflu)/Thickness%dht(icflu,jcflu,kcflu,tau)
  cflvm  = Adv_vel%vh_nt(icflv,jcflv,kcflv)/Thickness%dht(icflv,jcflv,kcflv,tau)
  cflwtm = Adv_vel%w_bt(icflwt,jcflwt,kcflwt)

  cflup   = cflup*fudge
  cflvp   = cflvp*fudge
  cflwtp  = cflwtp*fudge
  cflup0  = cflup
  cflvp0  = cflvp
  cflwtp0 = cflwtp
  call mpp_max(cflup)
  call mpp_max(cflvp)
  call mpp_max(cflwtp)

  if (cflup == cflup0) then
    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f7.2,a,f7.2,a,f7.0,a)')&
    ' u_et (',cflum,' m/s) is ',cflup/fudge,' % of the local CFL limit (',abs(100.0*cflum/(cflup/fudge)),&
    ' m/s) at (i,j,k) = (',icflu+Dom%ioff,',',jcflu+Dom%joff,',',kcflu,'),'&
    ,' (lon,lat,dpt) = (',Grd%xu(icflu,jcflu),',',Grd%yu(icflu,jcflu),',',  Grd%zt(kcflu),' m)'
  endif

  if (cflvp == cflvp0) then
    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f7.2,a,f7.2,a,f7.0,a)')&
    ' v_nt (',cflvm,' m/s) is ',cflvp/fudge,' % of the local CFL limit (',abs(100.0*cflvm/(cflvp/fudge)),&
    ' m/s) at (i,j,k) = (',icflv+Dom%ioff,',',jcflv+Dom%joff,',',kcflv,'),'&
    ,' (lon,lat,dpt) = (',Grd%xu(icflv,jcflv),',',Grd%yu(icflv,jcflv),',',  Grd%zt(kcflv),' m)'
  endif

  if (cflwtp == cflwtp0) then
    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f7.2,a,f7.2,a,f7.0,a)')&
    ' w_bt (',cflwtm,' m/s) is ',cflwtp/fudge,' % of the local CFL limit (',abs(100.0*cflwtm/(cflwtp/fudge)),&
    ' m/s) at (i,j,k) = (',icflwt+Dom%ioff,',',jcflwt+Dom%joff,',',kcflwt,'),'&
    ,' (lon,lat,dpt) = (',Grd%xt(icflwt,jcflwt),',',Grd%yt(icflwt,jcflwt),',',  Grd%zw(kcflwt),' m)'
  endif


end subroutine cfl_check1
! </SUBROUTINE> NAME="cfl_check1"


!#######################################################################
! <SUBROUTINE NAME="cfl_check2">
!
! <DESCRIPTION>
! Perform the second of two CFL checks.  
! </DESCRIPTION>
subroutine cfl_check2(Time, Thickness, Adv_vel)

 type(ocean_time_type), intent(in)      :: Time
 type(ocean_thickness_type), intent(in) :: Thickness
 type(ocean_adv_vel_type), intent(in)   :: Adv_vel

  real :: dtmax, cflu, cflv, cflw
  real :: umax, pcflu, vmax, pcflv, wmax, pcflw, scl
  real, dimension(isd:ied,jsd:jed) :: u, v, w
  integer :: i, j, k, is, ie, js, je, numcfl, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (cfl_check2): module needs initialization ')
  endif 

  tau = Time%tau   

  numcfl = 0   ! counter for CFL violations 
  write (stdout(),'(/60x,a/)') ' CFL summary II:  '
  write (stdout(),'(1x,a,f5.2/)')'Locations (if any) where advective velocity exceeds CFL limit by a factor of', cflcrt
  dtmax = max(dtime_t, dtime_u)  
  do k=1,nk
    u(:,:) = Adv_vel%uh_et(:,:,k)/Thickness%dht(:,:,k,tau)
    v(:,:) = Adv_vel%vh_nt(:,:,k)/Thickness%dht(:,:,k,tau)
    w(:,:) = Adv_vel%w_bt(:,:,k)
    do j=jsc,jec
      do i=isc,iec
        cflu  = abs((dtmax/Grd%dxte(i,j))*u(i,j))    
        cflv  = abs((dtmax/Grd%dytn(i,j))*v(i,j))
        cflw  = abs((dtmax/Thickness%dhwu(i,j,k))*w(i,j))*Grd%tmask(i,j,k)
        if (cflu >= cflcrt .or. cflv >= cflcrt .or. cflw >= cflcrt) then
          write (stdout(),'(/,a,i4,a1,i3,a,i3,a,f6.3,/)') ' Note: CFL velocity limit exceeded at coordinate (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',cflcrt
          umax  = Grd%dxte(i,j)/dtmax            
          pcflu = abs(100.0*u(i,j)/umax)
          vmax  = Grd%dytn(i,j)/dtmax
          pcflv = abs(100.0*v(i,j)/vmax)
          wmax  = Thickness%dhwu(i,j,k)/dtmax
          pcflw = abs(100.0*w(i,j)/wmax)
          write (stdout(),'(a,f8.2,a,g15.8,a)') ' u reached   ', pcflu,' % of the local CFL limit (',umax,' m/s)'
          write (stdout(),'(a,f8.2,a,g15.8,a)') ' v reached   ', pcflv,' % of the local CFL limit (',vmax,' m/s)'
          write (stdout(),'(a,f8.2,a,g15.8,a)') ' w_bt reached', pcflw,' % of the local CFL limit (',wmax,' m/s)'
          is = max(isc,i-3)
          ie = min(iec,i+3)
          js = max(jsc,j-3)
          je = min(jec,j+3)
          scl = 0.0
          write (stdout(),9100)'u_et (m/s)', Time%itt, k, Grd%zt(k), Grd%xt(is,j), Grd%xt(ie,j), &
                               Grd%yt(i,js), Grd%yt(i,je), scl
          call matrix (u(is:ie,js:je), is, ie, js, je, scl)

          write (stdout(),9100) 'v_nt (m/s)', Time%itt, k, Grd%zt(k), Grd%xt(is,j), Grd%xt(ie,j), &
                               Grd%yt(i,js), Grd%yt(i,je), scl
          call matrix (v(is:ie,js:je), is, ie, js, je, scl)

          write (stdout(),9100)  'w_bt (m/s)', Time%itt, k, Grd%zw(k), Grd%xt(is,j), Grd%xt(ie,j), &
                               Grd%yt(i,js), Grd%yt(i,je), scl
          call matrix (w(is:ie,js:je), is, ie, js, je, scl)

          numcfl = numcfl + 1

          if (numcfl > max_cfl_violations) then
            call mpp_error(FATAL,'==>Error in ocean_adv_vel_diag_mod: max cfl violations exceeded... Terminating.')
          endif

        endif
      enddo
    enddo
  enddo

9100 format(1x,a12,1x,'ts=',i10,1x,',k=',i3,', z(m)=',f8.1,', lon(deg):',f6.2,' --> ',f6.2,&
            ', lat(deg):',f6.2,' --> ',f6.2,', scaling=',1pg10.3)

end subroutine cfl_check2
! </SUBROUTINE> NAME="cfl_check2"


!#######################################################################
! <SUBROUTINE NAME="maximum_bottom_w">
!
! <DESCRIPTION>
! Compute maximum vertical velocity on the bottom of tracer and velocity cells.
! The vertical velocity at bottom of a column of tracer cells should be roundoff.  
! For flat bottom simulations, the vertical velocity on the bottom of the 
! velocity cell column should also be roundoff.  For simulations with topography,
! the vertical velocity on the bottom of a velocity cell column will not vanish
! due to the effects of topography.  
! </DESCRIPTION>
!
subroutine maximum_bottom_w(Adv_vel)

  type(ocean_adv_vel_type), intent(in) :: Adv_vel

  real    :: wtbot, wubot, wtbot0, wubot0, fudge
  integer :: i, j, k, iwtbot, jwtbot, kwtbot, iwubot, jwubot, kwubot

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (maximum_bottom_w): module needs initialization ')
  endif 

  ! Find Max error in "w_bt" at bottom

  wtbot=0.0; iwtbot=isc; jwtbot=jsc; kwtbot=1
  do j=jsc,jec
    do i=isc,iec
      k = Grd%kmt(i,j)
      if (k /= 0 .and. (abs(Adv_vel%w_bt(i,j,k)) > abs(wtbot))) then
        wtbot  = Adv_vel%w_bt(i,j,k)
        iwtbot = i
        jwtbot = j
        kwtbot = k
      endif
    enddo
  enddo
  wtbot = abs(wtbot)

  ! Find Max "w_bu" at bottom (not an error: part of the slope velocity)

  wubot=0.0; iwubot=isc; jwubot=jsc; kwubot=1
  do j=jsc,jec
    do i=isc,iec
      k = Grd%kmu(i,j)
      if (k /= 0 .and. (abs(Adv_vel%w_bu(i,j,k)) > abs(wubot))) then 
        wubot  = Adv_vel%w_bu(i,j,k)
        iwubot = i
        jwubot = j
        kwubot = k
      endif
    enddo
  enddo
  wubot = abs(wubot)

  write (stdout(),'(//60x,a/)') ' Bottom vertical velocity summary:'
  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  wtbot  = wtbot*fudge
  wubot  = wubot*fudge
  wtbot0 = wtbot
  wubot0 = wubot
  wtbot  = abs(wtbot)
  wubot  = abs(wubot)
  call mpp_max(wtbot)
  call mpp_max(wubot)

  if (abs(wtbot0) == wtbot) then
    wtbot = wtbot0/fudge
    write (unit,9112) wtbot, iwtbot+Dom%ioff, jwtbot+Dom%joff, kwtbot, Grd%xt(iwtbot,jwtbot), Grd%yt(iwtbot,jwtbot), Grd%zt(kwtbot)
  endif
  
  if (abs(wubot0) == wubot) then
    wubot = wubot0/fudge
    write (unit,9113) wubot, iwubot+Dom%ioff, jwubot+Dom%joff, kwubot, Grd%xu(iwubot,jwubot), Grd%yu(iwubot,jwubot), Grd%zt(kwubot)
  endif

9112  format(/' Maximum T-cell bottom velocity (',es10.3,' m/s){error}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)')
9113  format(' Maximum U-cell bottom velocity (',es10.3,' m/s){slope}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)'/)

end subroutine maximum_bottom_w
! </SUBROUTINE> NAME="maximum_bottom_w"


!#######################################################################
! <SUBROUTINE NAME="max_continuity_error">
!
! <DESCRIPTION>
! Compute continuity error. Should be roundoff if all is working well.  
! </DESCRIPTION>
!
subroutine max_continuity_error(Adv_vel, Dens)

  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  type(ocean_density_type), intent(in) :: Dens
  
  real, dimension(isd:ied,jsd:jed,nk) :: divgt, divgu
  real                                :: bigt, bigu, bigt0, bigu0, fudge
  integer                             :: i, j, k, iu, ju, ku, it, jt, kt

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (max_continuity_error): module needs initialization ')
  endif 
  
  do k=1,nk
     divgt(:,:,k) = (BDX_ET(Adv_vel%uh_et(:,:,k)) + &
          BDY_NT(Adv_vel%vh_nt(:,:,k)) + &
          Adv_vel%w_bt(:,:,k-1)-Adv_vel%w_bt(:,:,k))*Grd%tmask(:,:,k)
     divgu(:,:,k) = (BDX_EU(Adv_vel%uh_eu(:,:,k)) + &
          BDY_NU(Adv_vel%vh_nu(:,:,k)) + &
          Adv_vel%w_bu(:,:,k-1)-Adv_vel%w_bu(:,:,k))*Grd%umask(:,:,k)
  enddo

  bigt = epsln !find location of largest continuity error
  bigu = epsln 
  it = isc; jt = jsc; kt = 1
  iu = isc; ju = jsc; ku = 1
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (abs(divgt(i,j,k)) > abs(bigt)) then
          bigt = divgt(i,j,k)
          it = i
          jt = j
          kt = k
        endif 
        if (abs(divgu(i,j,k)) > abs(bigu)) then
          bigu = divgu(i,j,k)
          iu = i
          ju = j
          ku = k
        endif 
      enddo
    enddo
  enddo
  
  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  bigt = bigt*fudge
  bigu = bigu*fudge
  bigt0 = bigt
  bigu0 = bigu
  bigt  = abs(bigt)
  bigu  = abs(bigu)
  call mpp_max(bigt)
  call mpp_max(bigu)

  write (stdout(),'(/60x,a/)') ' Continuity error summary:'

  if (abs(bigu0) == bigu) then
    bigu = bigu0/fudge
    write (unit,'(/,a,es10.3,a,i4,a1,i3,a1,i3,a,f7.2,a,f7.2,a,f7.0,a)') ' Maximum U-cell Continuity Error (',bigu,&
    ' m/s) is at (i,j,k) = (',iu+Dom%ioff,',',ju+Dom%joff,',',ku,'),  (lon,lat,dpt) = ('&
     ,Grd%xu(iu,ju),',',Grd%yu(iu,ju),',',  Grd%zt(ku),' m)'
  endif
  
  if (abs(bigt0) == bigt) then
    bigt = bigt0/fudge
    write (unit,'(a,es10.3,a,i4,a1,i3,a1,i3,a,f7.2,a,f7.2,a,f7.0,a/)') ' Maximum T-cell Continuity Error (',bigt,&
    ' m/s) is at (i,j,k) = (',it+Dom%ioff,',',jt+Dom%joff,',',kt,'),  (lon,lat,dpt) = ('&
     ,Grd%xt(it,jt),',',Grd%yt(it,jt),',',  Grd%zt(kt),' m)'
  endif

end subroutine max_continuity_error
! </SUBROUTINE> NAME="max_continuity_error"


!#######################################################################
! <SUBROUTINE NAME="transport_on_z">
!
! <DESCRIPTION>
! Compute volume transports on z-levels and send to diag_manager.
! </DESCRIPTION>
!
subroutine transport_on_z(Time, Adv_vel)

  type(ocean_time_type), intent(in)    :: Time
  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  integer :: i, j, k
  
  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (transport_on_z): module needs initialization ')
  endif 

  ! volume transports leaving faces of grid cells
  wrk1=0.0

  if (id_tx_trans > 0) then 
    do k=1,nk
         do j=jsc,jec 
            do i=isc,iec
               wrk1(i,j,k) = Adv_vel%uh_et(i,j,k)*Grd%dyte(i,j)*1e-6
            enddo
         enddo
    enddo 
    used = send_data (id_tx_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_ty_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%vh_nt(i,j,k)*Grd%dxtn(i,j)*1e-6 
         enddo 
      enddo
    enddo 
    used = send_data (id_ty_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_tz_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%w_bt(i,j,k)*Grd%dxt(i,j)*Grd%dyt(i,j)*1e-6 
         enddo 
      enddo
    enddo 
    used = send_data (id_tz_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_ux_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%uh_eu(i,j,k)*Grd%dyue(i,j)*1e-6
         enddo 
      enddo
    enddo 
    used = send_data (id_ux_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%umask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_uy_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%vh_nu(i,j,k)*Grd%dxun(i,j)*1e-6 
         enddo 
      enddo
    enddo 
    used = send_data (id_uy_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%umask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_uz_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%w_bu(i,j,k)*Grd%dxu(i,j)*Grd%dyu(i,j)*1e-6
         enddo 
      enddo 
    enddo 
    used = send_data (id_uz_trans, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%umask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

end subroutine transport_on_z
! </SUBROUTINE> NAME="transport_on_z"


!#######################################################################
! <SUBROUTINE NAME="transport_on_rho">
!
! <DESCRIPTION>
! Classify horizontal transport of volume according to potential 
! density classes. 
!
! Diagnostic makes sense when potrho is monotonically 
! increasing with depth, although the algorithm does 
! not explicitly make this assumption.  
!
! Stephen.Griffies@noaa.gov
! Zhi.Liang@noaa.gov
!
! </DESCRIPTION>
!
subroutine transport_on_rho (Time, Dens, Adv_vel)

  type(ocean_time_type), intent(in)    :: Time
  type(ocean_density_type), intent(in) :: Dens
  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  type(time_type)                      :: next_time 

  integer :: i, j, k, k_rho, potrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (transport_on_rho): module needs initialization ')
  endif 

  potrho_nk = size(Dens%potrho_ref(:))
  work(:,:,:,:) = 0.0

  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_tx_trans_rho,next_time) .or. need_data(id_ty_trans_rho,next_time)) then

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(Dens%potrho(i,j,k) <= Dens%potrho_ref(1)) then 
                   work(i,j,1,1) =  work(i,j,1,1) + Adv_vel%uh_et(i,j,k)
                   work(i,j,1,2) =  work(i,j,1,2) + Adv_vel%vh_nt(i,j,k)
               endif
               do k_rho=2,potrho_nk-1
                  if(     Dens%potrho_ref(k_rho) <  Dens%potrho(i,j,k)  ) then
                      if( Dens%potrho(i,j,k)     <= Dens%potrho_ref(k_rho+1)) then 
                          work(i,j,k_rho,1) =  work(i,j,k_rho,1) + Adv_vel%uh_et(i,j,k)
                          work(i,j,k_rho,2) =  work(i,j,k_rho,2) + Adv_vel%vh_nt(i,j,k)
                      endif
                  endif
               enddo
               if(Dens%potrho_ref(potrho_nk) < Dens%potrho(i,j,k)) then 
                   work(i,j,potrho_nk,1) =  work(i,j,potrho_nk,1) + Adv_vel%uh_et(i,j,k)             
                   work(i,j,potrho_nk,2) =  work(i,j,potrho_nk,2) + Adv_vel%vh_nt(i,j,k)
               endif
            enddo
         enddo
      enddo

      do k_rho=1,potrho_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_rho,1) = work(i,j,k_rho,1)*Grd%dyte(i,j)*1e-6*Grd%tmask(i,j,1)
               work(i,j,k_rho,2) = work(i,j,k_rho,2)*Grd%dxtn(i,j)*1e-6*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_rho > 0) then 
          used = send_data (id_tx_trans_rho, work(:,:,:,1), &
               Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)
      endif
      if (id_ty_trans_rho > 0) then 
          used = send_data (id_ty_trans_rho, work(:,:,:,2), &
               Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)
      endif

  endif

end subroutine transport_on_rho
! </SUBROUTINE> NAME="transport_on_rho"


!#######################################################################
! <SUBROUTINE NAME="transport_on_theta">
!
! <DESCRIPTION>
! Classify horizontal transport of volume according to potential 
! temperature classes.  This diagnostic is useful to deduce the 
! heat that is transported between potential temperature classes. 
!
! Diagnostic makes sense when theta is monotonically decreasing
! with depth, although the algorithm does not explicitly make
! this assumption.  
!
! Stephen.Griffies@noaa.gov
! Zhi.Liang@noaa.gov
!
! </DESCRIPTION>
!
subroutine transport_on_theta (Time, Dens, Theta, Adv_vel)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_density_type), intent(in)     :: Dens
  type(ocean_prog_tracer_type), intent(in) :: Theta
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_theta, theta_nk, tau
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (transport_on_theta): module needs initialization ')
  endif 

  tau = Time%tau

  theta_nk = size(Dens%theta_ref(:))
  work(:,:,:,:)  = 0.0

  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_tx_trans_theta,next_time) .or. need_data(id_ty_trans_theta,next_time)) then

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
              if(Theta%field(i,j,k,tau) <= Dens%theta_ref(1)) then 
                 work(i,j,1,1) =  work(i,j,1,1) + Adv_vel%uh_et(i,j,k)
                 work(i,j,1,2) =  work(i,j,1,2) + Adv_vel%vh_nt(i,j,k)
                  endif
      do k_theta=2,theta_nk-1
                  if(    Dens%theta_ref(k_theta) <  Theta%field(i,j,k,tau) ) then
                     if( Theta%field(i,j,k,tau)  <= Dens%theta_ref(k_theta+1)) then 
                         work(i,j,k_theta,1) =  work(i,j,k_theta,1) + Adv_vel%uh_et(i,j,k)
                         work(i,j,k_theta,2) =  work(i,j,k_theta,2) + Adv_vel%vh_nt(i,j,k)
                     endif
                  endif
               enddo
              if(Dens%theta_ref(theta_nk) < Theta%field(i,j,k,tau)) then 
                 work(i,j,theta_nk,1) =  work(i,j,theta_nk,1) + Adv_vel%uh_et(i,j,k)
                 work(i,j,theta_nk,2) =  work(i,j,theta_nk,2) + Adv_vel%vh_nt(i,j,k)
                  endif
               enddo
            enddo
         enddo

      do k_theta=1,theta_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_theta,1) = work(i,j,k_theta,1)*Grd%dyte(i,j)*1e-6*Grd%tmask(i,j,1)
               work(i,j,k_theta,2) = work(i,j,k_theta,2)*Grd%dxtn(i,j)*1e-6*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_theta > 0) then 
          used = send_data (id_tx_trans_theta, work(:,:,:,1), &
                 Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)
      endif
      if (id_ty_trans_theta > 0) then 
          used = send_data (id_ty_trans_theta, work(:,:,:,2), &
                 Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)
      endif

  endif

end subroutine transport_on_theta
! </SUBROUTINE> NAME="transport_on_theta"


!#######################################################################
! <SUBROUTINE NAME="vertical_reynolds_check">
!
! <DESCRIPTION>
! This subroutine computes the Reynolds number associated with vertical
! friction visc_cbu and vertical velocity w_bu.
! </DESCRIPTION>
!
subroutine vertical_reynolds_check(Adv_vel, Thickness, visc_cbu)

  type(ocean_adv_vel_type), intent(in)            :: Adv_vel
  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: visc_cbu

  real    :: ramb, reyz
  real    :: fudge
  real    :: reynz0, reynz, reynw, reynmw
  integer :: ireynz, jreynz, kreynz
  integer :: i, j, k, kk

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error ocean_adv_vel_diag_mod (vertical_reynolds_check): module needs initialization')
  endif 

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynz=isc; jreynz=jsc; kreynz=1; reynz=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        ramb = 1.0/(visc_cbu(i,j,k) + epsln)

        kk = min(k+1,nk)
        if (k >= Grd%kmu(i,j)) then
          reyz = 0.0
        else 
          reyz = Grd%umask(i,j,kk)*abs(Adv_vel%w_bu(i,j,k)*Thickness%dhwu(i,j,k))*ramb
        endif
        if (reyz > reynz) then
          ireynz = i
          jreynz = j
          kreynz = k
          reynz  = reyz
          reynw  = Adv_vel%w_bu(i,j,k)
          reynmw = 1.0/ramb
        endif
      enddo
    enddo
  enddo
  write (stdout(),'(/60x,a/)') ' Vertical Laplacian Grid Reynolds number summary'


  fudge = 1 + 1.e-12*mpp_pe()  ! used to distinguish processors when result is independent of processor
  reynz = reynz*fudge 
  if(reynz == 0.0) then 
    reynz = reynz + 1.e-12*mpp_pe() 
  endif 

  reynz0 = reynz
  call mpp_max(reynz)
  
  if (reynz == reynz0) then
     write(unit,'(/,a,es9.2,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f8.2,a,/,1x,a,es9.2,a,es9.2,/)') &
           ' Maximum vertical grid Re # = ', reynz,&
           ' at (i,j,k) = (',ireynz+Dom%ioff,',',jreynz+Dom%joff,',',kreynz,&
           '),  (lon,lat,dpt) = (',Grd%xu(ireynz,jreynz),',', &
           Grd%yu(ireynz,jreynz),',', Grd%dzt(kreynz),' m). ', 'Vertical velocity "w_bu" is (m/s) ',&
            reynw, ' and vertical viscosity (m^2/s) "visc_cbu" is ', reynmw
  endif

end subroutine vertical_reynolds_check
! </SUBROUTINE> NAME="vertical_reynolds_check"


end module ocean_adv_vel_diag_mod


