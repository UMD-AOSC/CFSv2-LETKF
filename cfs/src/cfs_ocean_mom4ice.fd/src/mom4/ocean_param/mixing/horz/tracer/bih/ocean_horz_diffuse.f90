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
module ocean_horz_diffuse_mod
! 
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! Thickness weighted time tendency for tracer 
! from horizontal biharmonic diffusion  
!</OVERVIEW>
!
!<DESCRIPTION>
! The diffusivity  used to determine the strength of the tendency can be 
! a general function of space yet it is constant in time.  A namelist 
! option exists that determines this diffusivity as a local function 
! of the grid spacing. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_horz_diffuse_bih_nml">
!  <DATA NAME="horz_diffuse_on" TYPE="logical">
!  Must be true to use this module
!  </DATA> 
!  <DATA NAME="abih" UNITS="m^4/sec" TYPE="real">
!  This is the value for the space-time constant biharmonic diffusivity. 
!  </DATA> 
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the diffusivity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: FATAL, NOTE, stdout, stdlog
use fms_mod,          only: write_version_number, open_namelist_file, close_file, check_nml_error
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE
use mpp_mod,          only: mpp_error

use ocean_domains_mod,   only: get_local_indices, set_ocean_domain
use ocean_operators_mod, only: FMX, FMY, FDX_PT, FDY_PT, BDX_ET, BDY_NT
use ocean_obc_mod,       only: store_ocean_obc_tracer_flux
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type, ocean_time_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_prog_tracer_type
use ocean_types_mod,     only: missing_value
use ocean_workspace_mod, only: wrk1 


implicit none

private

public ocean_horz_diffuse_init
public horz_diffuse
private del2_tracer

real ::    abih      = 0.0e10        ! constant horz biharmonic diffusion coeff for tracers
real ::    vel_micom = 0.0           ! constant velocity scale (m/s) for setting micom diffusivity  
logical :: tracer_mix_micom =.false. ! for spatial dependent diffusivity 

real, dimension(:,:,:), allocatable :: diff_cet    ! horizontal diffusivity for eastern face of T cell (m^2/s)
real, dimension(:,:,:), allocatable :: diff_cnt    ! horizontal diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:,:), allocatable :: del2_tracer ! Laplacian operator acting on tracers

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xflux_diff
integer, dimension(:), allocatable  :: id_yflux_diff
integer, dimension(:), allocatable  :: id_h_diffuse

character(len=128) :: version=&
     '=>Using: /bih/ocean_horz_diffuse.F90 ($Id$)'
character (len=128) :: tagname = &
     '$Name$'

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

real, dimension(:,:,:), allocatable :: fx, fy
integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: num_prog_tracers      = 0
logical :: horz_diffuse_on       = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc=.false.

namelist /ocean_horz_diffuse_bih_nml/ horz_diffuse_on, abih, tracer_mix_micom, vel_micom

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_horz_diffuse_init">
!
! <DESCRIPTION>
! Initialize the horizontal biharmonic diffusion module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that diffusivity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_horz_diffuse_init(Grid, Domain, Time, T_prog, dtime, obc, tmask_sigma)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, intent(in)                            :: dtime
  logical, intent(in)                         :: obc
  real, intent(in), optional                  :: tmask_sigma(Domain%isc:Domain%iec,Domain%jsc:Domain%jec)

  real    :: dxdymn, ah_crit
  integer :: ioun, io_status, ierr
  integer :: i, j, k, n, num
 
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_horz_diffuse_mod (ocean_horz_diffuse_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  have_obc = obc

  num_prog_tracers = size(T_prog(:))

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun  = open_namelist_file()
  read (ioun,ocean_horz_diffuse_bih_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_horz_diffuse_bih_nml)  
  write (stdlog(),ocean_horz_diffuse_bih_nml)
  ierr = check_nml_error(io_status, 'ocean_horz_diffuse_bih_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(horz_diffuse_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_horz_diffuse_mod from bih.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_horz_diffuse_mod.')
    return
  endif 

  if(abih > 0 .or. (tracer_mix_micom .and. vel_micom > 0)) then 
     write (stdout(),'(/1x,a)') ' ==> NOTE: USING horizontal biharmonic tracer diffusion.'
     write (stdout(),'(1x,a)') ' ==>WARNING: horizontal biharmonic diffusion will generally contribute to unphysically'
     write (stdout(),'(7x,a)') 'large diapycnal mixing. Such may be damaging for longer simulations.'
     write (stdout(),'(7x,a)') 'It is recommended that one instead uses ocean_neutral_physics.F90'
  endif

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_horz_diffuse_mod: using forward time step of (secs)', dtime 

  Dom => Domain
  Grd => Grid
  
  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux')

  allocate (diff_cet(isd:ied,jsd:jed,nk))
  allocate (diff_cnt(isd:ied,jsd:jed,nk))
  allocate (del2_tracer(isd:ied,jsd:jed,nk))
  allocate (fx(isd:ied,jsd:jed,nk))
  allocate (fy(isd:ied,jsd:jed,nk))  
  
  diff_cet    = abih
  diff_cnt    = abih
  del2_tracer = 0.0
  fx          = 0.0
  fy          = 0.0
  
  ! Micom background diffusivity
  ! Space scale set by grid size
  ! Velocity scale input via namelist 
  if(tracer_mix_micom) then
    do k=1,nk
      diff_cet(:,:,k) = vel_micom*(2.0*Grd%dxt(:,:)*Grd%dyte(:,:)/(Grd%dxt(:,:)+Grd%dyte(:,:)))**3
      diff_cnt(:,:,k) = vel_micom*(2.0*Grd%dxtn(:,:)*Grd%dyt(:,:)/(Grd%dxtn(:,:)+Grd%dyt(:,:)))**3
    enddo
    do j=jsc,jec
      write (stdout(),'(a,i4,a,e14.7,a)') ' MICOM biharmonic diffusivity at (isc,',j,',1) = ',diff_cet(isc,j,1),' m^4/s'
    enddo
  endif 

  if(PRESENT(tmask_sigma)) then 
     do j=jsc,jec
        do i=isc,iec
           if(tmask_sigma(i,j) > 0.0) then   
              k = Grd%kmt(i,j)
              diff_cet(i,j,k) = diff_cet(i,j,k)*(1.0-tmask_sigma(i,j))
              diff_cnt(i,j,k) = diff_cnt(i,j,k)*(1.0-tmask_sigma(i,j))
           endif
        enddo
     enddo
  endif
  
  if(PRESENT(tmask_sigma) .or. tracer_mix_micom) then
    call mpp_update_domains (diff_cet(:,:,:), Dom%domain2d)
    call mpp_update_domains (diff_cnt(:,:,:), Dom%domain2d)
  endif

  if (dtime /= 0.0) then
    num = 0
    do j=jsc,jec
      do i=isc,iec
        dxdymn = (2.0/(1.0/(Grd%dxt(i,j))**2 + 1.0/Grd%dyt(i,j)**2))**2
        ah_crit = 0.0625*dxdymn/dtime
        if (diff_cet(i,j,1)*Grd%tmask(i,j,1) > ah_crit .and. num <= 10 .and. i == isc) then
          num = num + 1
          if (num == 1) write (stdout(),'(/,(1x,a))')&
          '==> Warning: horz diffusive criteria exceeded for "ah". use a smaller "dtts", "ah". Show first 10 violations' 
          write (stdout(),'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') ' at (i,j)= (',i,',',j,'), (lon,lat)= (', &
                Grd%xt(i,j),',',Grd%yt(i,j), '),  "ah" = ', diff_cet(i,j,1),' m^4/s. the critical value =',ah_crit,' m^4/s'
        endif
      enddo
    enddo
  endif

  ! register for diagnostics manager 
  allocate (id_xflux_diff(num_prog_tracers))
  allocate (id_yflux_diff(num_prog_tracers))
  allocate (id_h_diffuse(num_prog_tracers))
  id_xflux_diff=-1
  id_yflux_diff=-1
  id_h_diffuse =-1

  do n=1, num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_xflux_diff(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd i-diffusive heat flux', 'Watts/m',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_yflux_diff(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd j-diffusive heat flux', 'Watts/m',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd horz-diffusion heating', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
     else
         id_xflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd i-diffusive flux of '//trim(T_prog(n)%name), 'm^2*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))
         id_yflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd j-diffusive flux of '//trim(T_prog(n)%name), 'm^2*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd horz-diffusion of '//trim(T_prog(n)%name), 'm*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
     endif

  enddo

end subroutine ocean_horz_diffuse_init
! </SUBROUTINE>  NAME="ocean_horz_diffuse_init"


!#######################################################################
! <SUBROUTINE NAME="horz_diffuse">
!
! <DESCRIPTION>
! This function computes the thickness weighted time tendency for tracer 
! from horizontal biharmonic diffusion. 
! </DESCRIPTION>
!
subroutine horz_diffuse (Time, Thickness, Tracer, ntracer, diag_flag)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer, intent(in)                         :: ntracer
  logical, intent(in), optional               :: diag_flag 
  logical                                     :: send_diagnostics 

  real, dimension(isd:ied,jsd:jed)    :: tchg
  real, dimension(isd:ied,jsd:jed,2)  :: tracr
  integer                             :: i, j, k
  integer                             :: tau, taum1
 
  if(.not. horz_diffuse_on) then 
    return 
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! fx = flux component through "eastern"  face of T-cells at level k
  ! fy = flux component through "northern" face of T-cells at level k
  ! tracr(:,:,1) = tracer at level k-1
  ! tracr(:,:,2) = tracer at level k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_horz_diffuse_mod (horz_diffuse): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  call delsq_tracer (Thickness, Tracer%field(:,:,:,taum1), tau)
  if(Dom%xhalo==1 .and. Dom%yhalo==1) call mpp_update_domains (del2_tracer, Dom%domain2d)
  tracr(:,:,1) = -del2_tracer(:,:,1)
  do k=1,nk
     tracr(:,:,2) = -del2_tracer(:,:,k)
     fx(:,:,k)    = diff_cet(:,:,k)*FDX_PT(tracr(:,:,1:2),k)*FMX(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
     fy(:,:,k)    = diff_cnt(:,:,k)*FDY_PT(tracr(:,:,1:2),k)*FMY(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
     tracr(:,:,1) = tracr(:,:,2)
  enddo

  !  set redundancies for tripolar
  if(Grd%tripolar) call mpp_update_domains(fx, fy, Dom_flux%domain2d, gridtype=CGRID_NE)

  do k=1,nk
     wrk1(:,:,k) = (BDX_ET(fx(:,:,k)) + BDY_NT(fy(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Grd%tmask(i,j,k)*wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! send fluxes to diag_manager 
  ! absence of minus sign relative to laplacian due to minus sign used for biharmonic.
  if(send_diagnostics) then 

      if (id_xflux_diff(ntracer) > 0) then
          used = send_data(id_xflux_diff(ntracer), fx(:,:,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)           
      endif
      if (id_yflux_diff(ntracer) > 0) then
          used = send_data(id_yflux_diff(ntracer), fy(:,:,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)           
      endif
      if (id_h_diffuse(ntracer) > 0) then
          used = send_data(id_h_diffuse(ntracer), wrk1(isc:iec,jsc:jec,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))           
      endif

  endif

  if (have_obc) call store_ocean_obc_tracer_flux(Tracer,-1.0*fx, -1.0*fy)
  
end subroutine horz_diffuse
! </SUBROUTINE> NAME="horz_diffuse"


!#######################################################################
! <SUBROUTINE NAME="delsq_tracer">
!
! <DESCRIPTION>
! Subroutine computes the laplacian operator acting on tracer with unit 
! diffusivity. Units of del2_tracer are tracer/length^2
! </DESCRIPTION>
!
subroutine delsq_tracer (Thickness, tracer, index)

  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: tracer
  integer, intent(in)                             :: index

  real, dimension(isd:ied,jsd:jed)                :: fe, fn
  real, dimension(isd:ied,jsd:jed,2)              :: tracr
  integer                                         :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_horz_diffuse_mod (delsq_tracer): module needs to be initialized')
  endif 


  ! fe = Diffusive flux across east face of T cells. Use unit diffusivity
  ! fn = Diffusive flux across north face of T cells. Use unit diffusivity

  tracr(:,:,1) = tracer(:,:,1)
  do k=1,nk
    tracr(:,:,2) = tracer(:,:,k)
    fe(:,:) = FDX_PT(tracr(:,:,1:2),k)*FMX(Thickness%dht(:,:,k,index)*Grd%tmask(:,:,k))
    fn(:,:) = FDY_PT(tracr(:,:,1:2),k)*FMY(Thickness%dht(:,:,k,index)*Grd%tmask(:,:,k))
    del2_tracer(:,:,k) = Grd%tmask(:,:,k)*(BDX_ET(fe(:,:)) + BDY_NT(fn(:,:)))/Thickness%dht(:,:,k,index)
    tracr(:,:,1) = tracr(:,:,2)
  enddo

end subroutine delsq_tracer
! </SUBROUTINE> NAME="delsq_tracer"


end module ocean_horz_diffuse_mod
      
      
