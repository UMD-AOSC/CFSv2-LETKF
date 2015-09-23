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
module ocean_rivermix_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Keith.Dixon@noaa.gov">  K.W. Dixon 
!</CONTACT>
!
!<OVERVIEW>
! Tracer source from discharging river with depth or 
! mixing rivers with depth. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted tendency [tracer*meter/sec]
! associated with discharge of river tracer content 
! over a user defined column of ocean points. Points are
! selected based on whether river flow into a point is nonzero.
! Contribution added to tracer source array.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Guide to MOM4 (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <NOTE>
! Algorithm ensures total tracer is conserved.  Note that volume/mass is 
! modified by river water within the eta-equation using the big leap-frog.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_rivermix_nml">
!  <DATA NAME="river_insertion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to insert tracers from 
!  rivers. 
!  </DATA>
!  <DATA NAME="river_diffusion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to diffuse tracers from 
!  rivers. 
!  </DATA> 
!  <DATA NAME="river_diffusivity" TYPE="real" UNITS="m^2/s">
!  Vertical diffusivity enhancement at river mouths which is applied 
!  to a depth of river_diffusion_thickness, with linear tapering to zero
!  enhancement from the ocean surface to river_diffusion_thickness. 
!  </DATA> 
!  <DATA NAME="river_diffuse_temp" TYPE="logical">
!  Logical to determine if enhance vertical diffusion of temp at river mouths 
!  </DATA> 
!  <DATA NAME="river_diffuse_salt" TYPE="logical">
!  Logical to determine if enhance vertical diffusion of salt and all other 
!  passive tracers at river mouths 
!  </DATA> 
!  <DATA NAME="debug_river" TYPE="logical">
!  For debugging 
!  </DATA> 
!</NAMELIST>
!

use axis_utils_mod,      only: frac_index, nearest_index
use constants_mod,       only: epsln, rho0, rho0r
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog 
use mpp_domains_mod,     only: domain2D
use mpp_mod,             only: mpp_error

use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_thickness_type 
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_external_mode_type
use ocean_workspace_mod, only: wrk1 

implicit none

private 

public ocean_rivermix_init
public rivermix
private river_discharge_tracer
private river_kappa

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

real    :: river_insertion_thickness=0.0 ! static thickness (m) of ocean column where discharge rivers.  
                                         ! actual thickness is based on model grid spacing. min thickness=dtz(1).
real    :: river_diffusion_thickness=0.0 ! static thickness (m) of ocean column where diffuse tracer at river mouths.    
                                         ! actual thickness is based on model grid spacing. min thickness=dtz(1).      
real    :: river_diffusivity=0.0         ! enhancement to the vertical diffusivity (m^2/s) at river mouths 
logical :: river_diffuse_temp=.false.    ! to enhance diffusivity of temp at river mouths over river_thickness column
logical :: river_diffuse_salt=.false.    ! to enhance diffusivity of salt at river mouths over river_thickness column

logical :: river_discharge=.true.        ! internally set for discharge of river tracers into a vertical column
logical :: debug_river=.false.           ! for debugging

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_rivermix
integer :: id_diff_cbt_river_t=-1
integer :: id_diff_cbt_river_s=-1
integer :: unit=6       ! processor zero writes to unit 6

! time step (sec)
real    :: dtime 

integer :: num_prog_tracers=0

character(len=128) :: version=&
       '$Id$'
character (len=128) :: tagname=&
     '$Name$'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

logical :: module_is_initialized = .FALSE.

namelist /ocean_rivermix_nml/ river_diffuse_temp, river_diffuse_salt, &
                              river_insertion_thickness, river_diffusion_thickness, &
                              river_diffusivity, debug_river

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_rivermix_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers at river mouths. 
! </DESCRIPTION>
!
  subroutine ocean_rivermix_init(Grid, Domain, Time, Time_steps, T_prog, debug)

    type(ocean_grid_type), target               :: Grid
    type(ocean_domain_type), target             :: Domain
    type(ocean_time_type), target               :: Time 
    type(ocean_time_steps_type), intent(in)     :: Time_steps
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    logical, intent(in), optional               :: debug

    integer :: n, nz_insert, nz_diffuse, model
    integer :: io_status, ioun, ierr
    real    :: ztop=0.0

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error from ocean_rivermix_mod (ocean_rivermix_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk

    Dom => Domain
    Grd => Grid

    dtime = Time_steps%dtts

    num_prog_tracers = size(T_prog(:))

    if (PRESENT(debug)) debug_river = debug

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_rivermix_nml,IOSTAT=io_status)
    write (stdlog(),ocean_rivermix_nml)
    write (stdout(),'(/)') 
    write (stdout(),ocean_rivermix_nml)
    ierr = check_nml_error(io_status,'ocean_rivermix_nml')
    call close_file (ioun)

    ! get nominal thickness over which to distribute river runoff 
    nz_insert  = nearest_index(river_insertion_thickness,Grd%zw)
    nz_diffuse = nearest_index(river_diffusion_thickness,Grd%zw)    
    nz_insert  = max(1,nz_insert)
    nz_diffuse = max(1,nz_diffuse)    
   

    if(river_insertion_thickness==0.0) then 
       call mpp_error(NOTE, '==>Note: resetting river_insertion_thickness to dzt(1).  Discharge all river contents into top cell')  
    else 
      write(stdout(),'(a,i3,a)')'==>Note: if using waterflux, discharge river tracer over ',nz_insert,' grid points in vertical'
    endif

    if(river_diffuse_temp) then 
        if(Time_steps%aidif==0.0) then
            call mpp_error(FATAL, ' ==>Error: aidif=1.0 must be set to allow stable time stepping with river_diffuse_temp=.true.')  
        endif
        if(river_diffusion_thickness==0.0) then 
            call mpp_error(NOTE, ' ==>Note: resetting river_diffusion_thickness=dzt(1) for enhanced diff_cbt of temperature.')  
        else  
            write(stdout(),'(a,i3,a)')' ==>Note: enhance temp diff_cbt at rivers over ',nz_diffuse,' points in vertical'
            write(stdout(),'(a,e12.5,a)')'          Do so by adding diffusivity of ',river_diffusivity,' m^2/s to diff_cbt(temp)'
        endif
    endif

    if(river_diffuse_salt) then 
        if(Time_steps%aidif==0.0) then
            call mpp_error(FATAL, ' ==>Error: aidif=1.0 must be set to allow stable time stepping with river_diffuse_salt=.true.')  
        endif
        if(river_diffusion_thickness==0.0) then 
            call mpp_error(NOTE, ' ==>Note: resetting river_diffusion_thickness=dzt(1) for enhanced diff_cbt of salt and passive')  
        else 
            write(stdout(),'(a,i3,a)')' ==>Note: enhance salt/passive diff_cbt at rivers over ',nz_diffuse,' points in vertical'
            write(stdout(),'(a,e12.5,a)')'          Do so by adding diffusivity of ',river_diffusivity,' m^2/s to diff_cbt(salt)'
        endif
    endif

    ! register for diag_manager 
    allocate (id_rivermix(num_prog_tracers))
    id_rivermix = -1

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_rivermix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rivermix', &
                Grd%tracer_axes(1:3), Time%model_time, &
                'rho_cp*rivermix*dht*temp', &
                'Watt/m^2', missing_value=-1.e10, range=(/-1.e10,1.e10/))
       else
           id_rivermix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rivermix', &
                Grd%tracer_axes(1:3), Time%model_time, &
                'rho0*rivermix*dht*tracer for'//trim(T_prog(n)%name),&
                trim(T_prog(n)%units)//' m*kg/sec', missing_value=-1.e10, range=(/-1.e10,1.e10/))
       endif
    enddo

   id_diff_cbt_river_t = -1
   id_diff_cbt_river_t = register_diag_field ('ocean_model', 'diff_cbt_river_t', &
                            Grid%tracer_axes(1:3), Time%model_time, &
                            'diff_cbt(temp) enhancement at rivers', 'm^2/s',missing_value=-10.0, &
                            range=(/-10.0,1e6/))

   id_diff_cbt_river_s = -1
   id_diff_cbt_river_s = register_diag_field ('ocean_model', 'diff_cbt_river_s', &
                         Grid%tracer_axes(1:3), Time%model_time, &
                         'diff_cbt(salt) enhancement at rivers', 'm^2/s',missing_value=-10.0, &
                         range=(/-10.0,1e6/))

  end subroutine ocean_rivermix_init
! </SUBROUTINE> NAME="ocean_rivermix_init"




!#######################################################################
! <SUBROUTINE NAME="rivermix">
!
! <DESCRIPTION>
! This subroutine computes one or all of the following: 
!
! (1) Thickness weighted tracer source associated with 
! river tracer content discharged into a vertical column of ocean 
! tracer cells. This is done if river_discharge=.true.
!
! (2) Enhance vertical diffusivity at river mouths. 
! This is done if river_diffuse_temp=.true. or 
! river_diffuse_salt=.true. 
!
! Doing one or both are useful for models with fine vertical  
! resolution, where discharging river content to top cell 
! is often not numerically suitable nor physically relevant.
!
! </DESCRIPTION>
!
subroutine rivermix (Time, Thickness, T_prog, Ext_mode, river, diff_cbt, index_temp, index_salt)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_thickness_type), intent(in)               :: Thickness
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_external_mode_type), intent(in)           :: Ext_mode
  integer, intent(in)                                  :: index_temp
  integer, intent(in)                                  :: index_salt
  real, intent(in), dimension(isd:ied,jsd:jed)         :: river
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt
  integer :: i,j

  if(.not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_rivermix_mod (rivermix): module must be initialized')
  endif 

  river_discharge=.false. 
  do j=jsc,jec
     do i=isc,iec
       if(river(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) river_discharge=.true.
     enddo
  enddo
  if(river_discharge) call river_discharge_tracer(Time, Thickness, T_prog(1:num_prog_tracers), Ext_mode, river)
  if(river_diffuse_temp) call river_kappa(Time, T_prog(index_temp), diff_cbt(isd:ied,jsd:jed,:,1))
  if(river_diffuse_salt) call river_kappa(Time, T_prog(index_salt), diff_cbt(isd:ied,jsd:jed,:,2) )

end subroutine rivermix
! </SUBROUTINE> NAME="rivermix"


!#######################################################################
! <SUBROUTINE NAME="river_discharge_tracer">
!
! <DESCRIPTION>
! Compute thickness weighted tracer source [tracer*m/s]
! associated with the discharge of tracer from a river over 
! a vertical column whose thickness is set by River_insertion_thickness 
! and whose horizontal location is given by the river array. 
! </DESCRIPTION>
!
subroutine river_discharge_tracer (Time, Thickness, T_prog, Ext_mode, river)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_thickness_type), intent(in)        :: Thickness
  type(ocean_prog_tracer_type), intent(inout)   :: T_prog(:)
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed)  :: river

  integer :: i, j, k, n, nz
  integer :: tau
  real    :: depth, thkocean
  real    :: delta(nk), delta_rho_tocean(nk), delta_rho0_triver(nk)
  real    :: zextra, zinsert, tracerextra, tracernew(nk)
  
  
  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_rivermix_mod (river_discharge_tracer): module must be initialized')
  endif 

  tau = Time%tau
  
  do n = 1, num_prog_tracers  
    T_prog(n)%wrk1(:,:,:) = 0.0
  enddo  
  delta             = 0.0
  delta_rho_tocean  = 0.0
  delta_rho0_triver = 0.0
  wrk1              = 0.0
               
  do j=jsc,jec
     do i=isc,iec

        if (river(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then


!! river contains the volume rate (m/s) of fluid with tracer 
!! concentration triver to be distributed to underlying levels
            

            depth    = min(Grd%ht(i,j),river_insertion_thickness)        ! be sure not to discharge river content into rock 
            nz       = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! number of k-levels into which discharge rivers
            nz       = max(1,nz)                                         ! make sure have at least one cell to discharge into

            thkocean = 0.0
            do k=1,nz
               thkocean = thkocean + Thickness%dht(i,j,k,tau)
            enddo
            do k=1,nz
               delta(k) = Thickness%dht(i,j,k,tau)/(epsln+thkocean)
            enddo


            do n=1,num_prog_tracers

! mix a portion of the river volume inserted over the course of a tracer timestep
! into each model level from k=1:nz in proportion to the fractional cell thickness.
! The mixing scheme used here follows scheme devised by Keith Dixon to handle hydraulic
! control across xland mixing points
!
!
! The following schematic illustrates the scheme assuming constant layer thickness
! neglecting free surface height variations and assuming the river insertion depth
! extends to the base of the third model level               
!               
!                  Total depth = H; cell thickness = dz; River volume flux/timestep = R
!  x--------x
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of Triver + 2*zinsert meters
!  !        !        of water at T*(2)
!  !        !        T*(1) = (T(1).dz + 3.T*(2).zinsert + Triver.zinsert ) / (dz + 3.zinsert)               
!  x--------x               
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of Triver + zinsert meters
!  !        !        of water at T*(3)           
!  !        !        T*(2) = (T(2).dz + 2.T*(3).zinsert + Triver.zinsert ) / (dz + 2.zinsert)               
!  !        !
!  x--------x                              
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of Triver
!  !        !        T*(3) = (T(3).dz + Triver.zinsert ) / (dz + zinsert)
!  !        !
!  xxxxxxxxxx               
               
               zextra=0.0
               do k=nz,1,-1
                  tracernew(k) = 0.0

                  if (k.eq.nz) then
                      tracerextra=0.0
                  else
                      tracerextra = tracernew(k+1)
                  endif

                  zinsert = river(i,j)*dtime*delta(k)
                  tracernew(k) = (tracerextra*zextra + T_prog(n)%field(i,j,k,tau)*Thickness%dht(i,j,k,tau) + &
                                  T_prog(n)%triver(i,j)*zinsert) / (zextra+Thickness%dht(i,j,k,tau)+zinsert)
                  
                  zextra=zextra+zinsert
               enddo


               k=1
               T_prog(n)%wrk1(i,j,k) = rho0*(tracernew(k)*(Thickness%dht(i,j,k,tau)+river(i,j)*dtime) -&
                             T_prog(n)%field(i,j,k,tau)*Thickness%dht(i,j,k,tau))/dtime
               
               do k=2,nz
                  T_prog(n)%wrk1(i,j,k) = rho0*Thickness%dht(i,j,k,tau)*(tracernew(k) - T_prog(n)%field(i,j,k,tau))/dtime
               enddo

               do k=1,nz
                  T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)*rho0r
               enddo

               if (debug_river) then
                   write(6,*) 'i,j,n,river(m/sec),triver= ',i,j,n,river(i,j),T_prog(n)%triver(i,j)
                   do k=1,nz
                      write(6,*) 'i,j,k,tracer(k),tracernew(k),tendency(tracer*thickness/sec)= ',&
                           i,j,k,T_prog(n)%field(i,j,k,tau),tracernew(k),T_prog(n)%wrk1(i,j,k)*rho0r
                   enddo
               endif
                
               
               
            enddo ! n end  

        endif ! river > 0

     enddo   ! i end
  enddo      ! j end


  ! fill source array and send diagnostics

  do n=1,num_prog_tracers 
    if(id_rivermix(n) > 0) then 
      used = send_data (id_rivermix(n), T_prog(n)%wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 
  enddo 


end subroutine river_discharge_tracer
! </SUBROUTINE> NAME="river_discharge_tracer"



!#######################################################################
! <SUBROUTINE NAME="river_kappa">
!
! <DESCRIPTION>
! This subroutine enhances the vertical diffusivity kappa over 
! a vertical column whose thickness is set by river_diffusion_thickness 
! and whose horizontal location is given by the rmask array.
! Note that rmask can be > 0 even if river=0 in the case when 
! use virtual salt flux.   
! The enhanced diffusivity is maximum at the top cell and is linearly 
! interpolated to the normal diffusivity at the depth set by 
! river_diffusion_thickness
! </DESCRIPTION>
!
subroutine river_kappa (Time, Tracer, kappa)

  type(ocean_time_type), intent(in)                  :: Time
  type(ocean_prog_tracer_type), intent(in)           :: Tracer
  real, intent(inout), dimension(isd:ied,jsd:jed,nk) :: kappa

  integer :: i, j, k, nz
  real    :: depth

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_rivermix_mod (river_kappa): module must be initialized')
  endif 

  wrk1=0.0

  do j=jsc,jec
     do i=isc,iec
        if (Tracer%riverdiffuse(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then
            depth = min(Grd%ht(i,j),river_diffusion_thickness)
            nz    = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw)))
            nz    = max(1,nz)                                
            do k=1,nz 
              wrk1(i,j,k)  = river_diffusivity*(1.0 - Grd%zt(k)/Grd%zt(nz))
              kappa(i,j,k) = kappa(i,j,k) + wrk1(i,j,k)
            enddo
        endif
     enddo
  enddo

  if (id_diff_cbt_river_t > 0 .and. Tracer%name=='temp') then 
      used = send_data (id_diff_cbt_river_t, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_diff_cbt_river_s > 0 .and. Tracer%name=='salt') then 
      used = send_data (id_diff_cbt_river_s, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

end subroutine river_kappa
! </SUBROUTINE> NAME="river_kappa"


end module ocean_rivermix_mod
      
      




