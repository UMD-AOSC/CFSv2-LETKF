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
#include <fms_platform.h>

module ocean_sponge_mod
!
!<CONTACT EMAIL="Bonnie.Samuels@noaa.gov"> Bonnie Samuels 
!</CONTACT>
!
!<CONTACT EMAIL="Robert.Hallberg@noaa.gov"> R.W. Hallberg 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="swathi@cmmacs.ernet.in"> P. S. Swathi  
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted tracer tendency [tracer*meter/sec] from sponges.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module applies sponges to tracers. The sponges
! can occur at any location and with any distribution in the domain, and
! with any time step and damping rate.  Sponges occur where positive
! inverse restore times occur in the field passed to sponge_init.  An
! array of tracer tendencies due to the sponges is augmented through a
! call to sponge_tracer_source.  The array of tracer tendencies must be
! reset to zero between calls.
!
! Different damping rates can be specified for each field by making
! calls to register_sponge_rate - no sponges are applied to fields for
! which uniformly zero inverse damping rates are set with a call to
! register_sponge_rate.  The value towards which a field is damped is
! set with calls to register_sponge_field; successive calls are used to
! set up linear interpolation of this restore rate.
!
! The user is responsible for providing (and registering) the data on
! the model grid of values towards which the tracers are being driven.
!</DESCRIPTION>
!
use fms_mod,                  only: write_version_number, open_namelist_file, close_file, file_exist, check_nml_error
use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
use mpp_mod,                  only: mpp_sum, mpp_chksum, mpp_error
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_time_type 
use ocean_workspace_mod,      only: wrk1

implicit none

private

#include <ocean_memory.h>

type ocean_sponge_type
   integer :: id ! time_interp_external index
   character(len=32) :: name ! tracer name corresponding to sponge
   real, dimension(:,:), _ALLOCATABLE :: damp_coeff _NULL ! the inverse damping rate (tracer units/ sec)
end type ocean_sponge_type


type(ocean_sponge_type), allocatable, dimension(:) :: Sponge
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_sponge_init
public sponge_tracer_source

character(len=126)  :: version = '$Id$'
character (len=128) :: tagname = '$Name$'

integer :: num_prog_tracers = 0

logical :: module_is_initialized = .FALSE.

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_sponge_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the sponges.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_sponge_init(Grid, Domain, T_prog, dtime)

  type(ocean_grid_type), target            :: Grid
  type(ocean_domain_type), target          :: Domain
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  real,intent(in)                          :: dtime

  logical :: error=.false.
  integer :: i, j, k, nt, tr_no, ioun, io_status, ierr, n
  real    :: dtimer

  character(len=128) :: name

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_sponge_mod (ocean_sponge_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  num_prog_tracers = size(T_prog(:))

  allocate( Sponge(num_prog_tracers) )

  call write_version_number( version, tagname )

! provide for namelist over-ride of defaults 
!  ioun =  open_namelist_file()
!  read  (ioun, ocean_sponge_nml,iostat=io_status)
!  write (stdlog(), ocean_sponge_nml)
!  ierr = check_nml_error(io_status, 'ocean_sponge_nml')
!  call close_file(ioun)

  Dom => Domain
  Grd => Grid

#ifndef STATIC_MEMORY    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  do n = 1, num_prog_tracers
     Sponge(n)%id = -1
  enddo

  dtimer = 1.0/dtime

  do n = 1, num_prog_tracers

     ! read damping rates (1/sec) 

     name = 'INPUT/'//trim(T_prog(n)%name)//'_sponge_coeff.nc'
     if (file_exist(trim(name))) then
        write(stdout(),*) '==> Using sponge restoring values specified from file '//trim(name) 
        allocate(Sponge(n)%damp_coeff(isd:ied,jsd:jed))
        call read_data(name,'coeff',Sponge(n)%damp_coeff,domain=Domain%domain2d,timelevel=1)

        where (Grd%tmask(:,:,1) == 0.0)
           Sponge(n)%damp_coeff = 0.0
        end where

        ! modify damping rates to allow restoring to be solved implicitly
        ! note: test values between zero and 4.0e-8 revert to damping rates defined above

        do j=jsc,jec
           do i=isc,iec
!              do k=1,nk
                 if (dtime*Sponge(n)%damp_coeff(i,j) > 4.0e-8) then
                    if (dtime*Sponge(n)%damp_coeff(i,j) > 37.0) then
                       Sponge(n)%damp_coeff(i,j) = dtimer
                    else
                       Sponge(n)%damp_coeff(i,j) = (1.0 - exp(-dtime*Sponge(n)%damp_coeff(i,j))) * dtimer
                    endif
                 else if (dtime*Sponge(n)%damp_coeff(i,j) <= 0.0) then
                    Sponge(n)%damp_coeff(i,j) = 0.0
                 endif
!              enddo
           enddo
        enddo

        ! read restoring data

        name = 'INPUT/'//trim(T_prog(n)%name)//'_sponge.nc'
        Sponge(n)%id = init_external_field(name,T_prog(n)%name,domain=Domain%domain2d)
        if (Sponge(n)%id < 1) then 
          call mpp_error(FATAL,'==>Error: in ocean_spong_mod: damping rates are specified but sponge values are not')
        endif 
        write(stdout(),*) '==> Using sponge data specified from file '//trim(name) 
     else
        write(stdout(),*) '==> '//trim(name)//' not found.  Sponge not being applied '
     endif
  enddo


end subroutine ocean_sponge_init
! </SUBROUTINE> NAME="ocean_sponge_init"


!#######################################################################
! <SUBROUTINE NAME="sponge_tracer_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted time tendencies of the
! tracers due to the sponges.
! </DESCRIPTION>
!
subroutine sponge_tracer_source(Time, Thickness, T_prog)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: n, taum1, tau, k

  taum1 = Time%taum1
  tau   = Time%tau
  
  do n = 1, size(T_prog(:))
     if (Sponge(n)%id > 0) then
        call time_interp_external(Sponge(n)%id, Time%model_time, wrk1)  ! get sponge value for current model time
        do k=1, nk
           T_prog(n)%th_tendency(isc:iec,jsc:jec,k) = T_prog(n)%th_tendency(isc:iec,jsc:jec,k) + &
             Thickness%dht(isc:iec,jsc:jec,k,tau)*Sponge(n)%damp_coeff(isc:iec,jsc:jec)*(wrk1(isc:iec,jsc:jec,k) - &
             T_prog(n)%field(isc:iec,jsc:jec,k,taum1))
        enddo
     endif
  enddo

  return

end subroutine sponge_tracer_source
! </SUBROUTINE> NAME="sponge_tracer_source"


end module ocean_sponge_mod
