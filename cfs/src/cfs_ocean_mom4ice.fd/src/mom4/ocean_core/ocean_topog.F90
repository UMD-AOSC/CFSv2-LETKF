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
module ocean_topog_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Set up ocean bottom topography. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Set up ocean bottom topography. Reads information from grid specification file. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_topog_nml">
!  <DATA NAME="flat_bottom" TYPE="logical">
!  For debugging, it is often useful to over-ride the grid spec file
!  and simply make the domain flat bottom. 
!  </DATA> 
!  <DATA NAME="flat_bottom_kmt" TYPE="integer">
!  Number of depth levels to use for the flat_bottom option.
!  </DATA> 
!  <DATA NAME="flat_bottom_ht" TYPE="real">
!  Depth to make the flat_bottom. 
!  </DATA> 
! </NAMELIST>
!
use fms_mod,         only: open_namelist_file, close_file, check_nml_error
use mpp_domains_mod, only: mpp_update_domains
use mpp_io_mod,      only: mpp_open, mpp_read, mpp_get_info, mpp_get_atts, mpp_get_fields
use mpp_io_mod,      only: MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
use mpp_io_mod,      only: fieldtype
use mpp_mod,         only: mpp_error, mpp_min, mpp_pe, FATAL, NOTE, stdout, stdlog

use ocean_domains_mod, only: get_local_indices, get_domain_offsets
use ocean_types_mod,   only: ocean_grid_type, ocean_domain_type

implicit none
private

character(len=256) :: version='CVS $Id$'
character(len=256) :: tagname='Tag $Name$'

#include <ocean_memory.h>

! for output
integer :: writeunit=6

logical :: flat_bottom=.false.
integer :: flat_bottom_kmt=50
real    :: flat_bottom_ht=5500.0

namelist /ocean_topog_nml/ flat_bottom, flat_bottom_kmt, flat_bottom_ht

public ocean_topog_init

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_topog_init">
!
! <DESCRIPTION>
! Initialize the ocean bottom topography.
! </DESCRIPTION>
!
subroutine ocean_topog_init (Domain, Grid, grid_file)

  type(ocean_domain_type), intent(inout) :: Domain
  type(ocean_grid_type), intent(inout)   :: Grid
  character(len=*), intent(in), optional :: grid_file

  real, dimension(:,:), allocatable          :: tmp
  type(fieldtype), allocatable, dimension(:) :: fields
  character(len=128)                         :: name, grd_file
  logical                                    :: land_jeq1
  logical                                    :: found_ht, found_kmt
  logical                                    :: is_new_grid
  real                                       :: min_depth, min_depth0, fudge 
  integer                                    :: imin, jmin, kmin
  integer                                    :: unit, ioun, io_status, ierr
  integer                                    :: n, nvar, ndim, natt, ntime
  integer                                    :: i, j, ioff, joff

  grd_file = 'INPUT/grid_spec'

  write( stdlog(),'(/a/)') trim(version)

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (Grid%kmt(isd:ied,jsd:jed))
  allocate (Grid%ht(isd:ied,jsd:jed))
  allocate (Grid%kmu(isd:ied,jsd:jed))
  allocate (Grid%hu(isd:ied,jsd:jed))
#else
  call get_domain_offsets(Domain,ioff,joff)
#endif
  

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_topog_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_topog_nml)  
  write (stdlog(),ocean_topog_nml)
  ierr = check_nml_error(io_status, 'ocean_topog_nml')
  call close_file(ioun)
  
  if (PRESENT(grid_file)) grd_file = trim(grid_file)
  call mpp_open(unit,trim(grd_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
  allocate(fields(nvar))
  call mpp_get_fields(unit, fields)

  Grid%kmt=0;Grid%kmu=0;Grid%ht=0;Grid%hu=0
  found_ht = .false. ;  found_kmt = .false.

  !--- determine if new grid or old grid
  do i=1, nvar
     call mpp_get_atts(fields(i),name=name)
     select case (name)
        case ('depth_t')
           is_new_grid = .true.
        case ('ht')
           is_new_grid = .false.
     end select
  enddo

  do i=1, nvar
     call mpp_get_atts(fields(i),name=name)
     if(is_new_grid) then
        select case (name)
           case ('depth_t')
              call mpp_read(unit,fields(i),Domain%domain2d,Grid%ht)
              found_ht = .true.
           case ('num_levels') 
              allocate(tmp(isd:ied,jsd:jed));tmp=0
              call mpp_read(unit,fields(i),Domain%domain2d,tmp)
              found_kmt = .true.
        end select
     else
        select case (name)
           case ('ht')
              call mpp_read(unit,fields(i),Domain%domain2d,Grid%ht)
              found_ht = .true.
           case ('kmt') 
              allocate(tmp(isd:ied,jsd:jed));tmp=0
              call mpp_read(unit,fields(i),Domain%domain2d,tmp)
              found_kmt = .true.
        end select
     endif
  enddo
  if(is_new_grid) then
     if(.not. found_ht)  call mpp_error(FATAL,'ocean_topog_mod: field depth_t is not in file '//trim(grd_file) )  
     if(.not. found_kmt) call mpp_error(FATAL,'ocean_topog_mod: field num_levels is not in file '//trim(grd_file) )
  else
     if(.not. found_ht)  call mpp_error(FATAL,'ocean_topog_mod: field ht is not in file '//trim(grd_file) )  
     if(.not. found_kmt) call mpp_error(FATAL,'ocean_topog_mod: field kmt is not in file '//trim(grd_file) )
  endif

  call mpp_update_domains(Grid%ht,Domain%domain2d)
  call mpp_update_domains(tmp,Domain%domain2d)
  Grid%kmt=tmp

  if(flat_bottom) then 
      call mpp_error(NOTE,&
      '==>ocean_topog_init: flat_bottom=.true. so will set all cells that are not land equal to their max depth.')
      do j=jsc,jec
         do i=isc,iec
            if(Grid%kmt(i,j) > 0) then  
               Grid%kmt(i,j) = flat_bottom_kmt
               Grid%ht(i,j)  = flat_bottom_ht
            endif  
         enddo
      enddo
  endif


! construct depths at u cell points as min of surrounding t cells

  do j=jsc,jec
     do i=isc,iec
        Grid%kmu(i,j) = min(Grid%kmt(i,j), Grid%kmt(i+1,j), Grid%kmt(i,j+1), Grid%kmt(i+1,j+1))
        Grid%hu(i,j)  = min(Grid%ht(i,j), Grid%ht(i+1,j), Grid%ht(i,j+1), Grid%ht(i+1,j+1))
     enddo
  enddo
  call mpp_update_domains(Grid%kmu,Domain%domain2d)
  call mpp_update_domains(Grid%hu,Domain%domain2d)

  deallocate (tmp)


! check whether any ocean occupies j=1 
! The FMS Sea Ice Simulator (SIS) requires land to be present
! for all points at jsc+joff=1.  If this is not the case, then 
! the ocean model cannot be coupled to SIS.  The SIS requirement 
! of all land at jsc+joff=1 is related to the use of a northeast 
! B-grid convention.  To couple the models, the ocean grid
! must be generated with fill_first_row=.true. 
  if (jsc+joff==1) then
    land_jeq1 = .false. 
    do i=isc,iec
      if(Grid%kmt(i,1) > 0) land_jeq1 = .true. 
    enddo
    if(land_jeq1) then 
      call mpp_error(NOTE,&
      '==>ocean_topog_init: FMS/Sea Ice (SIS) needs land at j=1. Presently have ocean at j=1, so cannot couple to SIS.')
    endif 
  endif 

! find shallowest water and provide caveat for overly shallow regions
  imin= 0 ; jmin=0 ; kmin=0 ; min_depth = 1.e10
  do j=jsc,jec
    do i=isc,iec
      if(Grid%kmt(i,j) > 0 .and. Grid%ht(i,j) < min_depth) then 
        imin      = i
        jmin      = j
        kmin      = Grid%kmt(i,j)
        min_depth = Grid%ht(i,j)
      endif 
    enddo
  enddo 
  fudge      = 1 + 1.e-12*mpp_pe() ! to distinguish processors when min_depth is independent of processor
  min_depth  = min_depth*fudge 
  min_depth0 = min_depth
  call mpp_min(min_depth)  

  if(min_depth0 == min_depth) then 
      write(writeunit,'(/a,f12.5)')' The minimum ocean depth is (meters) ', min_depth  
      write(writeunit,'(a,i3,a,i3,a,i3,a)')' and this occurs at (i,j,k) = (',&
                      imin+Domain%ioff,',',jmin+Domain%joff ,',',kmin,')'  
      write(writeunit,'(a,f10.4,a,f10.4,a,f12.5,a/)')' which has (long,lat,depth) =  (' &
           ,Grid%xt(imin,jmin),',',Grid%yt(imin,jmin), ',', Grid%ht(imin,jmin),')'
      if(min_depth < 50.0) then 
          write(writeunit,'(a)')' Beware that overly shallow regions (e.g., those shallower than 50m) may be subject' 
          write(writeunit,'(a)')' to numerical problems if strong surface forcing is not adequately mixed vertically.' 
          write(writeunit,'(a)')' Such problems may occur especially in shallow regions with kmt==2.  A symptom is the'
          write(writeunit,'(a)')' the current speed and/or tracer deviations becoming overly large due to the deposition'
          write(writeunit,'(a)')' of wind and/or buoyancy over just a small upper ocean region. Such problems can be'
          write(writeunit,'(a)')' resolved by adding sufficient vertical mixing in these regions, as indeed happens in'
          write(writeunit,'(a)')' Nature due to such phenomena as tides and breaking surface waves.'
      endif
  endif


end subroutine ocean_topog_init
! </SUBROUTINE> NAME="ocean_topog_init"


end module ocean_topog_mod

