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
module idealized_ic_mod 
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  !
  !<CONTACT EMAIL="Zhi.Liang@noaa.gov"> Z. Liang  </CONTACT>
  !<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies </REVIEWER>

  !<OVERVIEW>
  ! For preparing idealized mom4 initial conditions
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This program prepares idealized initial conditions (active and passive tracer). 
  ! Results are netcdf files that are then read into mom4.
  ! Various idealized options are available as selected by a namelist.
  !</DESCRIPTION>
  !

  use fms_mod,         only : write_version_number, stdout, open_namelist_file
  use fms_mod,         only : file_exist, check_nml_error, read_data
  use mpp_mod,         only : mpp_chksum, mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, NOTE 
  use mpp_io_mod,      only : mpp_open, mpp_close, mpp_write_meta, mpp_write, mpp_get_axis_data
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_get_atts, mpp_get_axis_data, axistype
  use mpp_io_mod,      only : MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_OVERWR, MPP_NETCDF, fieldtype
  use mpp_domains_mod, only : domain2d, mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_global_field, mpp_domains_set_stack_size
  use constants_mod,   only : RADIAN

  implicit none
  private

  integer, parameter :: max_tracer = 10

  !--- namelist ----------------------------------------------------------
  !<NAMELIST NAME="idealized_ic_nml">
  !  <DATA NAME="temp_type" TYPE="character(len=128)" >
  !  Control the iealized initial temperature condition options. 
  !  There are five options available and the default value is 
  !  "constant_temp_ic". When temp_type is
  !  1. = "constant_temp_ic", use spatially constant initial potential temperature.
  !  2. = "exponential_temp_ic", use initial potential temperature that is 
  !        exponential in the vertical.
  !  3. = "equatorial_temp_ic", use initial temp condition for idealized equatorial studies.
  !  4. = "shelfbowl_temp_ic", Initial conditions for Winton etal shelf-bowl test case.  
  !        Use tanh transition between cold shelf and warm bowl waters, instead of Heaviside. 
  !  5. = "zonal_levitus_temp_ic", use zonal average Levitus temp as initial conditions. 
  !  </DATA>
  !  <DATA NAME="salt_type" TYPE="character(len=128)" >
  !  Control the idealized initial salinity condition options. 
  !  There are four options available and the default value is 
  !  "constant_salt_ic". When salt_type is 
  !  1. = "constant_salt_ic", use spatially constant initial salinity.
  !  2. = "exponential_salt_ic", use initial salinity that is exponential in the vertical.
  !  3. = "salinity_profile_ic", use initial salinity condition as set by profile in function salt0
  !  </DATA>
  !  <DATA NAME="age_type" TYPE="character(len=128)" >
  !  Control the idealized initial age condition options. 
  !  There is only one option now available, and the default value is 
  !  "constant_age_ic". 
  !  </DATA>
  !  <DATA NAME="constant_temp_value" TYPE="real">
  !  Value for uniform initial temp.
  !  </DATA>
  !  <DATA NAME="constant_salt_value" TYPE="real">
  !  Value for uniform initial salinity.
  !  </DATA>
  !  <DATA NAME="constant_age_value" TYPE="real">
  !  Value for uniform initial age.
  !  </DATA>
  !  <DATA NAME="num_active_tracer" TYPE="integer" >
  !  Number of active tracers will be generated. Its value should be 0, 1, 2.
  !  Its default value is 2. ( temp and salt )
  !  </DATA>
  !  <DATA NAME="active_tracer" TYPE="character(len=128),dimension(2)" >
  !  name array of the active tracers will be generated. 
  !  Its element value should be 'temp' or 'salt'
  !  </DATA>
  !  <DATA NAME="active_tracer_file" TYPE="character(len=128)" >
  !  active tracer output file.
  !  </DATA>
  !  <DATA NAME="num_passive_tracer" TYPE="integer" >
  !  Number of passive tracers will be generated. Its value should be no more than max_tracer.
  !  Its default value is 1 (age). 
  !  </DATA>
  !  <DATA NAME="" TYPE="character(len=128),dimension(max_tracer)" >
  !  name array of the passive tracers will be generated. 
  !  Its element value should be 'age'.
  !  </DATA>
  !  <DATA NAME="passive_tracer_file" TYPE="character(len=128)" >
  !  passive tracer output file.
  !  </DATA>
  !  <DATA NAME="grid_file" TYPE="character(len=128)" >
  !  grid descriptor file.
  !  </DATA>
  !  <DATA NAME="t1" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="t0" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="z0" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="thk" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  ! </NAMELIST>

  character(len=128) :: temp_type = 'constant_temp_ic'
  character(len=128) :: salt_type = 'constant_salt_ic'
  character(len=128) :: age_type  = 'constant_age_ic'

  real               :: constant_temp_value    = 12.0
  real               :: constant_salt_value    = 34.72
  real               :: constant_age_value     = 0.0
  real               :: t0=7.5, t1=10.0, z0=30.0, thk=80.0
  integer            :: num_active_tracer = 2   
  character(len=128) :: grid_file = 'grid_spec.nc'
  character(len=128) :: active_tracer(2)
  character(len=128) :: active_tracer_file = 'ocean_temp_salt.res.nc'
  integer            :: num_passive_tracer = 1
  character(len=128) :: passive_tracer(max_tracer)
  character(len=128) :: passive_tracer_file = 'ocean_age.res.nc'

  namelist /idealized_ic_nml/ temp_type, salt_type, constant_temp_value, constant_salt_value, &
       t1, t0, z0, thk,  num_active_tracer, active_tracer, active_tracer_file,           &
       num_passive_tracer, passive_tracer, passive_tracer_file, grid_file

  !--- version information 
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$' 

  !--- output data -------------------------------------------------------
  real, dimension(:,:,:), allocatable :: temp, salt, age
  real, dimension(:,:,:), allocatable :: global_temp, global_salt, global_age
  !--- other variables
  logical                           :: generate_temp_ic = .FALSE.
  logical                           :: generate_salt_ic = .FALSE.
  logical                           :: generate_age_ic  = .FALSE.
  real, dimension(:,:), allocatable :: xt, yt, kmt 
  real, dimension(:),   allocatable :: zt, grid_x_t, grid_y_t
  logical                           :: module_is_initialized = .FALSE.
  integer                           :: ni, nj, nk, isc, iec, jsc, jec
  type(domain2d),save               :: domain
  !--- public interface --------------------------------------------------
  public  idealized_ic_init, write_idealized_ic_data, idealized_ic_end

contains

  !#######################################################################
  ! <SUBROUTINE NAME="idealized_ic_init">
  !
  ! <DESCRIPTION>
  ! Initialize the module generating ideal initial conditions.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine idealized_ic_init

    !--- local variables -------------------------------------------------
    integer :: unit, io_status, ierr, n
    integer :: chksum_t, chksum_s, chksum_a
    !---------------------------------------------------------------------

    if ( module_is_initialized ) call mpp_error(FATAL, 'idealized_ic_mod: attempted reinitialization')

    module_is_initialized = .TRUE.

    !--- default value of tracer_name
    active_tracer(1)  = 'temp'
    active_tracer(2)  = 'salt'
    passive_tracer(1) = 'age'

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, idealized_ic_nml,iostat=io_status, end = 10)
          ierr = check_nml_error(io_status,'idealized_ic_nml')
       enddo
10     call mpp_close (unit)
    else
       call mpp_error(NOTE,'idealized_ic_mod: file input.nml does not exist, use default namelist option')  
    endif

    if(num_active_tracer .gt. 2) call mpp_error(FATAL, &
            'idealized_ic_mod: num_active_tracer should be the value of 0, 1 or 2')

    if(num_passive_tracer .gt. max_tracer) call mpp_error(FATAL, &
            'idealized_ic_mod: num_passive_tracer should be no greater than max_tracer')   
 
    do n = 1, num_active_tracer
       select case(active_tracer(n))
       case('temp')
          generate_temp_ic = .true.
       case('salt')
          generate_salt_ic = .true.
       case default
          call mpp_error(FATAL, trim(active_tracer(n))//' is not a valid option of nml active_tracer')
       end select
    enddo

    do n = 1, num_passive_tracer
       select case(passive_tracer(n))
       case('age')
          generate_age_ic = .true.
       case default
          call mpp_error(FATAL, trim(passive_tracer(n))//' is not a valid option of nml passive_tracer')
       end select
    enddo
    !--- at least one of generate_temp_ic and generate_salt_ic should be true
    if(.not.generate_temp_ic .and. .not.generate_salt_ic .and. .not.generate_age_ic)  &
      call mpp_error(FATAL, &
         'idealized_ic_mod: nml "generate_temp_ic"="generate_salt_ic"="generate_salt_ic"=false.  At least one should be true')

    !--- write out version information and namelist option ---------------
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
       write (stdout(), nml= idealized_ic_nml)
    endif

    !--- get grid information
    call get_grid

    !--- allocate memory to module variables
    allocate(temp(isc:iec,jsc:jec,nk), salt(isc:iec,jsc:jec,nk), age(isc:iec,jsc:jec,nk) )
    allocate(global_temp(ni,nj,nk), global_salt(ni,nj,nk), global_age(ni,nj,nk) )
    ! idealized temperature initial conditions 
    if(generate_temp_ic) then
       call idealized_temp_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,temp,global_temp)
       chksum_t = mpp_chksum(temp)
       write (stdout(),'(/(a,I20/))')' Idealized_ic Checksum: T=', chksum_t 
    endif

    ! idealized salinity initial conditions 
    if(generate_salt_ic) then 
       call idealized_salt_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,salt,global_salt)
       chksum_s = mpp_chksum(salt)
       write (stdout(),'(/(a,I20/))') ' Idealized_ic Checksum: S=', chksum_s
    endif

    ! idealized age initial conditions 
    if(generate_age_ic) then 
       call idealized_age_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,age,global_age)
       chksum_a = mpp_chksum(age)
       write (stdout(),'(/(a,I20/))') ' Idealized_ic Checksum: age=', chksum_a
    endif

    return

  end subroutine idealized_ic_init

  !#######################################################################
  !--- idealized temperature initial condition
  subroutine idealized_temp_ic

    real    :: t_value, fact
    integer :: i, j, k, kb

    select case (trim(temp_type))
    case('constant_temp_ic')
       t_value = constant_temp_value
       temp =  t_value
       write (stdout(),'(/a,(f6.2,a)/)') ' Note: constant temp IC is being set. T=',t_value,' deg C.'
    case('equatorial_temp_ic')
       write (stdout(),'(/a/)') ' Note: idealized equatorial temp IC is being set.'
       do k=1,nk
          if(z0 == 0.0) call mpp_error(FATAL,'idealized_ic_mod: nml "z0"= 0.0, but should be nonzero')
          t_value = t0*(1.0-tanh((zt(k)-thk)/z0)) + t1*(1.0-zt(k)/zt(nk))
          temp(:,:,k) = t_value
          write (stdout(),*) ' k=',k,' T=',t_value,' deg C.'
       enddo
    case('exponential_temp_ic')
       write (stdout(),'(/a/)') ' Note: exponential temperature IC is being set.'
       do k=1,nk
          t_value = 22.0*exp(-zt(k)/1000.0)
          temp(:,:,k) = t_value
          write (stdout(),*) ' k=',k,' T=',t_value,' deg C.'
       enddo
    case('shelfbowl_temp_ic')
       write (stdout(),'(/a/)') ' Note: shelfbowl temperature ic being set, assuming south-edge of shelf is at 70N.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                t_value = 5.0*(1.0 + 0.5*(1.0-tanh(0.5*(yt(i,j)-70.0))) )
                temp(i,j,k)  = t_value
             enddo
          enddo
       enddo
    case('zonal_levitus_temp_ic')
       write (stdout(),'(/a/)') ' Note: zonal averaged Levitus temperature IC is being set for depth < 2000m.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                temp(i,j,k)  = theta0 (yt(i,j), zt(k))
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(temp_type)//' is not a valid option of nml "temp_type" ')
    end select

    ! zero out temperature in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                temp(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

  end subroutine idealized_temp_ic

!#######################################################################
  ! idealized salinity initial conditions 
  subroutine idealized_salt_ic

    real    :: s_value, fact
    integer :: i,j,k,kb

    select case (trim(salt_type))
    case('constant_salt_ic')
       s_value = constant_salt_value
       salt = s_value
       write (stdout(),'(/a,2(f6.2,a)/)') ' Note: constant salinity IC being set. s=',s_value,' psu.'
    case('exponential_salt_ic')
       write (stdout(),'(/a/)') ' Note: exponential salinity IC is being set.'
       do k=1,nk
          s_value = 22.0*exp(-zt(k)/1000.0)
          salt(:,:,k) = s_value
          write (stdout(),*) ' k=',k,' S=',s_value, ' psu.'
       enddo
    case('salinity_profile_ic')
       write (stdout(),'(/a/)') ' Note: general salinity IC profile from salt0 is being set.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                salt(i,j,k)  = salt0 (zt(k))
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(salt_type)//' is not a valid option of nml "salt_type" ')
    end select

    ! zero out salinity in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                salt(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

end  subroutine idealized_salt_ic

!#######################################################################
  ! idealized age initial conditions 
  subroutine idealized_age_ic

    real    :: age_value, fact
    integer :: i,j,k,kb

    select case (trim(age_type))
    case('constant_age_ic')
       age_value = constant_age_value
       age = age_value
       write (stdout(),'(/a,2(f6.2,a)/)') ' Note: constant age IC being set. age=',age_value,' seconds.'
    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(age_type)//' is not a valid option of nml "age_type" ')
    end select

    ! zero out age in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                age(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

end  subroutine idealized_age_ic

!#######################################################################
  !--- read grid from grid_spec file.
  subroutine get_grid

    integer, dimension(2) :: layout  = (/ 1, 0 /)
    integer :: unit, npes, ndim, nvar, natt, ntime, len, i, j, k, kb, mk
    type(axistype), dimension(:), allocatable :: axes
    real, dimension(:,:), allocatable :: tmp, ht
    character(len=128) :: name

    !--- get grids information -------------------------------------------
    npes = mpp_npes()
    if(file_exist(trim(grid_file))) then
       call mpp_open(unit,trim(grid_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
       call mpp_get_info(unit, ndim, nvar, natt, ntime)
       allocate(axes(ndim))
       call mpp_get_axes(unit,axes)
       do i=1, ndim
          call mpp_get_atts(axes(i),name=name,len=len)
          select case(trim(name))
          case ('grid_x_T')
             ni = len
             allocate(grid_x_t(ni))
             call mpp_get_axis_data(axes(i),data=grid_x_t)
          case ('grid_y_T')
             nj = len
             allocate(grid_y_t(nj))
             call mpp_get_axis_data(axes(i),data=grid_y_t)        
          case ( 'zt' )
             nk = len
             allocate(zt(nk))
             call mpp_get_axis_data(axes(i),data=zt ) 
          end select
       enddo
       call mpp_close(unit)
       if(ni == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis grid_x_T')
       if(nj == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis grid_y_T')
       if(nk == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis zt')

       !--- define domain ------------------------------------------------
       call mpp_define_layout((/1,ni,1,nj/),npes,layout)
       call mpp_define_domains((/1,ni,1,nj/),layout, domain)
       call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
       allocate(xt(ni,nj), yt(ni,nj), kmt(isc:iec,jsc:jec),  &
            tmp(isc:iec,jsc:jec), ht(isc:iec,jsc:jec) )
       call read_data(trim(grid_file), 'x_T', xt)
       call read_data(trim(grid_file), 'y_T', yt)
       call read_data(trim(grid_file), 'depth_t', ht, domain)
       call read_data(trim(grid_file), 'num_levels', tmp,domain)
       kmt = tmp
       deallocate(tmp, axes)
    else
       call mpp_error(FATAL, 'idealized_ic_mod: file '//trim(grid_file)//' does not exist in your work directory')
    endif

  end subroutine get_grid

  !#######################################################################
  ! <SUBROUTINE NAME="write_idealized_ic_data">
  !
  ! <DESCRIPTION>
  !    Write out tracer data to netcdf file 
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine write_idealized_ic_data

    real, dimension(:,:,:,:),         allocatable :: tracer
    character(len=128), dimension(:), allocatable :: fld_name, fld_unit, fld_longname
    integer                                       :: n

    if(num_active_tracer .gt. 0) then
       allocate(tracer(ni,nj,nk,num_active_tracer), fld_name(num_active_tracer) )
       allocate( fld_unit(num_active_tracer), fld_longname(num_active_tracer) )
       n = 0
       if(generate_temp_ic) then
          n = n + 1
          tracer(:,:,:,n) = global_temp
          fld_name(n)     = 'temp'
          fld_unit(n)     = 'deg_c'
          fld_longname(n) = 'initial potential temp'
       endif
       if(generate_salt_ic) then
          n = n + 1
          tracer(:,:,:,n) = global_salt
          fld_name(n)     = 'salt'
          fld_unit(n)     = 'psu'
          fld_longname(n) = 'initial salinity'
       endif
       call write_tracer_data(active_tracer_file, tracer, fld_name, fld_unit, fld_longname )
       deallocate(tracer, fld_name, fld_unit, fld_longname)
    endif

    if(num_passive_tracer .gt. 0) then
       allocate(tracer(ni,nj,nk,num_passive_tracer), fld_name(num_passive_tracer) )
       allocate(fld_unit(num_passive_tracer), fld_longname(num_passive_tracer) )
       n = 0
       if(generate_age_ic) then
          n = n + 1
         tracer(:,:,:,n) = global_age
          fld_name(n)     = 'age_global'
          fld_unit(n)     = 'sec'
          fld_longname(n) = 'initial age'
       endif
       call write_tracer_data(passive_tracer_file, tracer, fld_name, fld_unit, fld_longname )
       deallocate(tracer, fld_name, fld_unit, fld_longname)
    endif

    return

  end subroutine write_idealized_ic_data


  !#######################################################################
  !--- write tracer data to netcdf file, including writing out meta.

  subroutine write_tracer_data(file, tracer, fld_name, fld_unit, fld_longname)

    character(len=*),               intent(in)   :: file
    real, dimension(:,:,:,:),    intent(inout)   :: tracer
    character(len=*), dimension(:), intent(in)   :: fld_name, fld_unit, fld_longname

    integer                                      :: unit, num_tracers, n
    type(axistype)                               :: axis_xt, axis_yt, axis_zt
    type(fieldtype)                              :: fld_x_T, fld_y_T
    type(fieldtype),dimension(size(fld_name(:))) :: fld_tracer

    num_tracers = size(fld_name(:))

    call mpp_open(unit,trim(file), MPP_OVERWR, MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)

    !--- write out meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_xt,'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = grid_x_t )
    call mpp_write_meta(unit, axis_yt,'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = grid_y_t )
    call mpp_write_meta(unit,axis_zt,'zt','meters','zt',&
         data=zt,cartesian='Z',sense=-1)
    call mpp_write_meta(unit, fld_x_T, (/axis_xt, axis_yt/), 'x_T', 'degree_east',  &
         'Geographic longitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_T, (/axis_xt, axis_yt/), 'y_T', 'degree_north',  &
         'Geographic latitude of T_cell centers', pack=1)
    do n = 1, num_tracers
       call mpp_write_meta(unit, fld_tracer(n), (/axis_xt, axis_yt, axis_zt/), trim(fld_name(n)), &
            trim(fld_unit(n)), trim(fld_longname(n)), pack=1)
    enddo
    !--- write out data --------------------------------------------------
    call mpp_write(unit, axis_xt)
    call mpp_write(unit, axis_yt)
    call mpp_write(unit, axis_zt)
    call mpp_write(unit, fld_x_T, xt)
    call mpp_write(unit, fld_y_T, yt)
    do n = 1, num_tracers    
       call mpp_write(unit, fld_tracer(n), tracer(:,:,:,n) )
    enddo
    call mpp_close(unit)

  end subroutine write_tracer_data

  !#######################################################################
 ! <SUBROUTINE NAME="idealized_ic_end">
  !
  ! <DESCRIPTION>
  !  Release memory.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine idealized_ic_end

    deallocate(temp, salt, age, xt, yt, kmt, zt, grid_x_t, grid_y_t)

    module_is_initialized = .false.
    return

  end subroutine idealized_ic_end


  !#######################################################################
  ! This subroutine returns estimates of global mean potential
  ! temperature for model initialization as a function of depth.
  ! it is used to produce a reference thermal stratification for the
  ! upper 2000m of the OCEAN`s test case.  below 2000m, the
  ! potential temperature returned is 2.0 degrees C.  surface
  ! values are set slightly above 18.4 degrees C at the reference
  ! latitude "reflat".
  ! the estimates are produced from a 7th order ploynomial fit to
  ! the annual mean world ocean potential temperature observations
  ! of Levitus (1982).
  !
  ! input [units]:
  !
  !   a latitdue (ydeg): [degrees] <BR/>
  !   a zt value (depth): [meters] <BR/>
  !
  ! output [units]:
  !
  !   potential temperature estimate (est): [degrees centigrade]
  !
  ! variables:
  !
  !   coeft     = coefficients for the polynomial fit of potential
  !               temperature vs. depth
  !
  !   reflat    = reference latitude at which observed surface
  !               temperatures approximately equal coeft(1)
  !
  !   factor    = the ratio of the cosine of the latitude requested
  !               ("ydeg") to the reference latitude ("reflat")
  !               used to scale the upper 2000 meters of the vertical
  !               temperature profile
  !
  !   tmin,tmax = the minumum and maximum potential temperatures
  !               allowed at the time of model initialization
  !
  ! reference:
  !   Levitus, S., Climatological atlas of the world ocean, NOAA
  ! Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.
  !
  ! Originally coded by Keith Dixon (Keith.Dixon@noaa.gov)
  function theta0 (ydeg, depth)

    real, intent(in)   :: ydeg, depth
    integer, parameter :: ndeg=7
    real, save, dimension(ndeg+1) :: coeft = (/ 0.184231944E+02,-0.430306621E-01, 0.607121504E-04,-0.523806281E-07,&
                                                0.272989082E-10,-0.833224666E-14,0.136974583E-17,-0.935923382E-22 /)
    real, save :: tmin=2.0, tmax=25.0, reflat=34.0
    real       :: refcos, coslat, factor, z, est, bb, theta0, theta, t0, t1, h, z0
    integer    :: nn

    refcos = abs(cos(reflat/RADIAN))

    coslat = abs(cos(ydeg/RADIAN))
    factor = coslat/refcos

    if (depth > 2000.) then
       est = 2.0
    else
       est = 0.0
       bb = 1.0
       do nn=1,ndeg+1
          est = est + coeft(nn)*bb
          bb = bb*depth
       enddo
       est = est * factor
    endif

    if (est > tmax) est = tmax
    if (est < tmin) est = tmin

    theta0 = est

  end function theta0

  !#######################################################################
  ! This function returns a salinity in psu
  function salt0 (depth)

    real, intent(in) :: depth
    real :: salt0

    salt0 = 32.0 + .001*depth 

  end function salt0

  !#######################################################################

end module idealized_ic_mod
