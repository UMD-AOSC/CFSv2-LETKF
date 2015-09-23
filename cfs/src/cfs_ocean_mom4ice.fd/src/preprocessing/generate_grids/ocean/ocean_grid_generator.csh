#!/bin/tcsh -f
#=======================================================================
#      preprocessing : ocean_grid_generator
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp can be used to generate horizontal grid and/or vertical grid
#  and/or topography. It creates Makefile, compile if necessary and  saves output.
#  The generated data sets are grid_spec.nc (default value of namelist varible
#  output_file of grid_generator_nml) and are saved in the directory $work_dir.
#  When you want to generate topography (not idealized or special topography), a
#  topography source data file is needed, which is specified by the namelist
#  variable topog_file of topog_nml.
#  *******Important notice: You need to configure hgrid_nml, vgrid_nml, topog_nml
#  and grid_generator_nml, which are listed in this runscripts, to obtain desired grid.
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================
  set echo
# set DEBUG                                                 # uncomment this to debug your run with totalview
  set platform     = sgi                                    # A unique identifier for your platform
  set name         = "ocean_grid"                           # name of the grid file will be generated
  set npes         = 1                                      # number of processors
  set root         = $cwd:h:h:h:h                           # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include                          # fms include directory  
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/ocean_grid_generator.exe # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform      # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/ocean_grid_generator.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_libMPI" )        # list of cpp #defines to be passed to the source files
  if($npes == 1) set cppDefs  = ( "-Duse_netCDF" )

# list the source code
  
  set CORE         = "$tooldir/{grids_type.f90,grids_util.f90,ocean_grid_generator.f90}"
  set CORE         = "$CORE $tooldir/{hgrid.f90,vgrid.f90,topog.f90}"
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist      =  ( $CORE $UTILITIES )                       # list of the source code

# the following compiler choice are for GFDL user only. If the following command
# does not apply to you, comments or remove them
  if( $platform == "sgi" ) then
     source /opt/modules/modules/init/csh
     module switch mpt mpt_1900  
     module switch mipspro mipspro_741n
  endif
  if( $platform == "ia64" ) then
     source /opt/modules/default/init/tcsh
     module switch ifort.8.0.044 ifort.8.1.023
     setenv MALLOC_CHECK_ 0
  endif


# compile the code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist /usr/local/include

  make $executable:t

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist

    cat >input.nml <<!
     &ocean_grid_generator_nml
       grid_type   = 'hgrid_vgrid_topog'
       output_file = '$name.nc'  /    
    &hgrid_nml
       nxlons=2,x_lon=0.,360.,dx_lon=2.5,2.5,
       nylats=2,y_lat=-90.,90.,dy_lat=2.0,2.0,
       tripolar_grid=.false.,lat_join=65, 
       debug = .true. /
    &vgrid_nml
       nzdepths=3,z_depth=0.0,100.0,5600.0,dz_depth=25.0,25.0,975.0
       debug = .true. /
    &topog_nml
       topography = 'rectangular_basin' 
       topog_depend_on_vgrid = .TRUE.
       topog_file='$root/data/topog/OCCAM_p5degree.nc',
       topog_field='topo',
       fill_first_row = .true.
       kmt_min=2,
       filter_topog=.false.,
       num_filter_pass=5,
       scale_factor=-1,
       interp_method = "spherical"
       debug = .true. /
    &fms_nml
       domains_stack_size = 103680 /
!

#  run the executable
  if ( $?DEBUG ) then
     if($npes == 1) then
         totalview $executable:t > fms.out
     else
         totalview mpirun -a -np $npes $executable:t > fms.out
     endif
  else
     if($npes == 1) then
        $executable:t >fms.out
     else
        mpirun -np $npes $executable:t >fms.out
     endif
  endif

  cat fms.out     

#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile.out $name.logfile.out

  unset echo

