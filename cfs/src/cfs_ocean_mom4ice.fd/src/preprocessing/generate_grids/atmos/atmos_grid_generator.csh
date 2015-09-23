#!/bin/tcsh 
#=======================================================================
#      preprocessing : atmos_grid_generator
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp can be used to generate bgrid or spectral horizontal grid
#  for atmosphere model or land model. It creates Makefile, compile if necessary
#  and saves output. The generated data sets are atmos_grid.nc (default value of
#  namelist varible output_file of atmos_grid_generator_nml) and are saved in
#  the directory $output_dir/data. The standard output is stored at
#  $output_dir/ascii/.
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================
  set echo
# set DEBUG                                                 # uncomment this to debug your run with totalview
  set platform     = sgi                                    # A unique identifier for your platform
  set name         = "atmos_grid"                           # name of the grid file will be generated
  set root         = $cwd:h:h:h:h                           # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include                          # fms include directory  
  set atmtool      = $root/src/atmos_spectral/tools         # directory contains atmos_spectral tool code.
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/atmos_grid_generator.exe # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform      # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/atmos_grid_generator.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF" )        # list of cpp #defines to be passed to the source files

# list the source code
  set CORE         = " $tooldir/{atmos_grid.f90,atmos_grid_generator.f90} "
  set CORE         = " $CORE $atmtool"
  set UTILITIES    = "$root/src/shared/{axis_utils,constants,fms,mpp,fft,platform,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir "
  set srclist      =  ( $CORE $UTILITIES )

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
  
      
# compile the model code and create executable
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
     &atmos_grid_generator_nml
       output_file = '$name.nc'  /    
    &atmos_grid_nml
       grid_type = 'bgrid'
       num_lon=144, num_lat=90   /
!

#  run the executable
  if ( $?DEBUG ) then
     totalview $executable:t > fms.out     
  else
     $executable:t >fms.out
  endif
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile.out $name.logfile.out

  unset echo

