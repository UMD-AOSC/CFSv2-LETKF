#!/bin/tcsh 
#=======================================================================
#      idealized_ic : generate idealized initial condition
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscript can be used to generate tracer initial conditions for mom4.
#  It creates Makefile, compiles if necessary and  saves output.The generated
#  data sets are saved in the directory $output_dir/data. A grid_spec file
#  is needed for this program. 
#  ******Important notice: to obtain desired initial condition, you need to
#  specify initial_nml which is listed in this runscript.
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================

  set echo
# set DEBUG                                                 # uncomment this to debug your run with totalview
  set platform     = sgi                                    # A unique identifier for your platform
  set name         = "idealized_ic"                         # name of the tool
  set npes         = 1                                      # number of processors
  set root         = $cwd:h:h:h:h                           # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include                          # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/idealized_ic.exe            # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform      # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/idealized_ic.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_libMPI" )        # list of cpp #defines to be passed to the source files
  if($npes == 1) set cppDefs  = ( "-Duse_netCDF" )

# destination grid
  set grid_spec_file  = $tooldir:h:h/generate_grids/ocean/workdir/ocean_grid.nc  # destination grid

# list the source code
  set CORE      = $tooldir/{idealized_ic.f90,idealized_ic_driver.f90}
  set UTILITIES = "$sharedir/{fms,mpp,platform,constants,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist   = ( $CORE $UTILITIES )

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
    &idealized_ic_nml
       temp_type = 'zonal_levitus_temp_ic',
       salt_type = 'salinity_profile_ic'
       constant_temp_value = 15.1 ,
       constant_salt_value = 32. 
       passive_tracer_file = 'ocean_age.res.nc'
       grid_file           = '$grid_spec_file' /
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

