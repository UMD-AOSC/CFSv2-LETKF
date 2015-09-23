#!/bin/tcsh 
#=======================================================================
#      grid_transer 
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp can be used to convert the old name convention grid specification
#  netcdf file to new name convention grid netcdf file. The user need to provide the
#  input grid file $old_grid and the program will generate the output grid file
#  $new_grid. The output file $new_grid contains both the new name convention grid
#  infofmation and all the fields in $old_grid.
#
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================

  set echo
# set DEBUG                                             # uncomment this to debug your run with totalview
  set platform     = sgi                                # A unique identifier for your platform
  set name         = "grid_transfer"                    # name of the grid file will be generated
  set root         = $cwd:h:h:h:h                       # The directory you created when you checkout
  set tooldir      = $cwd                               # directory of the tool
  set sharedir     = $root/src/shared                   # directory of the shared code.
  set includedir   = $root/include                      # fms include directory
  set workdir      = $tooldir/workdir                   # where the tool is run and output is produced
  set mkmf         = $root/bin/mkmf                     # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF" )                 # list of cpp #defines to be passed to the source files
  set executable   = $tooldir/exec/grid_transfer.exe    # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/grid_transfer.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif

# grid file to be converted
  set old_grid     = "/home/z1l/fms/src/preprocessing/generate_grids/grid_transfer/grid_spec_v7.nc"
                    # grid file you want to transfer to new format.

# list the source code
  
  set CORE         = "$tooldir/{grid_transfer.F90}"
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,platform,memutils}"
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
    &grid_transfer_nml
      old_grid = '$old_grid'
      new_grid = '$name.nc'
       /
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


