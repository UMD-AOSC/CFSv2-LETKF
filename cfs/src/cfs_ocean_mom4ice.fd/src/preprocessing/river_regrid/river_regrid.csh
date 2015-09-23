#!/bin/tcsh 
#=======================================================================
#      river_regrid : remap river network data.
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#  This program can remap river network data from spherical grid onto another
#  spherical grid. 
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.

  set echo
#  set DEBUG                                                 # uncomment this to debug your run with totalview
  set platform     = sgi                                    # A unique identifier for your platform
  set name         = "river_regrid"                         # name of the data file will be generated
  set root         = $cwd:h:h:h                             # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include                          # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/river_regrid.exe            # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform      # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/river_regrid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF" )                     # list of cpp #defines to be passed to the source files

# list the source code
  set CORE      = "$tooldir/{river_regrid.f90}"
  set UTILITIES = "$sharedir/{axis_utils,constants,fms,mpp,platform,memutils}"
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
     &river_regrid_nml
       river_input_file = '/home/z1l/fms/src/preprocessing/river_regrid/input/orig_runoff.rivers_init.nc'
       land_grid_file   = '/archive/kap/cmdt/restart/grid_spec.partial_areas.1deg.nc'/
!
#       land_grid_file   = '/archive/bw/fms/cold_start/bgrid/amip2/N45_OM2.grid_spec.nc' /

#  run the executable
  if ( $?DEBUG ) then
     totalview $executable:t > fms.out
  else
     $executable:t #>fms.out
     cat fms.out
  endif
  
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile.out $name.logfile.out

  unset echo
