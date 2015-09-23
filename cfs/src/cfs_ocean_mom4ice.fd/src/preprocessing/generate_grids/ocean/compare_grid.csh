#!/bin/tcsh 
#=======================================================================
#      compare_grid : compare topography of two grid files.
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp compares two grid descriptor files (generated via ocean_grid_generator) 
#  and creates a text file output listing line-by-line differences between the
#  two files. Output file format is the same as the grid_edits file used by
#  edit_grid.F90. These two compared files should have same grid size and same grid
#  resolution.
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================

  set echo
# set DEBUG                                            # uncomment this to debug your run with totalview
  set platform     = sgi                               # A unique identifier for your platform
  set name         = "compare_grid"                    # name of tool and name of the output
  set root         = $cwd:h:h:h:h                      # The directory you created when you checkout
  set tooldir      = $cwd                              # directory of the tool
  set sharedir     = $root/src/shared                  # directory of the shared code.
  set includedir   = $root/include                          # fms include directory
  set workdir      = $tooldir/workdir                  # where the tool is run and output is produced
  set executable   = $tooldir/exec/compare_grid.exe    # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/compare_grid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                    # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF" )                # list of cpp #defines to be passed to the source files

  # grid data will be compared
  set grid_file_1  = $workdir/ocean_grid.nc            # one of the grid file to be compared
  set grid_file_2  = $workdir/edit_grid.nc             # one of the grid file to be compared

# list the source code
  
  set CORE      = "$tooldir/compare_grid.f90"
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist   = "$CORE $UTILITIES"

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

#---------------------------------------------------------------------
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
    &compare_grid_nml
       grid_file_1   = '$grid_file_1' 
       grid_file_2   = '$grid_file_2'
       grid_edits    = 'grid_edits'
       mask_diff     = 'mask_diff'
       /
!

  if ( $?DEBUG ) then
     totalview $executable:t > fms.out
  else
     $executable:t >fms.out
  endif
  cat fms.out
#   --- rename ascii output file
  mv fms.out $name.fms.out
  mv logfile.out $name.logfile.out
  
  unset echo  
