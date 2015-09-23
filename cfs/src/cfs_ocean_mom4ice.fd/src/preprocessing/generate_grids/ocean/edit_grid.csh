#!/bin/tcsh 
#=======================================================================
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp can can edit the topography of input grid_spec file "orig_grid" 
#  according to the ascii input file "grid_edits". Then it will output the 
#  new grid_spec file "mod_grid". 
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================
  set echo
# set DEBUG                                            # uncomment this to debug your run with totalview
  set platform     = sgi                               # A unique identifier for your platform
  set name         = "edit_grid"                       # name of the grid file will be generated
  set root         = $cwd:h:h:h:h                      # The directory you created when you checkout
  set tooldir      = $cwd                              # directory of the tool
  set sharedir     = $root/src/shared                  # directory of the shared code.
  set includedir   = $root/include                     # fms include directory   
  set workdir      = $tooldir/workdir                  # where the tool is run and output is produced
  set executable   = $tooldir/exec/edit_grid.exe       # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/edit_grid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                    # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF" )                # list of cpp #defines to be passed to the source files

# grid to be edited and grid_edits test file
  set grid_edits = grid_edits.txt                      # text file to specify the edit region and new depth.
  set orig_grid  = $workdir/ocean_grid.nc         

# list the source code
  
  set CORE      = " $tooldir/{edit_grid.F90,topog.f90,grids_type.f90,grids_util.f90} "
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils}"
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

# if the grid_edits file does not exist, create one here.
  if( ! -f $grid_edits) then
     cat >$grid_edits <<EOF
     100:105, 50:65, 1000   #lon_start:lon_end, lat_start:lat_end, new height.
     65, 75, 200
     5,16, 0
EOF
  endif

# --- set up namelist  

  cat >input.nml <<!
    &edit_grid_nml
       orig_grid   = '$orig_grid' 
       mod_grid    = '$name.nc'
       grid_edits  = '$grid_edits'
       /
!

  if ( $?DEBUG ) then
     totalview $executable:t > fms.out
  else
     $executable:t >fms.out
  endif
  cat fms.out

# rename ascii output file
  mv fms.out $name.fms.out
  mv logfile.out $name.logfile.out

  unset echo
