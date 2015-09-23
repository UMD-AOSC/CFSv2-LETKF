#!/bin/tcsh 
#=======================================================================
#      regrid_3d : interp 3-D data
#         Contact : Zhi Liang   email : Zhi.Liang@noaa.gov
#
#  This runscritp can be used to regrid 3-d lat-lon gridded data to
#  logically rectangular grid described by grid descriptor file $dest_file.
#  Applies only to scalar fields.
#
#  This preprocessing program was tested on the sgi origin3000 system at GFDL.
#  In order to run on other system, some changes may be needed.
#=======================================================================

  set echo
# set DEBUG                                                 # uncomment this to debug your run with totalview
  set platform     = sgi                                    # A unique identifier for your platform
  set name         = "regrid_3d"                            # name of the data file will be generated
  set npes         = 1                                      # number of processors
  set root         = $cwd:h:h:h                             # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include                          # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/regrid_3d.exe            # executable created after compilation
  set mkmfTemplate = $root/bin/mkmf.template.$platform      # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/regrid_3d.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_libMPI" )        # list of cpp #defines to be passed to the source files
  if($npes == 1) set cppDefs  = ( "-Duse_netCDF" )

# input data file and destination grid
  set src_file        = /archive/fms/mom4/input_data/levitus_ewg.nc                     # source data file
  set dest_grid       = $tooldir:h/generate_grids/ocean/workdir/ocean_grid.nc # destination grid

# list the source code
  set CORE      = "$tooldir/{regrid_3d.f90}"
  set UTILITIES = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils}"
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
    &regrid_3d_nml
       src_file = '$src_file'
       numfields = 2
       src_field_name = 'temp', 'salt'
       dest_grid = '$dest_grid'
       dest_file = '$name.nc'
       num_nbrs = 4
       /
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

