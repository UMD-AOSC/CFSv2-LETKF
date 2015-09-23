#!/bin/tcsh 
#=======================================================================
#      regrid : interp 2-D/3-D logical rectangular grid to logical rectangular grid.
#         Contact : Zhi Liang   email : z1l@gfdl.noaa.gov
#
#  This runscritp can be used to regrid data $src_data on any grid $src_grid to logically
#  rectangular grid described by grid descriptor file $dst_grid. 
#
#  This preprocessing program was tested on the sgi origin3000 and altix system at GFDL.
#  In order to run on other system, some changes may be needed.

  set echo
#  set DEBUG                                                # uncomment this to debug your run with totalview
  set platform     = ia64                                   # A unique identifier for your platform
  set name         = "regrid"                               # name of the data file will be generated
  set npes         = 24                                     # number of processors
  set root         = $cwd:h:h:h                             # The directory you created when you checkout
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $root/include
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/regrid.exe               # executable created after compilation
  set mkmfTemplate = /home/fms/bin/mkmf.template.$platform  # path to template for your platform
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/regrid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = /home/fms/bin/mkmf                     # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_libMPI" )        # list of cpp #defines to be passed to the source files
  if($npes == 1) set cppDefs  = ( "-Duse_netCDF" )

# input data file and destination grid
  set src_data        =  /archive/z1l/regrid_data_set/ocean.000101-000512.uv.nc
  set src_grid        = /archive/ms/cm2/input/grid_spec_v5.nc
  set dst_grid       = /archive/z1l/grids/xgrid/N45_new_grid.nc

# list the source code
  set CORE      = "$tooldir/{regrid.F90}"
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

  if($status != 0) then
     unset echo
     echo " Error in compilation "
     exit 1
  endif

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist

   cat >input.nml <<!
   &regrid_nml
      src_data = '$src_data',
      src_grid = '$src_grid',
      dst_grid = '$dst_grid',
      dst_data = '$name.nc',
      num_flds = 2
      fld_name = 'u','v'
      fld_pos  =  'C','C'
      vector_fld = .true., .true.
      use_source_vertical_grid = .false.
      apply_mask = .true.
      debug      = .false. /
!

#  run the executable
  if ( $?DEBUG ) then
     if($npes == 1) then
         totalview $executable:t #> fms.out
     else
         totalview mpirun -a -np $npes $executable:t #> fms.out
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
