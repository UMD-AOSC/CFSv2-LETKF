#!/bin/csh -f
#
# This is a CSH program to test files used by time_interp_external_mod
# and/or data_override_mod.  This is intended to help debug runtime errors
# in FMS based codes where the above modules fail reading an INPUT file. 
# This is a quick way to make sure the file works without having to run the
# entire model.
#
# USER INPUT:
# filename = name of file to be tested
# fieldname = name of field in file
# year0 = initial year to use (typically start time of model)
# month0 = initial month to use
# day0 = initial day to use
# days_inc = increment in days
# ntime = number of time steps
# cal_type = calendar ('julian' , 'no_leap' or '360_day')
#
#
# Compilation:
#$MKMF -m Makefile -p test_time_interp_ext.exe -t $TEMPLATE -c -Dtest_time_interp_external -x  shared/{time_manager,fms,mpp,clocks,time_interp,axis_utils,platform,horiz_interp,constants,memutils} 

set NPES = 1
set FMS_BIN_DIR = /usr/local/bin

cat > input.nml <<EOF
 &test_time_interp_ext_nml
 filename='foo.nc'
 fieldname='foo'
 year0=1900
 month0=1
 day0=1
 days_inc=1
 ntime=1
 cal_type='julian'
 /
EOF

mpirun -np $NPES $FMS_BIN_DIR/test_time_interp_ext.exe
