#!/bin/bash

## load compiler and other modules
########################################
module load intel/2013.1.039
#module load openmpi/intel/1.8.1
#module load netcdf/4.3.2
module load netcdf-fortran
#module load hdf5
module load python/2.7.8

export HDFLIB=/cell_root/software/hdf/1.8.13/intel/2013.1.039/intel/shared/sys/lib
export LD_LIBRARY_PATH=$NETCDF_LIBDIR:$NETCDF_FORTRAN_LIBDIR:$HDFLIB:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/ofed/lib64"
ulimit -s unlimited

## need python's netcdf library, not on DT2, have to use my own version
export PYTHONPATH=/lustre/tsluka/software/python-netCDF4/lib:$PYTHONPATH

################################################################################
## determine the root directory of the CFS-LETKF,
## no need to modify this...
SOURCE="${BASH_SOURCE[0]}"
export CFS_LETKF_ROOT="$( cd -P "$( dirname "$SOURCE" )/.." && pwd )"
################################################################################

## set directories
########################################
export TMP_DIR_LOCAL=/dev/shm/$USER
export TMP_DIR_SHARED=/lustre/tsluka/tmp/cfs/

export FIX_DIR_AM=/lustre/tsluka/CFSv2-LETKF/support/fix/fix_am
export FIX_DIR_OM=/lustre/tsluka/CFSv2-LETKF/support/fix/fix_om
export CFSR_DIR=$CFS_LETKF_ROOT/DATA/CFSR
export OBS_ATM=$CFS_LETKF_ROOT/DATA/obs/atm_prepbufr
export OBSNCEP=$CFS_LETKF_ROOT/DATA

## other properties when running
export NPROC_OM=19      ## number of cores for ocean model
export NPROC_AM=20      ## number of cores for atmosphere model
                        ## plus 1 core for the coupler

export NPROC_LETKF=20   ## number of cores for the LETKF
