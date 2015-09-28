#!/bin/bash

## load compiler and other modules
########################################
module load intel/2013.1.039
module load openmpi/intel/1.8.1
module load netcdf/4.3.2
module load netcdf-fortran
module load hdf5
module load python/2.7.8
export LD_LIBRARY_PATH=$NETCDF_LIBDIR:$NETCDF_FORTRAN_LIBDIR:$HDFLIB:$LD_LIBRARY_PATH
ulimit -s unlimited


################################################################################
## determine the root directory of the CFS-LETKF,
## no need to modify this...
SOURCE="${BASH_SOURCE[0]}"
export CFS_LETKF_ROOT="$( cd -P "$( dirname "$SOURCE" )/.." && pwd )"
################################################################################

## set directories
########################################
export TMP_DIR_LOCAL=/dev/shm/$USER

export FIX_DIR_AM=/lustre/tsluka/CFSv2/support/fix/fix_am
export FIX_DIR_OM=/lustre/tsluka/CFSv2/support/fix/fix_om
export CFSR_DIR=$CFS_LETKF_ROOT/DATA/CFSR


## other properties when running
export NPROC_OM=8
export NPROC_AM=11
