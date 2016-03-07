#!/bin/bash

## ------------------------------------------------------------
## Define environment variables that are used throughough the
## CFSv2-LETKF system
## ------------------------------------------------------------

## load compiler and other modules
########################################
module load intel/2013.1.039
module load netcdf-fortran
module load python/2.7.8
#module load netcdf/4.3.2
#module load hdf5

## compiler
export F90=mpiifort
export FC=$F90
export F77=$F90

## other modules
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
## TMP_DIR_LOCAL is used for fast file IO, existing only on the remote
##  computational node
## TMP_DIR_SHARED must be accessible from multiple nodes
export TMP_DIR_LOCAL=/dev/shm/$USER
export TMP_DIR_SHARED=/lustre/tsluka/tmp/cfs/

## Fix files for the CFSv2, this must be optained elsewhere
export FIX_DIR_AM=/lustre/tsluka/CFSv2-LETKF/support/fix/fix_am
export FIX_DIR_OM=/lustre/tsluka/CFSv2-LETKF/support/fix/fix_om

##default location for data that is used in some of the scripts
export CFSR_DIR=$CFS_LETKF_ROOT/DATA/CFSR
export OBS_ATM=$CFS_LETKF_ROOT/DATA/obs/atm_prepbufr #TODO, remove this?
export OBSNCEP=$CFS_LETKF_ROOT/DATA  #TODO, remove this?


## Cluster configuration
################################
## The number of processors on each node
export NPROC_PERNODE=20

## number of cores for ocean model (OM) atmosphere model (AM) pluse
## one core (not defined) to be used for the coupler
export NPROC_OM=19        
export NPROC_AM=20        

## number of cores for the LETKF for the atmosphere (A) and ocean (O)
##  Should be a multiple of NPROC_PERNODE for best performance
export NPROC_LETKF_A=60
export NPROC_LETKF_O=120



## Append the path, sometimes it's usefule to have the
##  global_sighdr and such tools added to the path
export PATH=$CFS_LETKF_ROOT/cfs/bin:$CFS_LETKF_ROOT/util/bin:$PATH
