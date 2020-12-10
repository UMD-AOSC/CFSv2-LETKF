#!/bin/bash

#NETCDF_DIR=/cell_root/software/netcdf/4.3.2/intel/2013.1.039/openmpi/1.8.1/hdf5/1.8.13/hdf4/4.2.10/sys
#$NETCDFF_DIR=/cell_root/software/netcdf-fortran/4.4.1/netcdf/4.3.2/intel/2013.1.039/openmpi/1.8.1/sys
#NETCDF_LIB="-L$NETCDFF_DIR/lib -lnetcdff -Wl,-rpath,$NETCDFF_DIR/lib"
#NETCDF_INC="-I$NETCDFF_DIR/include -I$NETCDF_DIR/include"

NETCDF_ROOT=/gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.1.3
NETCDF_FORTRAN_ROOT=/gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.1.3
NETCDF_DIR=$NETCDF_ROOT
NETCDFF_DIR=$NETCDF_FORTRAN_ROOT
NETCDF_LIB="-L$NETCDF_DIR/lib -L$NETCDFF_DIR/lib -lnetcdff -lnetcdf"
NETCDF_INC="-I$NETCDFF_DIR/include -I$NETCDF_DIR/include -I/gpfs1/home/Libs/INTEL/ZLIB/zlib-1.2.8/include -I/gpfs1/home/Libs/INTEL/SZIP/szip-2.1/include -I/gpfs1/home/Libs/INTEL/HDF5/hdf5-1.8.11_new/include"
