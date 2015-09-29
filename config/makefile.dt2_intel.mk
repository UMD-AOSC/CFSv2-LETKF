export


## compilers
FC  = mpif90
F77 = $(FC)
F90 = $(FC)
CC  = mpicc

## path to the ESMF (earth system modeling framework) version v3.1.0rp5
## make sure that ESMF is compiled with the matching MPI settings for the
## MPI used to compile GFS (openmpi vs mpich, etc)
ESMFDIR   = /lustre/tsluka/CFSv2/external/esmf
ESMFMODS  = $(ESMFDIR)/mod/modO/Linux.intel.64.openmpi.default
ESMFLIBS  = $(ESMFDIR)/lib/libO/Linux.intel.64.openmpi.default

## path to NCEP libraries
NCEPLIBDIR = $(CFS_LETKF_ROOT)/lib-ncep

## paths to NETCDF (first 2 already defined by our environment)
# NETCDF_FORTRAN_INCDIR =
# NETCDF_INCDIR =
NETCDF_FORTRAN_LIBS = -L$(NETCDF_FORTRAN_LIBDIR) -L$(NETCDF_LIBDIR) -L$(HDF5LIB) -lnetcdf -lnetcdff -lhdf5
NETCDF_LIBS = -lnetcdf

## fortan compiler options dependant on compiler used
F_CNVTBE  = -convert big_endian  # big endian bit format
F_FREE    = -FR  # freeform fortran 90
F_STRICT  = -fp-model strict  # f77 style file format
F_R8 = -r8
F_R4 = -r4
F_I8 = -i8
F_I4 = -i4
FMCMODEL = -mcmodel=medium
LIBBLAS = -mkl
