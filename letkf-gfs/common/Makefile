#include ../configure.user
include ../../config/makefile.mk

#OBJS    = SFMT.o common.o common_mpi.o common_mtx.o netlib.o netlibblas.o common_letkf.o \
	  common_gfs.o common_gfs_pres.o common_obs_gfs.o common_mpi_gfs.o common_precip.o

libcommon := $(CFS_LETKF_ROOT)/common/common_cfs.a
libcommongfs := common_gfs.a
INCLUDE := -I$(CFS_LETKF_ROOT)/common
SRCS    := $(wildcard *.f90) $(wildcard *.f)
OBJS    := $(SRCS:%.f90=%.o)
OBJS    := $(OBJS:%.f=%.o)
#OBJS    = common.o common_mpi.o 
FOPTS   := $(INCLUDE) #-O2 -mcmodel=mediu

all: $(libcommongfs)

$(libcommongfs) : $(OBJS)
	ar -rcs $(libcommongfs) $(OBJS)

%.o : %.f90
	$(F90) $(FOPTS) -c $<
%.o : %.f
	$(F77) $(FOPTS) -c $<


# enssize:        enssize.o SFMT.o common.o common_mtx.o netlib.o common_letkf.o
# 	$(SFC) $(FOPTS) -o enssize enssize.o SFMT.o common.o common_mtx.o netlib.o common_letkf.o \
# 	$(LIBBLAS) -L$(SP_LIB) -lsp_4
# datetime:	datetime.o SFMT.o common.o
# 	$(SFC) $(FOPTS) -o datetime datetime.o SFMT.o common.o

common_mpi.o : $(libcommon)
common_gfs.o : $(libcommon)
common_mpi_gfs.o : $(libcommon) common_mpi.o common_gfs.o common_obs_gfs.o
common_obs_gfs.o : $(libcommon) common_gfs.o
common_letkf.o : $(libcommon) common_mtx.o
common_mtx.o : $(libcommon)
common_precip.o : $(libcommon) common_gfs.o common_letkf.o

clean:
	rm -f *.o *.mod enssize datetime
