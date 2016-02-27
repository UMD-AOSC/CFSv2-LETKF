include $(CFS_LETKF_ROOT)/config/makefile.mk

#COMMON_DIR = $(CFS_LETKF_ROOT)/letkf-gfs/common

# COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
#                $(COMMON_DIR)/common.o         \
#                $(COMMON_DIR)/common_mtx.o     \
#                $(COMMON_DIR)/netlib.o         \
#                $(COMMON_DIR)/common_letkf.o   \
#                $(COMMON_DIR)/common_gfs.o     \
#                $(COMMON_DIR)/common_obs_gfs.o


#LIBS := -L$(NCEPLIBDIR)/lib -lsp_4 $(LIBBLAS) -lbufrlib
LIBS := $(LIBBLAS) -L$(NCEPLIBDIR)/lib -lsp_4 -lbufrlib

#$(CFS_LETKF_ROOT)/common/common_cfs.a

#INCLUDE := -I$(CFS_LETKF_ROOT)/common
FOPTS= -O2
#$(INCLUDE) #-mcmodel=medium

.PHONY: all clean

all: dec_prepbufr

dec_prepbufr:dec_prepbufr.o common.o common_obs.o SFMT.o
	$(FC) $(FOPTS) $(COMMON_OBJS) $? -o $@ $(LIBS) 

dec_prepbufr.o: common.o common_obs.o
common.o:SFMT.o

%.o:%.f90
	$(F90) $(FOPTS) -c $<

clean:
	rm -f *.o *.mod dec_prepbufr
