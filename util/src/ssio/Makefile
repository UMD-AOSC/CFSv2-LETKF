
PROGS   = ss2grd grd2ss sscycle ss2grdp grdctl

OBJS    = ssio_tools.o physcons.o

COMMON_DIR = ./

INCS = -I$(NCEPLIBDIR)/incmod/sigio_4 -I$(NCEPLIBDIR)/incmod/sfcio_4
LIBS = -L$(NCEPLIBDIR)/lib -lsigio_4 -lbacio_4 -lsfcio_4 -lsp_4

FOPTS  = -O3 $(F_CNVTBE)
FOPTS2 = -O3

.PHONY: all clean

all:	$(PROGS)

commonlibs= common_gfs_pres.o common_gfs.o ssio_tools.o
ss2grd.o : common_gfs.o ssio_tools.o
common_gfs.o : common.o
ssio_tools.o: common_gfs_pres.o physcons.o


% : %.o $(commonlibs)
	$(FC) $(FOPTS) $^ $(LIBS) -o $@

%.o : %.f90
	$(FC) $(FOPTS) -c $<  -I$(COMMON_DIR) $(INCS) -o $@


clean:
	rm -f *.o *.mod $(PROGS)
