
FFLAGS  := -O3 $(F_R8) $(F_CNVTBE) -openmp -I$(NCEPLIBDIR)/incmod/sigio_4 -I$(NCEPLIBDIR)/incmod/gfsio_4 -I$(NCEPLIBDIR)/incmod/sfcio_4 -I$(NCEPLIBDIR)/incmod/landsfcutil_d
EXEC    := global_chgres
LIBS    := -L$(NCEPLIBDIR)/lib -lsfcio_4 -lgfsio_4 -lsigio_4 -llandsfcutil_d -lbacio_4 -lw3lib-2.0_d -lip_d -lsp_d -openmp

SOURCES := $(wildcard *.f) $(wildcard *.f90)
OBJECTS := $(SOURCES:%.f=%.o)
OBJECTS := $(OBJECTS:%.f90=%.o)

.PHONY: all clean

all:msg $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o $(EXEC)

%.o : %.f
	$(F77) $(FFLAGS) -c $<
%.o : %.f90
	$(F90) $(FFLAGS) -c $<

chgres.o   : nstio_module.o gfsio_rst.o funcphys.o surface_chgres.o
funcphys.o : physcons.o
physcons.o : machine_8.o

clean:msg
	rm -rf *.o *.a *.mod $(EXEC)
msg:
	@echo ""
	@echo "============================================================"
	@echo " CFS (chgres)"
	@echo "============================================================"    
