
NCEPLIBS := sfcio_4 bacio_4 sp_d w3lib-2.0_d
FFLAGS   := -O3 $(addprefix -I$(NCEPLIBDIR)/incmod/, $(NCEPLIBS))
EXEC     := global_cycle
LIBS     := -L$(NCEPLIBDIR)/lib $(addprefix -l,$(NCEPLIBS)) -openmp


SOURCES := $(wildcard *.f) $(wildcard *.f90)
OBJECTS := $(SOURCES:%.f=%.o)
OBJECTS := $(OBJECTS:%.f90=%.o)

.PHONY: all clean

all: msg $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o $(EXEC)

%.o : %.f
	$(F77) $(FFLAGS) -c $<
%.o : %.f90
	$(F90) $(FFLAGS) -c $<

cycle.o : machine.o 
clean: msg
	rm -rf *.o *.a *.mod $(EXEC)
msg:
	@echo ""
	@echo "============================================================"
	@echo " CFS (cycle)"
	@echo "============================================================"    
