
FFLAGS  := -O3 -I$(NCEPLIBDIR)/incmod/sigio_4 
EXEC    := global_sighdr
LIBS    := -L$(NCEPLIBDIR)/lib -lsigio_4 -lw3lib-2.0_4

SOURCES := $(wildcard *.f90)
OBJECTS := $(SOURCES:%.f90=%.o)

.PHONY: all clean msg

all: msg $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o $(EXEC)

%.o : %.f90
	$(F90) $(FFLAGS) -c $<

clean: msg
	rm -rf *.o *.a *.mod $(EXEC)

msg:
	@echo ""
	@echo "============================================================"
	@echo " CFS (sfchdr)"
	@echo "============================================================"    
