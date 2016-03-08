
FFLAGS  := -O3 -I$(NCEPLIBDIR)/incmod/sfcio_4 
EXEC    := global_sfchdr
LIBS    := -L$(NCEPLIBDIR)/lib -lsfcio_4 -lw3lib-2.0_4

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
