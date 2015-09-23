LIB := libsfcio_4.a

FFLAGS  := -O2 $(F_FREE) $(F_CNVTBE)
AR := ar
ARFLAGS := -ruv

SOURCES := $(wildcard *.f90)
OBJECTS := $(SOURCES:%.f90=%.o)

.PHONY: all clean

all: $(LIB)

$(LIB) : $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
%.o : %.f90
	$(F77) $(FFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.mod
	rm -rf $(LIB)
