

FFLAGS  := -O2 $(F_FREE) $(F_CNVTBE)
FFLAGB  := -O2 $(F_CNVTBE)

ARFLAGS := -ruv
AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f)
OBJECTS := $(SOURCES:%.f=%.o)

.PHONY: all clean

all: libsigio_4.a

libsigio_4.a: $(OBJECTS)
	$(AR) $@ $(OBJECTS)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

bafrio.o : bafrio.f
	$(FC) $(FFLAGB) -c $< -o $@

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.a
