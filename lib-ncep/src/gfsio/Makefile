LIB := libgfsio_4.a

FFLAGS  := -O3 $(F_STRICT) $(F_CNVTBE)
AR := ar
ARFLAGS := -ruv

SOURCES := $(wildcard *.f)  $(wildcard *.f90)
OBJECTS := $(SOURCES:%.f=%.o)
OBJECTS := $(OBJECTS:%.f90=%.o)

.PHONY: all clean

all: $(LIB)

$(LIB) : $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
%.o : %.f
	$(F77) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(F90) $(FFLAGS) -c $< -o $@

clean:
	rm -rf *.o
	rm -rf $(LIB)
