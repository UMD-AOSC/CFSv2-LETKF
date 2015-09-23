
FFLAGS  := -O3 $(F_STRICT) $(F_CNVTBE)
ARFLAGS := -ruv

SIZE_4 := 
SIZE_8 := $(F_R8) $(F_I8)
SIZE_D := $(F_R8)


AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f)
OBJECTS := $(SOURCES:%.f=%.o)



.PHONY: all clean

all: builddir libip_4/libip_4.a libip_8/libip_8.a libip_d/libip_d.a

builddir:
	@mkdir -p libip_4
	@mkdir -p libip_8
	@mkdir -p libip_d


## 4 byte version of libip
OBJECTS_4 := $(addprefix libip_4/,$(OBJECTS))
libip_4/libip_4.a: $(OBJECTS_4)
	$(AR)  $@ $(OBJECTS_4)
libip_4/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_4) -c $< -o $@


## 8 byte version of libip
OBJECTS_8 := $(addprefix libip_8/,$(OBJECTS))
libip_8/libip_8.a: $(OBJECTS_8)
	$(AR)  $@ $(OBJECTS_8)
libip_8/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_8) -c $< -o $@


## double precision (default integer, 8 byte real) version of libip
OBJECTS_D := $(addprefix libip_d/,$(OBJECTS))
libip_d/libip_d.a: $(OBJECTS_D)
	$(AR)  $@ $(OBJECTS_D)
libip_d/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_D) -c $< -o $@



clean:
	rm -rf libip_4
	rm -rf libip_8
	rm -rf libip_d
