
FFLAGS = -O3 $(F_STRICT) $(F_CNVTBE) -align dcommon
ARFLAGS = -ruv

SIZE_4 := 
SIZE_8 := $(F_R8) $(F_I8)
SIZE_D := $(F_R8)


AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f)
OBJECTS := $(SOURCES:%.f=%.o)



.PHONY: all clean

all: builddir libsp_4/libsp_4.a libsp_8/libsp_8.a libsp_d/libsp_d.a

builddir:
	@mkdir -p libsp_4
	@mkdir -p libsp_8
	@mkdir -p libsp_d


## 4 byte version of libip
OBJECTS_4 := $(addprefix libsp_4/,$(OBJECTS))
libsp_4/libsp_4.a: $(OBJECTS_4)
	$(AR)  $@ $(OBJECTS_4)
libsp_4/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_4) -c $< -o $@


## 8 byte version of libip
OBJECTS_8 := $(addprefix libsp_8/,$(OBJECTS))
libsp_8/libsp_8.a: $(OBJECTS_8)
	$(AR)  $@ $(OBJECTS_8)
libsp_8/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_8) -c $< -o $@


## double precision (default integer, 8 byte real) version of libip
OBJECTS_D := $(addprefix libsp_d/,$(OBJECTS))
libsp_d/libsp_d.a: $(OBJECTS_D)
	$(AR)  $@ $(OBJECTS_D)
libsp_d/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_D) -c $< -o $@



clean:
	rm -rf libsp_4
	rm -rf libsp_8
	rm -rf libsp_d
