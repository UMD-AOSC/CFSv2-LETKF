
FFLAGS  := -O3 $(F_FREE) $(F_CNVTBE)
ARFLAGS := -ruv

SIZE_4 := 
SIZE_8 := $(F_R8) $(F_I8)
SIZE_D := $(F_R8)


AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f)
OBJECTS := $(SOURCES:%.f=%.o)

LIB := liblandsfcutil

.PHONY: all clean

all: builddir lib_4/$(LIB)_4.a lib_8/$(LIB)_8.a lib_d/$(LIB)_d.a

builddir:
	@mkdir -p lib_4
	@mkdir -p lib_8
	@mkdir -p lib_d

## 4 byte version of lib
OBJECTS_4 := $(addprefix lib_4/,$(OBJECTS))
lib_4/$(LIB)_4.a: $(OBJECTS_4)
	$(AR)  $@ $(OBJECTS_4)
lib_4/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_4) -c $< -o $@ -module lib_4
lib_4/interp_utils.o : lib_4/ll2xy_utils.o
lib_4/interp_utils_nesdis.o : lib_4/ll2xy_utils.o


## 8 byte version of lib
OBJECTS_8 := $(addprefix lib_8/,$(OBJECTS))
lib_8/$(LIB)_8.a: $(OBJECTS_8)
	$(AR)  $@ $(OBJECTS_8)
lib_8/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_8) -c $< -o $@ -module lib_8
lib_8/interp_utils.o : lib_8/ll2xy_utils.o
lib_8/interp_utils_nesdis.o : lib_8/ll2xy_utils.o


## d byte version of lib
OBJECTS_D := $(addprefix lib_d/,$(OBJECTS))
lib_d/$(LIB)_d.a: $(OBJECTS_D)
	$(AR)  $@ $(OBJECTS_D)
lib_d/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_D) -c $< -o $@ -module lib_d
lib_d/interp_utils.o : lib_d/ll2xy_utils.o
lib_8/interp_utils_nesdis.o : lib_8/ll2xy_utils.o

clean:
	rm -rf lib_4
	rm -rf lib_8
	rm -rf lib_d
