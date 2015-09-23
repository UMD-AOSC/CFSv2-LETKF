
LIBINC  := -I../../incmod

FFLAGS  := -O3 $(F_STRICT) $(F_CNVTBE)
CFLAGS  := -O3 -DLINUX -Dfunder -DUNDERSCORE -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'
ARFLAGS := -ruv

SIZE_4 := 
SIZE_8 := $(F_R8) $(F_I8)
SIZE_D := $(F_R8)


AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f) $(wildcard *.c)
OBJECTS := $(SOURCES:%.f=%.o)
OBJECTS := $(OBJECTS:%.c=%.o)



all: builddir lib_4/libw3lib-2.0_4.a  lib_8/libw3lib-2.0_8.a lib_d/libw3lib-2.0_d.a

builddir:
	@mkdir -p lib_4
	@mkdir -p lib_8
	@mkdir -p lib_d


## 4 byte version
OBJECTS_4 := $(addprefix lib_4/,$(OBJECTS))
lib_4/libw3lib-2.0_4.a: $(OBJECTS_4)
	$(AR)  $@ $(OBJECTS_4)
lib_4/%.o : %.f
	$(FC) $(LIBINC)/sigio_4 $(FFLAGS) $(SIZE_4) -c $< -o $@
lib_4/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

## 8 byte version
OBJECTS_8 := $(addprefix lib_8/,$(OBJECTS))
lib_8/libw3lib-2.0_8.a: $(OBJECTS_8)
	$(AR)  $@ $(OBJECTS_4)
lib_8/%.o : %.f
	$(FC) $(LIBINC)/sigio_4 $(FFLAGS) $(SIZE_8) -c $< -o $@
lib_8/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

## double precision version
OBJECTS_D := $(addprefix lib_d/,$(OBJECTS))
lib_d/libw3lib-2.0_d.a: $(OBJECTS_D)
	$(AR)  $@ $(OBJECTS_D)
lib_d/%.o : %.f
	$(FC) $(LIBINC)/sigio_4 $(FFLAGS) $(SIZE_D) -c $< -o $@
lib_d/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -rf lib_?
