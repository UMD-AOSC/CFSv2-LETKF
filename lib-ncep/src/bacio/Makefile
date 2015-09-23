
FFLAGS = -O3 $(F_STRICT) $(F_CNVTBE)
CFLAGS = -O3 -DUNDERSCORE
ARFLAGS = -ruv

SIZE_4 := 
SIZE_8 := $(FR8) $(FI8)


AR := ar $(ARFLAGS)

SOURCES := $(wildcard *.f) bacio.v1.4.c
OBJECTS := $(SOURCES:%.f=%.o)
OBJECTS := $(OBJECTS:%.c=%.o) 



.PHONY: all clean


all: builddir libbacio_4/libbacio_4.a libbacio_8/libbacio_8.a 

builddir:
	@mkdir -p libbacio_4
	@mkdir -p libbacio_8


## 4 byte version of libbacio
OBJECTS_4 := $(addprefix libbacio_4/,$(OBJECTS))
libbacio_4/libbacio_4.a: $(OBJECTS_4)
	$(AR) $@ $(OBJECTS_4)
libbacio_4/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_4) -c $< -o $@
libbacio_4/%.o : %.c  libbacio_4/clib.h
	$(CC) $(CFLAGS) -Ilibbacio_4 -c $< -o $@
libbacio_4/clib.h: clib4.h
	cp clib4.h libbacio_4/clib.h


## 8 byte version of libbacio
OBJECTS_8 := $(addprefix libbacio_8/,$(OBJECTS))
libbacio_8/libbacio_8.a: $(OBJECTS_8)
	$(AR) $@ $(OBJECTS_8)
libbacio_8/%.o : %.f
	$(FC) $(FFLAGS) $(SIZE_8) -c $< -o $@
libbacio_8/%.o : %.c  libbacio_8/clib.h
	$(CC) $(CFLAGS) -Ilibbacio_8  -c $< -o $@
libbacio_8/clib.h: clib8.h
	cp clib8.h libbacio_8/clib.h



clean:
	rm -rf libbacio_4
	rm -rf libbacio_8

