.PHONY : all clean

platform   := cfg
name       := cfs_ocean_mom4ice
root       := $(shell pwd)
curdir     := $(root)/compile

code_dir   := $(root)/src
executable := $(curdir)/exec/$(name)
time_stamp := $(root)/bin/time_stamp.csh
mkmfTemplate := $(root)/bin/mkmf.template.$(platform)
mkmf       := $(root)/bin/mkmf
pathnames  := $(curdir)/path_names
cppDefs    := "-DENABLE_GDS -Duse_netCDF -Duse_libMPI"

mfile := $(executable)_h/Makefile



all : $(mfile)
	cd $(executable)_h ; make


$(mfile):
	mkdir -p $(executable)_h
	cd $(executable)_h ; $(mkmf) -f -a $(code_dir) -t $(mkmfTemplate) -p $(executable) -c $(cppDefs) $(pathnames) $(root)/include $(code_dir)/shared/mpp/include


clean:
	rm -rf $(curdir)/exec
