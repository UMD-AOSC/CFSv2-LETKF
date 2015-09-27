# CFSV2-LETKF
Strongly Coupled LETKF applied to the Climate Forecast System version 2 (CFSv2). Default settings use the T62/L64 version of GFS and the 0.5 degree MOM. Fix files must be obtained separately. Information about how to run most of the scripts can be shown by running with the `-h` flag, (e.g. `./util/get_cfsr -h`

## Compiling
### Config files
The first step in compiling is to setup the configuration files in the in the `config/` directory. The following files are needed:
* `makefile.mk` - specifies the compiler options used during the build process
* `setenv.sh` - specifies the environemnt and modules loaded, used during the build process as well as when running the model or scripts. user must run `source setenv.sh` in any new bash shell when run running or compiling the CFS.

These files need to exist by creating a symbolic link to one of the other existing machine depenent files (or by creating new ones using the exising files as a guide). E.g, on the Deepthought 2 computer, run the following from within `config/`:
```
ln -s makefile.dt2_intel.mk makefile.mk
ln -s makefile.dt2_intel.mk makefile.mk
```
### Building CFS
To build the CFS (assuming using a bash shell):
* setup the environment by running `source config/setenv.sh`
* from within `lib-ncep/` run `make`. This will produce several files within `lib-ncep/lib/` and several folders within `lib-ncep/incmod/`
* from within `cfs` run `make`

Check to make sure the following end files are produced in `cfs/bin/` (The first 3 being the main CFSv2 files, and the rest are supporting files used in post or pre processing)
* `cfs_ocean_mom4ice`  -  coupled MOM ocean model
* `global_fcst`  - Coupled GFS atmospheric model
* `cfs_mlc_coupler` - Model coupler
* `global_sfchdr`
* `global_sighdr`
* `global_cycle`
* `global_chgres`

Using `make` in parallel build mode (i.e. `make -j n`) should work for most of the components, the original make files for some of the components (especially the GFS) do not have the dependencies defined correctly, if there are build errors try running make without parallel mode.


## Preparing for experiments
CFSR data is needed as initial conditions for ensembles, verification, and as boundary conditions (ozone, etc) for model runs. To download one month of CFSR data and convert to T62 resolution, run:

```
./util/get_cfsr --start 20100101 --end 20100201
./util/chgres_cfsr --start 2010010100 --end 2010020100 --res 62
```

To download an initial conditions for the ocean, run:

```
./util/get_cfsr_ocnIC --start 20100101 --daily
```

