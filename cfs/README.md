This directory contains the CFSv2 code from NCEP. 

See the wiki entry on how to compile the CFS: (Compiling CFSv2)[https://github.com/UMD-AOSC/CFSv2-LETKF/wiki/Compiling]

####Directory Structure
Model Components
* `src/cfs_global_fcst.fd` - GFS atmospheric model
* `src/cfs_mlc_coupled.fd` - CFSv2 coupler
* `src/cfs_ocean_mom4ice.fd` - MOM ocean model

Utility programs
* `global_chgres.fd` - converts sig/sfc GFS files to different resolutions
* `global_cycle.fd`  - converts a GFS forecast into an analysis that can be used as initial conditions for a new model run.
* `global_sfchdr.fd` - prints header information about a GFS surface file
* `global_sighdr.fd` - prints header information about a GFS sigma file