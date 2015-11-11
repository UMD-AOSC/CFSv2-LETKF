# CFSV2-LETKF
Strongly Coupled LETKF applied to the Climate Forecast System version 2 (CFSv2). Default settings use the T62/L64 version of GFS and the 0.5 degree MOM. Fix files must be obtained separately. Information about how to run most of the scripts can be shown by running with the `-h` flag, (e.g. `./util/get_cfsr -h`

All further documentation, including how to compile and run, are contained on the wiki: [GitHub wiki for CFSv2-LETKF](https://github.com/UMD-AOSC/CFSv2-LETKF/wiki)

####Directory Structure
* `cfs` - unmodified NCEP code for the CFSv2
* `config` - site specific configuration for compiling and running
* `letkf-gfs` - GFS-LETKF, atmospheric data assimilation
* `letkf-mom` - MOM-LETKF, oceanic data assimilation
* `lib-ncep` - support libraries from NCEP required for compiling other executables
* `run` - scripts to initialize and run the CFS / data assimilation cycle
* `util` - other utility programs required by the CFS and/or GFS-LETKF code