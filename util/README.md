This folder contains utility programs required by the CFS. 


####CFSR files
* `get_cfsr` - Downloads CFSR initial conditions from the internet, run `get_cfsr -h` to get details on how to run the script.
* `chgres_cfsr` - Changes spectral resolution of CFSR atmospheric files. Run `chgres_cfsr -h` to get details on how to run the script.

####Grid / scpectral conversion
The following programs are created by running `make`
* `bin/grd2ss` - converts grided data into the spectral format required by the GFS
* `bin/grdctl` - creates an appropriate grads control file for gridded data
* `bin/ss2grd` - converts GFS spectral data to gridded format on model levels
* `bin/ss2grdp` - converts GFS spectral data to gridded format on pressure levels
* `bin/sscycle` - converts a forecast to an analysis by modifying the headers and replacing certain boundary field (e.g. ozone)

