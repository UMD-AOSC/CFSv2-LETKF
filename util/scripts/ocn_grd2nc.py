#!/usr/bin/env python

## built-in modules
import argparse
import os,shutil,sys
from glob import glob

## 3rd party modules
import netCDF4 as nc
import numpy as np

## custom modules
sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python')
import grdio


############################################################
## Get the command line arguemnts
parser = argparse.ArgumentParser(description=(
    "Converts the ocean grid from the LETKF's raw format (.grd file) to a netCDF format (.nc)"))
parser.add_argument('grdfile', metavar='GRD_IN', help=(
    "The source .grd file that we will be converting"))
parser.add_argument('infile_pfx', metavar='NC_IN', help=(
    "input netCDF files that are the same format as the ouptut this script produces (These "
    " files are not altered, they are copied and then the grd data is placed into the newly"
    " copied files"))
parser.add_argument('outfile_pfx', metavar='NC_OUT', help=(    
    "file prefix/directory location for the output files"))
parser.add_argument('--ores', required=False, choices=['05','1x1'], default='05', help=(
    "input/output grid resolution"))
args = parser.parse_args()

print "-------------------------------------------"
print "ocn_grd2nc.py"
args.gradCtlFilename = os.getenv("CFS_LETKF_ROOT")+'/letkf-mom/ocn_{}.ctl'.format(args.ores)


print args

## load the grd file
print "Loading grads ctl file: "+os.path.abspath(args.gradCtlFilename)
grdCtl = grdio.GradsCtl(args.gradCtlFilename)
print "Loading file: " + os.path.abspath(args.grdfile)
grdDat = grdio.readGrd(grdCtl, args.grdfile)


## make sure the required input netcdf files are there, and copy them to the destination
ncSfx=('temp_salt','velocity','sbc')
for s in ncSfx:
    ## make sure source file exists
    filePattern=args.infile_pfx+s+".*"
    files = glob(filePattern)
    if len(files) == 0:
        print ("ERROR: file {} not found".format(os.path.abspath(filePattern)))
        sys.exit(1)

    ## make sure destination directory exists
    dstDir = os.path.dirname(os.path.abspath(args.outfile_pfx))
    if not os.path.exists(dstDir):
        print 'making directory '+dstDir
        os.makedirs(dstDir)

    ## copy the nc file over
    shutil.copy(files[0], args.outfile_pfx+s+".res.nc")

## for each netcdf file
for s in ncSfx:
    ## open the netcdf file
    ncFilename = args.outfile_pfx+s+".res.nc"
    print "Editing: "+ncFilename
    ncDat = nc.Dataset(ncFilename,'r+')

    ## overwrite with the grd data
    if s == 'velocity':
        vvars = ( ('u','u') , ('v','v') )
    elif s == 'temp_salt':
        vvars = ( ('t','temp') , ('s','salt') )
    else:
        print "WARNING: sbc file is ignored for now"
        continue
        
    for v in vvars:
        ncDat.variables[v[1]][0] = np.nan_to_num(grdDat[v[0]][:])

    ncDat.close()
