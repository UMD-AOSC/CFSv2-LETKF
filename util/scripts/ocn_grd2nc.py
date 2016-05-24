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


gradCtlFilename = os.getenv("CFS_LETKF_ROOT")+'/letkf-mom/ocn_1x1.ctl'

############################################################
## Get the command line arguemnts
parser = argparse.ArgumentParser(description=(
    ""))
parser.add_argument('grdfile', metavar='GRD_IN', help=(
    ""))
parser.add_argument('infile_pfx', metavar='NC_IN', help=(
    ""))
parser.add_argument('outfile_pfx', metavar='NC_OUT', help=(
    ""))
args = parser.parse_args()
print args

## load the grd file
print "Loading grads ctl file: "+os.path.abspath(gradCtlFilename)
grdCtl = grdio.GradsCtl(gradCtlFilename)
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
