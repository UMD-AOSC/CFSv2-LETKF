#!/usr/bin/env python
import argparse
import datetime as dt
import subprocess as sp
import os, shutil, sys

import numpy as np
import netCDF4 as nc
import scipy.spatial

sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python')
import obsio



## get the command line arguments
parser = argparse.ArgumentParser(
    description=("Generate synthetic observations for the ocean using the lat/lon/depth/error information from "
                 "real obs, but the corresponding values from a nature run. To be used for perfect model experiments. "
                 "This script assumes the existing real observations are in daily binned files, and 6hr binned files "
                 "will be produced. For simplicity obs are generated at the closest grid points."))
parser.add_argument("nature", metavar="NATURE_PATH", help=(
    "path to the folder containing the nature run from which "
    "observations values will be generated"))
parser.add_argument("obs", metavar="REAL_OBS_PATH", help=(
    "path to the folder containing the actual observations. These "
    "existing observations will be used only in determining the "
    "synthetic obs locations and errors."))
parser.add_argument("startdate", metavar="START", help=(
    "Start date in YYYYMMDDHH format"))
parser.add_argument("enddate", metavar="END", help=(
    "End date in YYYYMMDDHH format"))
parser.add_argument("output", metavar="OUTPUT_PATH", help=(
    "Directory to place the created synthetic observations in"))
parser.add_argument("--split", action="store_true", default=False, help=(
    "by default, obs are left in 24 hour bins, this option randomly splits obs into 6 hour bins in the day"))

args = parser.parse_args()
args.startdate = dt.datetime.strptime(args.startdate, "%Y%m%d%H")
args.enddate = dt.datetime.strptime(args.enddate, "%Y%m%d%H")
print "Parameters: "+str(args)

## 
cdate = args.startdate
cdate = cdate.replace(hour=12)



## Read in the grid definition file, and save lat/lon pairs in a KD tree
## for easy index lookup,
print "Reading grid definition and proccessing lat/lons..."
ncdat = nc.Dataset('../../support/fix/fix_om/grid_spec_05.nc.T62')
lons = ncdat.variables["x_T"][:]
lats = ncdat.variables["y_T"][:]
depths = ncdat.variables["zt"][:].tolist()
xs,ys = np.meshgrid(np.arange(0,lons.shape[1]), np.arange(0,lons.shape[0]))
xs = xs.reshape(-1)
ys = ys.reshape(-1)
zs = np.arange(0,len(depths))
ll_tree = scipy.spatial.KDTree(zip(lons.reshape(-1),lats.reshape(-1)))
ncdat.close()


## produce the observations
ll_hr = None
ll_hr_date = None
obs = None
while cdate <= args.enddate:
    print ""
    print "Processing "+cdate.strftime("%Y%m%d%H")

    ## read in observations if they haven't been yet for this day
    ##  this block will create:
    ##  ll_hr (an array of lat lons used in each time)
    ##  ll_hr_date (the day that was loaded in)
    ##  obs (an array of all obs for the day
#    if ll_hr is None or ll_hr_date != cdate.strftime("%Y%m%d"):
#        ll_hr_date = cdate.strftime("%Y%m%d")
    obsfile = os.path.abspath(args.obs+cdate.strftime(
           "/%Y/%Y%m/%Y%m%d/%Y%m%d.dat"))
    print "  reading obsfile: "+obsfile
    obs = obsio.read(obsfile)
    ##   find the list of unique lat/lons
    ll = set([])
    for o in obs:
        ll.add((o[1],o[2]))        
    print "  read {} real obs with {} unique lat/lon".format(len(obs),len(ll))

                
    ## open the nature run
    dat_filename = args.nature+cdate.strftime(
        '/%Y/%Y%m/%Y%m%d/%Y%m%d%H/%Y%m%d%H.ocean_temp_salt.res.nc')
    dat = nc.Dataset(dat_filename).variables
    
    ## for this time bucket find the observations we are going to use
    obs2=[]
    prevll = (1e6,1e6)
    prevllres = None
    for o in obs:
        ll = (o[1],o[2])
        if ll == prevll:
            res = prevllres
        else:
            res = ll_tree.query(ll)
            prevllres=res
            prevll=ll
        y = ys[res[1]]
        x = xs[res[1]]
        z = min(range(len(depths)),key=lambda i: abs(depths[i]-o[3]))
        oid = o[0]
        err = o[5]            
        if oid == 2210:
            val = dat['temp'][0,z,y,x]
        elif oid == 2220:
            val = dat['salt'][0,z,y,x]
        else:
            print 'unexpected obs type ',oid
            sys.exit(1)
        val += np.random.normal(0, err)                
        obs2.append(
            (oid, lons[y][x], lats[y][x], depths[z], val, err, 0))

    ## write out the observations
    print "  writing {} synthetic obs".format(len(obs2))
    filename = args.output+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d.dat")
    obsio.write(obs2,filename)
    cdate += dt.timedelta(hours=24)
