#!/usr/bin/env python
import argparse
import datetime as dt
import subprocess as sp
import os, shutil, sys
from glob import glob


import numpy as np
import netCDF4 as nc
import scipy.spatial

sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python')
import obsio



## get the command line arguments
parser = argparse.ArgumentParser(
    description=("Generate synthetic observations for the ocean using the lat/lon/depth/error "
                 "information from real obs, but the corresponding values from a nature run. "
                 "To be used for perfect model experiments. This script assumes the existing "
                 "real observations are in daily binned files, and 6hr binned files "
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

args = parser.parse_args()
args.startdate = dt.datetime.strptime(args.startdate, "%Y%m%d%H")
args.enddate = dt.datetime.strptime(args.enddate, "%Y%m%d%H")
print "Parameters: "+str(args)

## 
cdate = args.startdate


## Read in the grid definition file, and save lat/lon pairs in a KD tree
## for easy index lookup,
print "Reading grid definition and proccessing lat/lons..."
ncdat = nc.Dataset('../../support/fix/fix_om/grid_spec_05.nc.T62')
lons = ncdat.variables["x_T"][:]
lons_min = np.min(lons)
lons_max = np.max(lons)
lats = ncdat.variables["y_T"][:]
depths = ncdat.variables["zt"][:].tolist()
xs,ys = np.meshgrid(np.arange(0,lons.shape[1]), np.arange(0,lons.shape[0]))
xs = xs.reshape(-1)
ys = ys.reshape(-1)
zs = np.arange(0,len(depths))
ll_tree = scipy.spatial.KDTree(zip(lons.reshape(-1),lats.reshape(-1)))
ncdat.close()

def getobs(ncdat, t_start, t_end, nightonly=True):
    obs = []
    times = ncdat.variables['sst_dtime'][0,0]
    olats = ncdat.variables['lat'][:]
    olons = ncdat.variables['lon'][:]
    err   = ncdat.variables['SSES_standard_deviation_error'][0] 
    rej_flag = ncdat.variables['rejection_flag'][0]

    t1 = dt.datetime.strptime(ncdat.start_date+" "+ncdat.start_time, "%Y-%m-%d %H:%M:%S UTC")
    t2 = dt.datetime.strptime(ncdat.stop_date+" "+ncdat.stop_time, "%Y-%m-%d %H:%M:%S UTC")    
    
    s1 = max((t_start-t1).total_seconds(),0)
    s2 = min((t_end-t1).total_seconds(), (t2-t1).total_seconds())
    usedLoc = set([])
    
    count = 0
    for i in range(1,len(times),10):
        if not (times[i] >= s1 and times[i] <= s2):            
            continue
        for j in range(1,olats.shape[0],10):
            if (nightonly and olats[j,i-1] < olats[j,i]):
                continue ## ignore the ascending leg 
            if rej_flag[j,i] == 0:
                l = olons[j,i]
                if l < lons_min:
                    l += 360
                if l > lons_max:
                    l -= 360
                
                res = ll_tree.query((l,olats[j,i]))
                y = ys[res[1]]
                x = xs[res[1]]
                ll = (lons[y,x], lats[y,x])
                if ll in usedLoc:
                   continue
                usedLoc.add(ll)
                obs.append((x,y, err[j,i]))
        count += 1
    return obs


## produce the observations
while cdate <= args.enddate:
    obs2 = []
    print ""
    print "Processing "+cdate.strftime("%Y%m%d%H")

    ## open the nature run
    dat_filename = args.nature+cdate.strftime(
        '/%Y/%Y%m/%Y%m%d/%Y%m%d%H/%Y%m%d%H.ocean_temp_salt.res.nc')
    dat = nc.Dataset(dat_filename).variables
    
    ## get the files available for this day (which includes the end of the previous day and
    ## beginning of next day, since 6 hour windows are used)
    c_start = cdate - dt.timedelta(hours=3)
    c_end   = cdate + dt.timedelta(hours=3)
    print "generating obs between ",c_start," and ", c_end
    infiles = glob(args.obs+c_start.strftime("/%Y/%Y%m/%Y%m%d/*"))
    if not c_start.day == c_end.day:
        infiles += glob(args.obs+c_end.strftime("/%Y/%Y%m/%Y%m%d/*"))

    ## start searching the files, processing any tracks that fall within the time we are looking at
    for f in sorted(infiles):
        ncdat = nc.Dataset(f)
        f_start = dt.datetime.strptime(ncdat.start_date+" "+ncdat.start_time, "%Y-%m-%d %H:%M:%S UTC")
        f_end   = dt.datetime.strptime(ncdat.stop_date +" "+ncdat.stop_time,  "%Y-%m-%d %H:%M:%S UTC")
        if (f_start < c_end and f_end > c_start):
            print "  using", f, f_start, f_end
            obsLoc = getobs(ncdat, c_start, c_end)
            for o in obsLoc:
                val = dat['temp'][0,0,o[1],o[0]]
                val += np.random.normal(0, o[2])
                obs2.append((2110, lons[o[1],o[0]], lats[o[1],o[0]], 0, val, o[2], 0))
        ncdat.close()
                                    

    # ## write out the observations
    print "  writing {} synthetic obs".format(len(obs2))
    filename = args.output+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d%H.dat")
    obsio.write(obs2,filename)
    cdate += dt.timedelta(hours=6)
