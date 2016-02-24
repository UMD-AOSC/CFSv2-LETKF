#!/usr/bin/env python
import argparse
import datetime as dt
import subprocess as sp
import os, shutil, sys
import hashlib

import numpy as np
import netCDF4 as nc


sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python')
import cfs
import obsio
import grdio


## get the command line arguments
parser = argparse.ArgumentParser(description=(
    "Generate synthetic observations for the atmosphere using the "
    "lat/lon/pressure/error information from real observations, but using"
    " the corresponding values from a nature run. To be used for "
    "perfect model experiments. For simplicity obs are generated "
    "at the closest grid points. Only the observations at the analysis "
    "time are used (t.dat)"))
parser.add_argument("nature", metavar="NATURE_PATH", help=(
    "Path to the folder containing the nature run from which observation"
    " values will be generated"))
parser.add_argument("obs", metavar="OBS_PATH", help=(
    "Path to the folder containing the actual observations. These "
    "existing observations will be used only in determining the "
    "synthetic obs locations and errors."))
parser.add_argument("startdate", metavar="START", help=(
    "Start date in YYYYMMDDHH format"))
parser.add_argument("enddate", metavar="END", help=(
    "End date in YYYYMMDDHH format"))
parser.add_argument("output", metavar="OUTPUT_PATH", help=(
    "Directory to place the created synthetic observations in"))
parser.add_argument("--force", "-f", action="store_true")

## parse the arguments
print "Synthetic atmosphere observation generation script.\n"
args = parser.parse_args()
args.startdate = dt.datetime.strptime(args.startdate, "%Y%m%d%H")
args.enddate = dt.datetime.strptime(args.enddate, "%Y%m%d%H")
args.tmpdir = os.getenv("TMP_DIR_LOCAL")+'/genobs_atm_'+hashlib.md5(args.output).hexdigest()[:6]
print "parameters: "+str(args)
cdate = args.startdate


## ------------------------------------------------------------


## setup a temporary work directory
if os.path.exists(args.tmpdir):
    print '[WARNING] temp directory already exists, deleting. '+args.tmpdir
    if (not args.force):
        if (raw_input('Are you sure you want to do this (y/n):') != 'y'):
            sys.exit(1)
    shutil.rmtree(args.tmpdir)
os.makedirs(args.tmpdir)
for f in ['ss2grd','grdctl']:
    shutil.copy(os.getenv("CFS_LETKF_ROOT")+'/util/bin/'+f, args.tmpdir+'/'+f)

    
##  read the spectral resolution from one of the nature files
natfile = os.path.abspath(args.nature+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d%H/%Y%m%d%H.sig"))
ares = int(sp.check_output("global_sighdr {} jcap".format(natfile), shell=True))
print "  Nature files are T{}".format(ares)


## write out the ss2grd namelist.
ares_x,ares_y = cfs.getAtmRes(ares)
print "  Nature files are a {} x {} grid".format(ares_x,ares_y)
with open(args.tmpdir+'/ssio.nml', 'w') as f:
    f.write("$common_gfs nlon={} nlat={} gfs_jcap={} nlev=64 /".format(
        ares_x,ares_y,ares))

    
## create the grd ctl file
sp.call('grdctl 0 0 0 0 0 x > grd.ctl', shell=True, cwd=args.tmpdir)
grdctl = grdio.GradsCtl(args.tmpdir+'/grd.ctl')


## obsid mappings
obsidmap = {
    1100 : 'ps',
    1210 : 't',
    1220 : 'q',
    1250 : 'u',
    1251 : 'v'}


## function to find the closes grid point dimension
##  caches the results in the hash table to speed things up
dim_cache = {}
def nearestDim(coord, dims, *dimV):
    ## if this dimension has a label (e.g. "x" or "y")
    ## check to see if we have seen this coordinate before,
    ## if so, just return its value we computed already
    if len(dimV) > 0:
        if dimV[0] not in dim_cache:
            dim_cache[dimV[0]] = {}
        if coord in dim_cache[dimV[0]]:
            return dim_cache[dimV[0]][coord]

    ## start searching for the closes lon/lat/height that matches this
    ## given coordinate
    nrIdx = 0
    nrDif = 1e20
    for d in range(len(dims)):
        v = abs(coord-dims[d])
        if v < nrDif:
            nrIdx = d
            nrDif = v

    ## cache the results in hash table for faster lookup next time
    if len(dimV) > 0:
        dim_cache[dimV[0]][coord] = nrIdx
    return nrIdx



## ------------------------------------------------------------
## create the observations

count = 0
while cdate <= args.enddate:
    print ""
    print "Processing " + cdate.strftime("%Y%m%d%H")
    badObs = {}  ## a list of unknown observation ids read in
                 ## and the number of those observations seen
    goodObs = {} ## a list of the observation ids that were
                 ## used and the number of those generated
    usedLoc = set([]) ## a hash of x/y/z/obid values so that
                      ## we only generate 1 ob per grid point    
    for x in obsidmap:
        goodObs[x] = 0

        
    ## read in the observation locations
    print "  reading obs file..."
    obsfile = os.path.abspath(args.obs+cdate.strftime(
        "/%Y/%Y%m/%Y%m%d/%Y%m%d%H/t.dat"))
    obs = obsio.read(obsfile)
    print "  {} observations loaded.".format(len(obs))

    
    ## generate an x/y grid coordinate for each observation
    print "  calculating grid x/y coords..."
    obs_xy = []
    for o in obs:
        x = nearestDim(o[1], grdctl.x, 'x')
        y = nearestDim(o[2], grdctl.y, 'y')
        obs_xy.append((x,y))    

        
    ## convert the nature run data into gridded format
    ## and then load it
    print "  Processing nature run..."
    sigfile = os.path.abspath(args.nature+cdate.strftime(
        "/%Y/%Y%m/%Y%m%d/%Y%m%d%H/%Y%m%d%H.sig"))
    sfcfile = os.path.abspath(args.nature+cdate.strftime(
        "/%Y/%Y%m/%Y%m%d/%Y%m%d%H/%Y%m%d%H.sfc"))
    os.symlink(sigfile, args.tmpdir+'/fort.11')
    os.symlink(sfcfile, args.tmpdir+'/fort.12')    
    sp.call('ss2grd',shell=True, cwd=args.tmpdir)
    print "  Loading nature run..."
    nat = grdio.readGrd(grdctl, args.tmpdir+'/fort.31')

    
    ## determine the surface pressure for each x/y, if given by a ps
    ## observation, this makes determining the z value more exact
    xy_ps = {}
    for i in range(len(obs)):
        o = obs[i]
        if o[0] == 1100: #Ps obs
            xy_ps[obs_xy[i]] = o[4]
        o = None

        
    ## for each observation, determine the z level
    obs_z = []
    for i in range(len(obs)):
        if obs[i][0] == 1100:
            ## PS always is at the surface... duh
            z = 0
        else:
            x,y = obs_xy[i]            
            if obs_xy[i] in xy_ps:
                ## we prevopusly found an observation of ps at this x/y
                ##   we can calculate P of the nature it should occur at.
                ## important to do for raobs
                ps = xy_ps[obs_xy[i]]
                o_p = obs[i][3]/ps*(nat['ps'][0,y,x]/100)
            else:
                ## otherwise, just use the P that is given by the real ob
                o_p = obs[i][3]
            ## find the z value based on the P
            vp = nat['p'][:,y,x]/100
            z = nearestDim(o_p, vp)
        obs_z.append(z)

    
    ## generate values for each observation
    synth_obs = []
    for i in range(len(obs)):
        o = obs[i]

        ## make sure this is an observation ID we recognize
        if not o[0] in obsidmap:
            if not o[0] in badObs:
                badObs[o[0]] = 0
            badObs[o[0]] += 1
            continue
        
        ## calculate the x/y/z coords
        x,y = obs_xy[i]
        z = obs_z[i]
        
        ## generate a value
        err = o[5]
        val = nat[obsidmap[o[0]]][z,y,x]
        if o[0] == 1100:
            val /= 100
        val += np.random.normal(0,err)
        
        if o[0] == 1100:  ## Ps            
            hgt = nat['orog'][0,y,x]
        else:
            hgt = nat['p'][z,y,x]/100
            
        o2 = (o[0], grdctl.x[x], grdctl.y[y], hgt, val, err, o[6])

        
        ## obs thinning, make sure only 1 observation of any given type
        ##  is generated at a given x/y/z
        hkey = (o[0], x, y, z)
        if hkey in usedLoc:
            continue
        usedLoc.add(hkey)

        ## add the observation
        goodObs[o[0]] += 1
        synth_obs.append(o2)

        
    ## diagnostic summary info
    for o in badObs:
        print "  [Error] Unkown obsid ({}) found {} times".format(o, badObs[o])
    total = 0
    for o in goodObs:
        print "  {:7d} {:2} obs".format(goodObs[o],obsidmap[o])
        total += goodObs[o]
    print "  {:7d} TOTAL obs".format(total)

    
    ## write the observatiosn out
    filename = args.output+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d%H.dat")
    obsio.write(synth_obs,filename)

    
    ## cleanup
    sp.call("rm fort.*",shell=True, cwd=args.tmpdir)

    
    ## increment the date
    cdate += dt.timedelta(hours=6)

## final cleanup
shutil.rmtree(args.tmpdir)
