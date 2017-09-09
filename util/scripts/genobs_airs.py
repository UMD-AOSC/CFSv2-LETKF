#!/usr/bin/env python

## built-in modules
import argparse
import datetime as dt
import subprocess as sp
import os, shutil, sys
import hashlib

## 3rd party modules
import numpy as np
import netCDF4 as nc

## my modules
sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python')
import cfs
import obsio
import grdio



################################################################################
################################################################################



## get the command line arguments
parser = argparse.ArgumentParser(description=(""))
parser.add_argument("nature", metavar="NATURE_PATH", help=(
    "Path to the folder containing the nature run from which observation"
    " values will be generated"))
parser.add_argument("startdate", metavar="START", help=(
    "Start date in YYYYMMDDHH format"))
parser.add_argument("enddate", metavar="END", help=(
    "End date in YYYYMMDDHH format"))
parser.add_argument("output", metavar="OUTPUT_PATH", help=(
    "Directory to place the created synthetic observations in"))
parser.add_argument("--force", "-f", action="store_true", help=(
    "Forces the removing of the destination directory and temp "
    "directory if they already exist"))

args = parser.parse_args()
args.is3d = True #TODO: allow this to be configured
args.startdate = dt.datetime.strptime(args.startdate, "%Y%m%d%H")
args.enddate = dt.datetime.strptime(args.enddate, "%Y%m%d%H")      
        
            
## generate a temporary directory that has a random hash at the end,
## so that the user can run more than one copy of this script at a time
## for different experiments
args.tmpdir = os.getenv("TMP_DIR_LOCAL")+'/genobs_airs_'+hashlib.md5(args.output).hexdigest()[:6]
print "parameters: "+str(args)
cdate = args.startdate



################################################################################
# get things ready
################################################################################


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
ares_x, ares_y = cfs.getAtmRes(ares)
print "  Nature files are T{}".format(ares)
print "  Nature files are a {} x {} grid".format(ares_x, ares_y)


## write out the ss2grd namelist, in preparation of running ss2grd (to convert
## the spectal nature run file into grided that I can read in).
with open(args.tmpdir+'/ssio.nml', 'w') as f:
    f.write("$common_gfs nlon={} nlat={} gfs_jcap={} nlev=64 /".format(
        ares_x,ares_y,ares))

    
## create the grd ctl file, for loading in grdded data later
sp.call('grdctl 0 0 0 0 0 x > grd.ctl', shell=True, cwd=args.tmpdir)
grdctl = grdio.GradsCtl(args.tmpdir+'/grd.ctl')


## obsid mappings,
## allows lookup of the variable name stored in the grd file
obsidmap = {
    1100 : 'ps',
    1210 : 't',
    1220 : 'q',
    1250 : 'u',
    1251 : 'v'}

## platformid mapping,
## allows lookup of observation platform name given the id
platformidmap = {
     1 : "ADPUPA",     2 : "AIRCAR",     3 : "AIRCFT",
     4 : "SATWND",     5 : "PROFLR",     6 : "VADWND",
     7 : "SATEMP",     8 : "ADPSFC",     9 : "SFCSHP",
    10 : "SFCBOG",    11 : "SPSSMI",    12 : "SYNDAT",
    13 : "ERS1DA",    14 : "GOESND",    15 : "QKSWND",
    16 : "MSONET",    17 : "SPGIPW",    18 : "RASSDA",
    19 : "WDSATR",    20 : "ASCATW",    21 : "TMPAPR",}


## read in the sigma level definitions
atmlvls = cfs.getAtmLevels()

## ------------------------------------------------------------
## function to find the closest grid point dimension.
##  caches the results in the hash table to speed things up
## ------------------------------------------------------------
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

    ## start searching for the closest lon/lat/height that matches this
    ## given coordinate
    nrIdx = 0       ## index of gridpoint dimension closest so far
    nrDif = 1e20    ## distance of closes point
    for d in range(len(dims)):
        v = abs(coord-dims[d])
        if v < nrDif:
            nrIdx = d
            nrDif = v

    ## cache the results in hash table for faster lookup next time
    if len(dimV) > 0:
        dim_cache[dimV[0]][coord] = nrIdx
        
    return nrIdx



################################################################################
## create the observations
################################################################################
basedate = dt.datetime(1993,1,1)
airs_lvls = (1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50)
q_err={
    1000:4,
    925:4,
    850:3,
    700:2,
    600:1,
    500:0.5,
    400:0.2,
    300:0.1}
count = 0
## for each 6 hour date in the list given
while cdate <= args.enddate:  
    print ""
    print "Processing " + cdate.strftime("%Y%m%d%H")
    
#    badObs = {}  ## a list of unknown observation ids read in
#                 ## and the number of those observations seen
    platforms={} ## a list of the platform types seen
                ## and the nubmer of those observations seen
    goodObs = dict( [ (x,0) for x in obsidmap ] )
                ## a list of the observation types that were
                ## used and the number of those generated
    usedLoc = set([]) ## a hash of x/y values so that
                     ## we only generate 1 ob per grid point

    ## determine what the day in track cycle this is
    trackFile = 'genobs_airs.dat/{:03d}_{:02d}.npy'.format(((cdate-basedate).days % 16)+1, cdate.hour)
    print "  using ",trackFile
    tracks = np.load(trackFile)

       
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

    
     ## generate values for each observation
    synth_obs = []
    print '  generating obs values...'
    for pt in tracks:
        ## get x/y point
        x = nearestDim(pt[1], grdctl.x, 'x')
        y = nearestDim(pt[0], grdctl.y, 'y')
        hkey = (x,y)
        if hkey in usedLoc:
            continue
        usedLoc.add(hkey)

        ##for each z point
        usedZLoc = set([])         
        for lvl in airs_lvls:
            z = nearestDim(lvl, nat['p'][:,y,x]/100)            
            hkey = (z,)
            if hkey in usedZLoc:
                continue
            usedZLoc.add(hkey)

            ## generate values
            err = 1
            val = nat['t'][z,y,x]            
            err_add = 1e10
            while abs(err_add) > err*3: # dont have really big errors...
                err_add = np.random.normal(0,err)
            val += err_add ## add gaussian noise
            hgt = nat['p'][z,y,x]/100
            newob = (1210, grdctl.x[x], grdctl.y[y], hgt, val, err, 14)

            ## add the observation
            goodObs[newob[0]] += 1
            if not newob[6] in platforms:
                platforms[newob[6]] = 1
            else:
                platforms[newob[6]] = platforms[newob[6]] + 1           
            synth_obs.append(newob)
            

            ## generate values
            if lvl >=  300:
                val = nat['q'][z,y,x]
                if val <= 0:
                    continue
                err = val/10.0
                err_add = 1e10
                while abs(err_add) > err*3: # dont have really big errors...
                    err_add = np.random.normal(0,err)
                val += err_add ## add gaussian noise
                if val <= 0:
                    continue
                newob = (1220, grdctl.x[x], grdctl.y[y], hgt, val, err, 14)
            
                ## add the observation
                goodObs[newob[0]] += 1
                if not newob[6] in platforms:
                    platforms[newob[6]] = 1
                else:
                    platforms[newob[6]] = platforms[newob[6]] + 1           
                synth_obs.append(newob)


        
    ## diagnostic summary info
    print "  ====== Observations variables ========"
    total = 0    
    for o in goodObs:
        print "  {:7d} {:2} obs".format(goodObs[o],obsidmap[o])
        total += goodObs[o]
    print "  {:7d} TOTAL obs".format(total)

    print ""
    print "  ====== Observations Platforms ========="
    for o in platforms:
        print "  {:7d} {:2}({:02d}) obs".format(platforms[o], platformidmap[int(o)], int(o))

    
    ## write the observatiosn out
    filename = args.output+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d%H/t.dat")
    obsio.write(synth_obs,filename)

    
    ## cleanup
    sp.call("rm fort.*",shell=True, cwd=args.tmpdir)

    
    ## increment the date
    cdate += dt.timedelta(hours=6)

    
## final cleanup
shutil.rmtree(args.tmpdir)
