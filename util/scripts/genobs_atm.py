#!/usr/bin/env python

################################################################################
## genobs_atm.py
## CFSv2-LETKF
##
## Generate synthetic observations for the atmosphere using the
## lat/lon/pressure/error information from real observations, but using
## the corresponding values from a nature run. To be used for perfect model
## experiments. For simplicity obs are generated at the closest grid points.
##
## Travis Sluka, 2016
## University of Maryland
################################################################################

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
parser = argparse.ArgumentParser(description=(
    "Generate synthetic observations for the atmosphere using the "
    "lat/lon/pressure/error information from real observations, but using"
    " the corresponding values from a nature run. To be used for "
    "perfect model experiments. For simplicity obs are generated "
    "at the closest grid points."))
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
parser.add_argument("--force", "-f", action="store_true", help=(
    "Forces the removing of the destination directory and temp "
    "directory if they already exist"))
parser.add_argument("--platforms","-p", help=(
    "a comma separated (no spaces!) list of observations platform "
    "IDs to use. See the wiki for further documentation. "
    "E.g., -p 1,2,3,8,9,19 is appropriate for perfect model experiments ("
    " most PREPBUFR obs except SATWND)"))


args = parser.parse_args()
args.startdate = dt.datetime.strptime(args.startdate, "%Y%m%d%H")
args.enddate = dt.datetime.strptime(args.enddate, "%Y%m%d%H")
if args.platforms is not None:
    args.platforms = [int(x) for x in args.platforms.split(',')]
## generate a temporary directory that has a random hash at the end,
## so that the user can run more than one copy of this script at a time
## for different experiments
args.tmpdir = os.getenv("TMP_DIR_LOCAL")+'/genobs_atm_'+hashlib.md5(args.output).hexdigest()[:6]
print "parameters: "+str(args)
cdate = args.startdate



################################################################################
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


## write out the ss2grd namelist.
with open(args.tmpdir+'/ssio.nml', 'w') as f:
    f.write("$common_gfs nlon={} nlat={} gfs_jcap={} nlev=64 /".format(
        ares_x,ares_y,ares))

    
## create the grd ctl file
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



## ------------------------------------------------------------
## create the observations

count = 0
## for each 6 hour date in the list given
while cdate <= args.enddate:
    
    print ""
    print "Processing " + cdate.strftime("%Y%m%d%H")
    
    badObs = {}  ## a list of unknown observation ids read in
                 ## and the number of those observations seen
    platforms={} ## a list of the platform types seen
                 ## and the nubmer of those observations seen
    goodObs = dict( [ (x,0) for x in obsidmap ] )
                 ## a list of the observation types that were
                 ## used and the number of those generated
    usedLoc = set([]) ## a hash of x/y/z/obid values so that
                      ## we only generate 1 ob per grid point

    ## obs, obs_xy, and obs_z will be generated below, all arrays
    ##  are the same length and corresponding elements go together
    ##  between the 3.
        
    ## read in the observation locations
    print "  reading obs file..."
    obs = []
    ## this will read in all timeslots for the 6 hour observation window
    for f in ["-3","-2","-1","","+1","+2","+3"]:
        obsfile = os.path.abspath(
            args.obs+cdate.strftime("/%Y/%Y%m/%Y%m%d/%Y%m%d%H/t")+f+".dat")
        obs += obsio.read(obsfile)
    print "  {} observations loaded.".format(len(obs))

    ## remove unwanted platforms if the user explicitly defined
    ## the list of platforms to use
    if args.platforms:
        obs = filter(lambda x: x[6] in args.platforms, obs)

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

        
    ## for each observation, determine the z level
    obs_z = []
    for i in range(len(obs)):
        if obs[i][0] == 1100:
            z = 0   ## PS always is at the surface... duh
        elif obs[i][6] in (8,9,19):
            ## ADPSFC, SFCSHP, and WDSATR are at the surface
            z = 0
        else:
            ## o_p will be the pressure level at which to generate the observation            
            x, y = obs_xy[i]
            if obs_xy[i] in xy_ps:
                ## we prevoiously found an observation of ps at this x/y
                ## we can calculate P of the nature it should occur at.
                ## important to do for raobs
                ps = xy_ps[obs_xy[i]]
                o_p = obs[i][3]/ps*(nat['ps'][0,y,x]/100)
            else:
                ## otherwise, just use the P that is given by the real ob
                o_p = obs[i][3]
                
            ## find the z value based on the P
            z = nearestDim(o_p, nat['p'][:,y,x]/100)
        obs_z.append(z)

    
    ## generate values for each observation
    synth_obs = []
    for i in range(len(obs)):
        o = obs[i]

        ## make sure this is an observation ID we recognize
        if o[0] not in obsidmap:
            if o[0] not in badObs:
                badObs[o[0]] = 0
            badObs[o[0]] += 1
            continue
       
        ## calculate the x/y/z coords
        x, y = obs_xy[i]
        z = obs_z[i]
        
        ## obs thinning, make sure only 1 observation of any given type
        ##  is generated at a given x/y/z
        hkey = (o[0], x, y, z)
        if hkey in usedLoc:
            continue
        usedLoc.add(hkey)

        
        ## generate a value
        err = o[5]
        val = nat[obsidmap[o[0]]][z,y,x]
        if o[0] == 1100: 
            val /= 100  ## convert pressure from Pa to hPa
        val += np.random.normal(0,err) ## add gaussian noise
        
        if o[0] == 1100:
            ## PS obs give height as the actual terrain height
            hgt = nat['orog'][0,y,x]
        else:
            ## other obs give height as pressure in mb
            hgt = nat['p'][z,y,x]/100
            
        newob = (o[0], grdctl.x[x], grdctl.y[y], hgt, val, err, o[6])

        
        ## add the observation
        goodObs[newob[0]] += 1
        if not newob[6] in platforms:
            platforms[newob[6]] = 1
        else:
            platforms[newob[6]] = platforms[newob[6]] + 1
            
        synth_obs.append(newob)

        
    ## diagnostic summary info
    for o in badObs:
        print "  [Error] Unkown obsid ({}) found {} times".format(o, badObs[o])
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
