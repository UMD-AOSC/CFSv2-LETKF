#!/usr/bin/env python
import argparse
import sys, os, shutil
sys.path.insert(1,os.getenv('CFS_LETKF_ROOT')+"/common/python")
import obsio
from glob import glob
import scipy.spatial



parser = argparse.ArgumentParser(description=(
    "Thins the atmospheric PREPBUFR observations. Thinning is done separately "
    "for each platform type / variable type combination"))
parser.add_argument('indir', help=(
    "input directory. Files below this directory should have "
    "the format YYYY/YYYYMM/YYYYMMDD/YYYYMMDDHH/t*.dat"))
parser.add_argument('outdir', help=(
    "output directory"))
parser.add_argument('--hz', metavar='deg', type=float, default=2, help=(
    "horizontal binning in degrees, (default: 2)"))
parser.add_argument('--hz_sfc', metavar='deg', type=float, default=3, help=(
    "horizontal binning, for the surface station obs (default: 3)"))
parser.add_argument('--vt', metavar='mb', type=float, default=50, help=(
    "vertical binning, in millibars (default: 50)"))
args = parser.parse_args()
print args

# create output directory
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)


    
# ------------------------------------------------------------
# ------------------------------------------------------------    
def hzthin(obs):
    if obs[0][6] in (8,):
        hzbin = args.hz_sfc ## extra thinning of surface obs
                  ## since we get these every 1 hour
    else:
        hzbin = args.hz #degrees

    if obs[0][6] in (4, 8, 9):
        vtbin = 2000 # only 1 ob for each hz location
           #for surface obs
    else:
        vtbin = args.vt

    obsout = []
    
    idx = set([])
    for o in obs:
        idxbin = ( o[0],
                   round(o[1] / hzbin),
                   round(o[2] / hzbin),
                   round(o[3] / vtbin))
        if idxbin not in idx:
            idx.add(idxbin)
            obsout += [o,]
    
    return obsout
    
    
# get list of obs to thin
infiles = sorted(glob(args.indir+'/*/*/*/*/*.dat'))
for f in infiles:
    # create directory for the output file
    outfile = args.outdir+'/'+'/'.join(f.split('/')[-5:])
    dirname = os.path.dirname(outfile)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    print f

    # load obs
    obsin = obsio.read(f)
    print "Input obs: ",len(obsin)

    #figure out which platform numbers are present
    pid = set( [ o[6] for o in obsin] )
    
    # horizontal thinning on each platform separately
    obsout = []
    for p in pid:
        obs = filter(lambda o: o[6] == p, obsin)
        obs2 = hzthin(obs)
        obsout += obs2
        
    print "Output obs: ", len(obsout)

    obsio.write(obsout, outfile)
    



