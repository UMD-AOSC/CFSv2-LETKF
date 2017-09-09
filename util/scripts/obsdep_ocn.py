#!/usr/bin/env python

import os,sys,shutil
import argparse
import subprocess as sp
from glob import glob
import hashlib
import datetime as dt

parser=argparse.ArgumentParser(description=(
    "Applies the observation operator stand-alone to the ocean, either "
    "the ensemble mean of the background or the analysis"))
parser.add_argument('path', help=(
    "Path to the experiment for which the observation operator should "
    "be applied"))
parser.add_argument('obsdir', help=(
    "directory to the observations that should be used for calculating"
    " location and obs departure"))
parser.add_argument('--start', help=(
    "start date, in YYYYMMDD format. If no date is given, all dates "
    "in the path are used"))
parser.add_argument('--end', help=(
    "end date, in YYYYMMDD format. If no date is given, only the --start "
    "date is used"))
parser.add_argument('-a',action="store_true", help=(
    "obs operator should be applied to the analysis. Either -a or -b "
    "needs to be specified"))
parser.add_argument('-b',action="store_true", help=(
    "obs operator should be applied to the background. Either -a or -b "
    "needs to specified"))
parser.add_argument('--ores', required=False, default='05', choices=['05','1x1'], help=(
    "ocean resolution, default: 05"))
parser.add_argument('--force', default=False, action="store_true", help=(
    "if set, forces the recalculation of any already existing obs "
    "departures, otherwise they will be skipped"))
args=parser.parse_args()


args.path = os.path.abspath(args.path)
args.tmpdir   = os.path.abspath(os.getenv("TMP_DIR_LOCAL")+'/ocn_omx_'+hashlib.md5(
    args.path+('a' if args.a else 'b')).hexdigest()[0:6])
args.letkfDir = os.path.abspath(os.getenv("CFS_LETKF_ROOT") + '/letkf-mom')
args.fixDir = os.path.abspath(os.getenv("FIX_DIR_OM"))
args.obsdir = os.path.abspath(args.obsdir)

if args.start is None and args.end is None:
    args.start = "00000000"
    args.end = "99999999"
elif args.start is None:
    print "incorrect specification of start/end dates"
    sys.exit(1)
elif args.end is None:
    args.end = args.start

if (not args.a and not args.b) or (args.a and args.b):
    print "must specify either -a or -b"
    sys.exit(1)

print args


    
# get list of files to process
bgfiles = sorted(glob(args.path+"/{}/mean/????/??????/*/*/*_ocn.grd".format(
    'gues' if args.b else 'anal')))

# filter out only those that have accompanying observations, and correct date
def filterfunc(d):
    d2 = d.split('/')[-1][:10]
    if d2 < args.start or d2 > args.end:
        return False
    return  d2[-2:] == '12'
bgfiles = sorted(filter(filterfunc, bgfiles))

# setup working directory
if os.path.exists(args.tmpdir):
    shutil.rmtree(args.tmpdir)
os.makedirs(args.tmpdir)
os.symlink(args.letkfDir+'/obs/obsop', args.tmpdir+'/obsop')
os.symlink(args.fixDir + '/grid_spec_{}.nc.T126'.format(args.ores),
           args.tmpdir+'/grid_spec.nc')
shutil.copy(args.path+"/cfg/letkf.nml", args.tmpdir+"/input.nml")

# grab a random ocean netcdf file to use for the grd2nc conversion
ncin = os.path.abspath(glob(args.path+'/gues/001/*ocean_temp_salt.res.nc')[0])[:-16]
sp.check_call('cp {}* {}'.format(ncin,args.tmpdir), shell=True)
ncin = os.path.abspath(glob(args.tmpdir+'/*ocean_temp_salt.res.nc')[0])[:-16]

# for each background state
for bg in bgfiles:

    # check to see if this date has already been done    
    cdate = bg.split('/')[-1][:10]
    outfile =  args.path+'/{}/{}/{}/{}/{}/{}_ocn.dat'.format(
        'gues' if args.b else 'anal',
        'omb' if args.b else 'oma',
        cdate[:4], cdate[:6], cdate[:8], cdate)

    print ""
    print "************************************************************"
    if os.path.exists(outfile) and not args.force:
        print cdate+" already exists, skipping"
        continue
    print "Calculating "+outfile

    
    # convert to netcdf, ugh
    script = os.getenv("CFS_LETKF_ROOT")+'/util/scripts/ocn_grd2nc.py'
    sp.check_call(script+' {} {} {} --ores {}'.format(
         bg, ncin, args.tmpdir+'/ocndat.ocean_', args.ores),
                  shell=True, cwd=args.tmpdir)

    # run the observation operator on it
    obsin = args.obsdir+ '/' + '/'.join(bg.split('/')[-5:-1])[:-2]+'.dat'    
    sp.check_call('obsop -obsin {} -gues ocndat -obsout out.dat'.format(obsin),
                  shell=True, cwd=args.tmpdir)

    # move the file
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)      
    shutil.move(args.tmpdir+'/out.dat', outfile)
    sp.check_call('rm ocndat.*', shell=True, cwd=args.tmpdir)

shutil.rmtree(args.tmpdir)
