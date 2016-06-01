#!/usr/bin/env python

import os,sys,shutil
import argparse
import subprocess as sp
from glob import glob
import hashlib

sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+"/common/python")
import cfs

parser=argparse.ArgumentParser()
parser.add_argument('path')
parser.add_argument('obsdir')
parser.add_argument('start')
parser.add_argument('end')
parser.add_argument('-a',action="store_true")
parser.add_argument('-b',action="store_true")
args=parser.parse_args()

args.path = os.path.abspath(args.path)
args.tmpdir   = os.path.abspath(os.getenv("TMP_DIR_LOCAL")+'/atm_omx_'+hashlib.md5(
    args.path+('a' if args.a else 'b')).hexdigest()[0:6])
args.letkfDir = os.path.abspath(os.getenv("CFS_LETKF_ROOT") + '/letkf-gfs')

args.obsdir = os.path.abspath(args.obsdir)
args.ares = 62
args.ares_x,args.ares_y = cfs.getAtmRes(args.ares)
print args

if (not args.a and not args.b) or (args.a and args.b):
    print "must specify either -a or -b"
    sys.exit(1)

    
#get list of files to process
bgfiles = sorted(glob(args.path+"/{}/mean/????/??????/*/*/*_atm.grd".format(
    'gues' if args.b else 'anal')))

#filter out incorrect dates
def filterfunc(d):
    d2= d.split('/')[-1][:10]
    return d2 >= args.start and d2 <= args.end
bgfiles = filter(filterfunc, bgfiles)

# setup working directory
if os.path.exists(args.tmpdir):
    shutil.rmtree(args.tmpdir)
os.makedirs(args.tmpdir)
shutil.copy(args.path+'/letkf.nml', args.tmpdir+'/letkf.nml')

# write spectral configuration file
with open(args.tmpdir+"/ssio.nml", 'w') as f:
    f.write("$common_gfs nlon={} nlat={} gfs_jcap={} nlev=64 /".format(
        args.ares_x, args.ares_y, args.ares))
    

# get a random sig/sfc file to use for later
sigfile = sorted(glob(args.path+'/anal/001/*.sig'))[-1]
sfcfile = sorted(glob(args.path+'/anal/001/*.sfc'))[-1]
print 'using '+sigfile
print 'using '+sfcfile
shutil.copy(sigfile, args.tmpdir+'/template.sig')
shutil.copy(sfcfile, args.tmpdir+'/template.sfc')
    
fcstHrStart = 6
fcstHrEnd = 6
baseSlot =1
window = 0
        
# for each background state
for bg in bgfiles:
    cdate = bg.split('/')[-1][:10]
    print cdate

    # link observations
    idx=1
    for slot in range(-window,window+1):
        slotPfx = '{:+d}'.format(slot) if slot!=0 else ''
        obsFile = (args.obsdir+'/{}/{}/{}/{}/t'+slotPfx+'.dat').format(
            cdate[:4], cdate[:6], cdate[:8], cdate)
        os.symlink(obsFile, args.tmpdir+'/obsin{:02d}.dat'.format(idx))

    #link the background
    shutil.copy(bg, args.tmpdir+'/guesin.grd')

    # convert from gridded to spectral
    os.symlink(args.tmpdir+'/guesin.grd', args.tmpdir+'/fort.41')
    os.symlink(args.tmpdir+'/template.sig', args.tmpdir+'/fort.11')
    os.symlink(args.tmpdir+'/template.sfc', args.tmpdir+'/fort.12')
    shutil.copy(args.tmpdir+'/template.sig', args.tmpdir+'/fort.21')
    shutil.copy(args.tmpdir+'/template.sfc', args.tmpdir+'/fort.22')
    cmd = ("$CFS_LETKF_ROOT/util/bin/grd2ss; "
           "mv fort.21 guesin.sig; "
           "mv fort.22 guesin.sfc; "
           "rm fort.*")
    sp.check_call(cmd, shell=True, cwd=args.tmpdir)

    
    # convert back to grided (annoying.. i know)
    os.symlink(args.tmpdir+'/guesin.sig', args.tmpdir+'/fort.11')
    os.symlink(args.tmpdir+'/guesin.sfc', args.tmpdir+'/fort.12')
    cmd = '$CFS_LETKF_ROOT/util/bin/ss2grd; mv fort.31 gues01.grd'
    sp.check_call(cmd, shell=True, cwd=args.tmpdir)

    # run the atm observation operator
    cmd = '$CFS_LETKF_ROOT/letkf-gfs/letkf/obsope'
    sp.check_call(cmd, shell=True, cwd=args.tmpdir)

    
    # move the file
    cdate = bg.split('/')[-1][:10]    
    outfile =  args.path+'/{}/{}/{}/{}/{}/{}_atm.dat'.format(
         'gues' if args.b else 'anal',
         'omb' if args.b else 'oma',
         cdate[:4], cdate[:6], cdate[:8], cdate)
    print outfile
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)      
    shutil.move(args.tmpdir+'/obsout.dat', outfile)

    #cleanup
    cmd = 'rm *.dat fort.* *.grd guesin.sfc guesin.sig'
    sp.check_call(cmd, shell=True, cwd=args.tmpdir)


shutil.rmtree(args.tmpdir)
