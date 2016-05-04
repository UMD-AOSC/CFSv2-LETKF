#!/usr/bin/env python
import argparse
import datetime as dt
import subprocess as sp
import os, shutil, sys


cfsroot = os.getenv("CFS_LETKF_ROOT")

## get the command line arguments
parser = argparse.ArgumentParser(
    description=("DownloadS the AMSR-E SST L2P (along track) data files"))

parser.add_argument('start', help=(
    'start date of data to get, format is "YYYYMMDD"'))
parser.add_argument('end', help=(
    'end date of data to get, format is "YYYYMMDD"'))

parser.add_argument('--dst', metavar="PATH", default=cfsroot+'/DATA/obs/ocn_amsre', help=(
    'directory to which the files are to be saved, the default is '+cfsroot+'/DATA/obs/ocn_amsre'))

args = parser.parse_args()
args.urlsite = "ftp.nodc.noaa.gov"
args.urldir = "/pub/data.nodc/ghrsst/L2P/AMSRE/REMSS/"
args.start = dt.datetime.strptime(args.start,"%Y%m%d")
args.end = dt.datetime.strptime(args.end,"%Y%m%d")
args.retries = 3

vargs = vars(args)
for v in vargs:
    print v,'=',vargs[v]

import ftplib
import re

ftp = ftplib.FTP(args.urlsite)
ftp.login("anonymous","")

cdate = args.start
while(cdate <= args.end):
    dy = (cdate - dt.datetime(cdate.year,1,1)).days + 1
    yr = cdate.year
    print ""
    print "processing: ",cdate,"  year/day: ",yr,"/",dy

    urldir=args.urldir+"/{}/{:03d}/".format(yr,dy)
    print urldir
    try:
        files = ftp.nlst(urldir)
    except:
        print "ERROR: directory not found:  "+urldir
        cdate += dt.timedelta(days=1)    
        continue
    
    i = 0
    for f in files:
        if re.match(".*AMSRE\-REMSS\-L2P\-amsr\_l2b\_v05", f):
            i += 1            
            outfilename = args.dst+('/{0}/{0}{1:02d}/{0}{1:02d}{2:02d}/'
                                    '{0}{1:02d}{2:02d}_{3:03d}.nc.gz').format(
                                        cdate.year,cdate.month,cdate.day,i)
            outpath = os.path.dirname(outfilename)
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            print "Downloading ",f
            retry = 0
            downloaded = False
            while (retry < args.retries and not downloaded):
                try:
                    with open(outfilename, 'wb') as output:
                        ftp.retrbinary('RETR %s' % f, output.write)
                        downloaded = True
                except:
                    print "retrying a download"
                    retry += 1
            if retry == args.retries:
                print "[ERROR] there was a problem downloading, giving up"
                sys.exit(1)
                        
    sp.call("gunzip *.gz -f", shell=True, cwd=outpath)

    cdate += dt.timedelta(days=1)
