#!usr/bin/env python
#===============================================================================
#
#       AUTHOR              : Kriti Bhargava
#       DATE WRITTEN        : Wed Aug  7 15:29:14 EDT 2019
#       LAST MODIFIED       : Wed Aug  7 15:29:14 EDT 2019
#       PURPOSE             : This code calculates the analysis increments 
#                             from the atmos and ocn grd files by first
#                             converting them to netcdf and then taking 
#                             the difference. It gives the user a choice 
#                             to save the analysis increments for all the 
#                             dates specified or just the average in case
#                             user wants to estimate the average and just 
#                             keep that.
#       USAGE               : .cal_ai.py [-h] [--cal_avg] [--store_all]
#                    startdate enddate exp {ocn,atm}
#   Claculate analysis increments for specific dates and if specified calculate
#   the average too.
#   
#   positional arguments:
#     startdate    START DATE in YYYYMMDDHH format
#     enddate      END DATE in YYYYMMDDHH format
#     exp          Experiment name for which analysis increments to be computed.
#     {ocn,atm}    Domain for calculating analysis increment ['ocn','atm']
#   
#   optional arguments:
#     -h, --help   show this help message and exit
#     --cal_avg    If specified calculates the temporal average analysis
#                  increments for the dates specified.
#     --store_all  If specified the option saves all the analysis increments file,
#                  if not only the average is kept if cal_avg specified.
#                
#===============================================================================

import argparse
import os
from tools_nc_plot import get_grd_to_nc, valid_date, get_mean_nc
import sys
import subprocess as sp
import datetime as dt
# get the CFS_LETKF_ROOT PATH (it's the environment variable set up already
cfsroot = os.getenv("CFS_LETKF_ROOT")

# Define argument parser
parser = argparse.ArgumentParser(description=("Claculate analysis increments for \
    specific dates and if specified calculate the average too."))
# Add arguments to the parser
parser.add_argument("startdate",type = valid_date,
    help=("START DATE in YYYYMMDDHH format"))
parser.add_argument("enddate", type= valid_date,
    help=("END DATE in YYYYMMDDHH format"))
parser.add_argument("exp",
    help=("Experiment name for which analysis increments to be computed."))
parser.add_argument("domain", choices=["ocn","atm"], 
    help=("Domain for calculating analysis increment ['ocn','atm']"))
parser.add_argument("--cycle", choices=[0,6,12,18],type=int,
    help=("Cycle for calculating analysis increment. Assumes all cycles if not specified"))
parser.add_argument("--cal_avg", action= "store_true", default = False,
    help=("If specified calculates the temporal average analysis increments for the dates specified."))
parser.add_argument("--store_all", action= "store_true", default = False,
    help=("If specified the option saves all the analysis increments file, \
        if not only the average is kept if cal_avg specified."))
args = parser.parse_args()

# Printing the argument list on the screen
print ("\n")
print("*"*31,"INPUT ARGUEMENTS","*"*31)
for arg in vars(args):
    print ("{0}\t\t\t{1}".format(arg.upper(),getattr(args,arg)))
print ("*"*80)

# Setting the date for proper cycle
cdate = args.startdate
if (args.cycle is None):
    time_inc = 6 #increase date by 6 hours for all 4 cycles
    cycle="all"
else:
    args.startdate = args.startdate.replace(hour=args.cycle)
    print ("UPDATING startdate in accordance with the cycle provided...")
    print ("\tUpdated start date is now {}.".format(cdate.strftime("%Y%m%d%H")))
    time_inc = 24
    cycle="{0}z".format((str(args.cycle)).zfill(2))
start_string = args.startdate.strftime("%Y%m%d%H")
end_string = args.enddate.strftime("%Y%m%d%H")
cdate=args.startdate
print ("\tTime increment set to {} hours".format(str(time_inc)))
# Converting to nc and calculating AI
if (args.cal_avg):
    anal_list = []
    gues_list = []
    ai_list = []
print ("*"*80)
print ("STARTING TO CONVERT TO NC AND CALCULATE ANALYSIS INCREMENT")
print ("*"*80)
while (cdate < args.enddate):
    # Getting the nc anal and guess files
    anal_out = get_grd_to_nc(domain=args.domain,cdate=cdate, exp_name=args.exp,
                  out_type="anal",ens="mean") #converts analysis grd file to nc
    gues_out = get_grd_to_nc(domain=args.domain,cdate=cdate, exp_name=args.exp,
                  out_type="gues",ens="mean") #converts guess grd file to nc
    anlinc_path = os.path.join(cfsroot,"DATA/{0}/nc_files/anlinc".format(args.exp))
    if not(os.path.exists(anlinc_path)):
        os.mkdir(anlinc_path)
        print ("\tCREATED DIRECTORY {}".format(anlinc_path))
    anlinc_file = os.path.join(anlinc_path,
                    "anlinc_{0}_{1}_{2}.nc".format(cdate.strftime("%Y%m%d%H"),args.domain,cycle))
    
    # Calulating analysis increment
    sp.call(["cdo sub {0} {1} {2}".format(anal_out, gues_out, anlinc_file)],
            shell = True) #UPDATE this
    print ("\t Analysis Increment calculated for {}.".format(cdate.strftime("%Y%m%d%H")))
    
    # Appending files to file list to get the average of
    if (args.cal_avg):
        anal_list.append(anal_out)
        gues_list.append(gues_out)
        ai_list.append(anlinc_file)
    cdate+=dt.timedelta(hours=time_inc)

# Calculating averages
if (args.cal_avg):
    print ("*"*80)
    print ("CALCULATING AVERAGES NOW")
    print ("*"*80)
    anal_mn_file = os.path.join(anlinc_path,
                    "anal_mn_{0}_{1}_{2}.nc".format(start_string,end_string,cycle))
    gues_mn_file = os.path.join(anlinc_path,
                    "gues_mn_{0}_{1}_{2}.nc".format(start_string,end_string,cycle))
    ai_mn_file = os.path.join(anlinc_path,
                    "ai_mn_{0}_{1}_{2}.nc".format(start_string,end_string,cycle))
    get_mean_nc(anal_list, anal_mn_file,args.exp) # Calculating mean analysis
    get_mean_nc(gues_list, gues_mn_file,args.exp) # Calculating mean analysis
    get_mean_nc(ai_list,   ai_mn_file  ,args.exp) # Calculating mean analysis

# Delete files:    
if ((args.cal_avg) and not(args.store_all)):
    print ("*"*80)
    print ("DELETING FILES NOW")
    print ("*"*80)
    for afile in anal_list:
        os.remove(afile)
    print ("\tAnalysis nc files removed for {0} to {1}.".format(start_string,end_string))
    for gfile in gues_list:
        os.remove(gfile)
    print ("\tGuess nc files removed for {0} to {1}.".format(start_string,end_string))
    #for aifile in ai_list:
    #    os.remove(aifile)

