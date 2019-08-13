#!/usr/bin/env python
#==============================================================================#
#
#           AUTHOR          : Kriti Bhargava
#           DATE WRITTEN    : Fri Jul 26 13:58:23 EDT 2019
#           LAST MODIFIED   : Fri Aug  9 14:16:19 EDT 2019
#           PURPOSE         : Plot ocean data from .nc files for a given variable
#                             at a given depth.
#
#==============================================================================#

import argparse
import subprocess as sp
import os
import datetime as dt
import re
import sys
import numpy as np
import cartopy as cart
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

#import tempfile
cfsroot = os.getenv("CFS_LETKF_ROOT")
#tmp_dir = os.path.join(cfsroot, "tmp")
def_date=dt.datetime(2006, 6, 1, 6)

def valid_date(s):
    try:
        return dt.datetime.strptime(s, "%Y%m%d%H")
    except ValueError:
        msg = "Not a valid date: '{0}'. Please use YYYYMMDDHH format.".format(s)
        print (msg)
        sys.exit()
        #raise TypeError(msg)

def get_grd_to_nc(domain="ocn",cdate=def_date,exp_name="exp01", store_ctl=False,out_type="anal",ens="mean"):
    # Create directory to store nc_files
    nc_dir =os.path.join(cfsroot,"DATA/{0}/nc_files/{1}".format(exp_name,out_type))
    if not(os.path.exists(nc_dir)):  #check if the directory exists
        os.makedirs(nc_dir)
    
    # Full path for the template ctl file and other filenames    
    date_string = cdate.strftime("%Y%m%d%H")  #date converted to string in YYYMMDDHH format
    path_ctl = os.path.join(cfsroot,
                "DATA/{0}/{1}/{2}".format(exp_name,out_type,ens)) #path to ctl template
    ctl_file=os.path.join(path_ctl,"{}.ctl".format(domain)) # full filename for ctl template
    out_ctl = os.path.join(path_ctl,
                "{0}_{1}.ctl".format(date_string,domain)) #ctl_file
    out_file = os.path.join(nc_dir,
                "{0}_{1}.nc".format(date_string,domain)) #output nc filename
    if (os.path.exists(out_file)):
        print ("\tThe nc file for {0} at {1} already exists.".format(out_type,date_string))
    else:
    # Write the ctlfile
        with open(ctl_file) as f:
            file_str = f.read()
        t_sub="1 linear {}".format(cdate.strftime("%d%b%Y"))
        file_str=re.sub("%y4",cdate.strftime("%Y"),file_str)
        file_str=re.sub("%m2",cdate.strftime("%m"),file_str)
        file_str=re.sub("%d2",cdate.strftime("%d"),file_str)
        file_str=re.sub("%h2",cdate.strftime("%H"),file_str)
        file_str=re.sub("OPTIONS template","",file_str)
        if domain ==("atm"):
            file_str=re.sub("10000 linear 1Jun2006",t_sub,file_str)
        else:
            file_str=re.sub("1000 LINEAR  1Jun2006",t_sub,file_str)
        f = open(out_ctl, "w")
        f.write(file_str)
        f.close()
        
        # call cdo and convert file
        status = sp.call(["cdo -f nc import_binary {0} {1}".format(out_ctl,out_file)],shell=True)
        if status ==0 and not(store_ctl):
           os.remove(out_ctl)
    return out_file

def get_mean_nc(file_list, out_file,exp_name="exp01"):
    # Change directory to directory containing nc files
    nc_dir =os.path.join(cfsroot,"DATA/{}/nc_files".format(exp_name))
    os.chdir(nc_dir)
    # Check if all files exist
    for f in file_list:
        if not(os.path.isfile(f)):
            print ("The file names {0} doesn't exist in the directory {1}".format(f,nc_dir))
            print ("Aborting now...")
            sys.exit()
    # Call cdo and convert file
    file_list=" ".join(file_list)
    status = sp.call(["cdo ensmean {0} {1}".format(file_list,out_file)],shell=True)

def plot_map(lons,lats,data, title,figname,proj="ccrs.PlateCarree()",cmap_type="default"):
    #TODO define function that creates subplots instead of the fig
    cmap, clev, tick_list,tick_label_list =get_cmap_tick(cmap_type)
    fig = plt.figure(figsize=(6,3))
    ax = plt.axes(projection = ccrs.Robinson(central_longitude=0))
    cf = plt.contourf(lons, lats,data,transform=ccrs.PlateCarree(),levels=clev,cmap=cmap,extend="both")
    ax.coastlines()
    gl=ax.gridlines(crs=ccrs.PlateCarree())
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0, 30, 60, 90])
    ax.set_title(title)
    ax.add_feature(cart.feature.LAND,zorder=100,edgecolor='k',facecolor="w")
    cb=fig.colorbar(cf, ax=ax,extend="both",ticks=tick_list,shrink=0.8)
    cb.ax.set_yticklabels(tick_label_list)
    plt.savefig(figname,dpi=600)
    plt.show()
    
def plot_zonal_mn(levs,lats,data, title,figname,cmap_type="default"):
    #TODO define function that creates subplots instead of the fig
    cmap, clev, tick_list,tick_label_list =get_cmap_ticks(cmap_type)
    fig = plt.figure(figsize=(6,3))
    ax =plt.axes()
    cf = plt.contourf(lats,levs,data,levels=clev,cmap=cmap,extend="both")
    ax.set_title(title)
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Ocean Depth")
    cb=fig.colorbar(cf, ax=ax,extend="both",ticks=tick_list,shrink=0.8)
    cb.ax.set_yticklabels(tick_label_list)
    plt.tight_layout()
    plt.savefig(figname,dpi=600)
    plt.show()

def plot_global_mn(levs,data, title,figname):
    #TODO define function that creates subplots instead of the fig
    print (levs)
    fig = plt.figure(figsize=(6,3))
    ax = plt.axes()
    pt = plt.plot(data,levs)
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Global Mean")
    ax.set_ylabel("Ocean Depth")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(figname,dpi=600)
    plt.show()

def get_cmap_ticks(cmap_type):
    if (cmap_type == "default"):
        cmap = plt.cm.jet
        clev = np.arange(-5,30.1,0.25)
        tick_list = [-5,0,5,10,15,20,25,30]
        tick_label_list = ["-5","0","5","10","15","20","25","30"]
    elif (cmap_type == "diverging"):
        num=1
        cmap = mcolors.LinearSegmentedColormap.from_list('err1',['#19204c', '#36459c','#a6d5e7','#ffffff','#f7844e','#a50026','#d405be'])
        cmap.set_under('#0c1026')
        cmap.set_over('#590014')
        clev = np.arange(-num,num+0.1,0.05)
        tick_list = np.arange(-num,num+1,0.5)
        #tick_label_list = ["-5","-3","-1","0","1","3","5"]
        tick_label_list = map(str,tick_list)
        #clev = np.arange(-5,5.1,0.25)
        #tick_list = [-5,-3,-1,0,1,3,5]
        #tick_label_list = ["-5","-3","-1","0","1","3","5"]
        return cmap, clev, tick_list,tick_label_list
