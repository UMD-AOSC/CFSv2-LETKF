################################################################################
################################################################################
import sys
import os
import netCDF4 as nc
import numpy as np
import scipy.spatial

_root = os.getenv("CFS_LETKF_ROOT")

## make sure this module is not being run from the command line
if __name__ == '__main__':
    logging.critical("This file should not be run from the command line")
    sys.exit(1)

## list of available GFS resolutions
aresList = [62,126,190,382,574,1148]
    
## Get the x/y atmosphere grid resolution base on the
## spectral resolution given
def getAtmRes(t_res):
    res={
        62   : ( '192',   '94'),
        126  : ( '384',  '190'),
        190  : ( '576',  '288'),
        382  : ('1152',  '576'),
        574  : ('1760',  '880'),
        1148 : ('2304', '1152')}
    assert(t_res in res)
    return res[t_res]

## gets the x/y coordinate on the ocean grid given the
## latitude and longitude
_ocnGrd_KDTree = None
_ocnGrd_lons = None
_ocnGrd_lats = None

def getOcnXY(lon, lat):
    global _ocnGrd_KDTree
    global _ocnGrd_lons
    global _ocnGrd_lats
    global _ocnGrd_xs
    global _ocnGrd_ys
    global _ocnGrd_lon_bounds
    if _ocnGrd_KDTree == None:
        print 'Loading ocean grid into KD Tree...'
        ncdat = nc.Dataset(_root+"/support/fix/fix_om/grid_spec_05.nc.T62")
        _ocnGrd_lons = ncdat.variables["x_T"][:]
        _ocnGrd_lats = ncdat.variables["y_T"][:]
        _ocnGrd_lon_bounds = ( _ocnGrd_lons[0][0],_ocnGrd_lons[0][-1])
        _ocnGrd_KDTree = scipy.spatial.KDTree(zip(_ocnGrd_lons.reshape(-1),
                                                  _ocnGrd_lats.reshape(-1)))
        xs, ys = np.meshgrid(np.arange(0,_ocnGrd_lons.shape[1]),
                             np.arange(0,_ocnGrd_lons.shape[0]))
        _ocnGrd_xs = xs.reshape(-1)
        _ocnGrd_ys = ys.reshape(-1)
        
        ncdat.close()

    if lon > _ocnGrd_lon_bounds[1]:
        lon -= 360
    if lon < _ocnGrd_lon_bounds[0]:
        lon += 360
        
    res = _ocnGrd_KDTree.query( (lon,lat) )
    y = _ocnGrd_ys[res[1]]
    x = _ocnGrd_xs[res[1]]
    return x,y
    

def getAtmLevels():
    lvlFile = _root+"/support/fix/fix_am/global_hyblev.l64.txt"
    if getAtmLevels.lvls == None:
        f = open(lvlFile,'r')
        d = f.read().splitlines()[1:]
        params = []
        for l in d:
            params.append([float(x) for x in l.split() ])
        getAtmLevels.lvls = params
    return getAtmLevels.lvls
getAtmLevels.lvls = None
