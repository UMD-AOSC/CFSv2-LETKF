################################################################################
################################################################################
import sys


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
