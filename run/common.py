## ------------------------------------------------------------
## ------------------------------------------------------------
## Common functions used by the CFSv2-LETKF control scripts.
## ------------------------------------------------------------
## ------------------------------------------------------------

import logging
import os, sys
from glob import glob


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



############################################################
## determine the number of ensemble members the experiment
## is using based on the contents of the experiment folder.
## simply counts the number of ensemble folders inside PATH/gues
############################################################
def getEnsMem(path):
    dirs = glob(path+'/gues/*/')
    dirs = [d.split('/')[-2] for d in dirs]
    dirs = sorted(filter(lambda x: x.isdigit(),dirs))
    mem = len(dirs)
    assert (int(dirs[-1]) == len(dirs))
    return mem



############################################################
## allow easy logging in the python scripts
## by using log.info(), log.warn(), etc
############################################################
def setupLog():
    log = logging.getLogger('')
    log.setLevel(logging.DEBUG)
    logFormat = logging.Formatter('[%(levelname)-5s %(asctime)s]  %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
    logScreen = logging.StreamHandler(sys.stdout)
    logScreen.setLevel(logging.INFO)
    logScreen.setFormatter(logFormat)
    log.addHandler(logScreen)
    logging.addLevelName(logging.INFO, "\033[01;37mINFO \033[00m")
    logging.addLevelName(logging.ERROR, "\033[01;31mERROR\033[00m")
    logging.addLevelName(logging.WARN, "\033[01;33mWARN \033[00m")
    logging.addLevelName(logging.CRITICAL, "\033[01;35mCRIT \033[00m")
    return log


############################################################
## sets up the logs set up in "setupLog" to be sent out to a file as well
############################################################
def addFileLog(log, path):
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    logFormat = logging.Formatter('[%(levelname)-5s %(asctime)s]  %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
    logFile = logging.FileHandler(filename=path, mode='w')
    logFile.setLevel(logging.DEBUG)
    logFile.setFormatter(logFormat)
    log.addHandler(logFile)
    logging.addLevelName(logging.INFO, "\033[01;37mINFO \033[00m")
    logging.addLevelName(logging.ERROR, "\033[01;31mERROR\033[00m")
    logging.addLevelName(logging.WARN, "\033[01;33mWARN \033[00m")
    logging.addLevelName(logging.CRITICAL, "\033[01;35mCRIT \033[00m")   
