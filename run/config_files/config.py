

## make sure this module is not being run from the command line
import sys
if __name__ == '__main__':
    print "ERROR: This file should not be run from the command line"
    sys.exit(1)

ares = ##atmospheric T resolution (62/126/etc)##
mem  = ##number of ensemble members##
oobs = [ ##list of directories to pull ocean obs from##
       ]
aobs = [ ##list of directories to pull atmospheric obs from##
       ]
    
