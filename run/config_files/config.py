

## make sure this module is not being run from the command line
import sys
if __name__ == '__main__':
    print "ERROR: This file should not be run from the command line"
    sys.exit(1)
