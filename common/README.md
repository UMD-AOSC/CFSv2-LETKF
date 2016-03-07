Common python and fortran code used throughout the system

#### Directory Structure
* `letkf` - fortran code common to observation operators, letkf, etc for the GFS-LETKF and MOM-LETKF domains
* `python` - common python modules, this can be included by specifying `sys.path.insert(1,os.getenv("CFS_LETKF_ROOT")+'/common/python'` within python code