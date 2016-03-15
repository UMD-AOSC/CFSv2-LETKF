import  numpy as np
import os
import re


_NUMBER = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'


class GrdDataDescrip:
    name=''
    levels=1
    offset=0
    description=''


class GradsCtl:
   
    def __init__(self, filename):
        self.undef = 9.99e33
        self.big_endian = False

        self.data = {}
        self.dataidx=[]
        self.x=None
        self.y=None
        self.z=None     
        self.filename = filename

        ## read in the file
        with open(filename,'r') as f:
            self.ctl = f.read()
        self.ctlU = self.ctl.upper()

        ## process options
        if "UNDEF" in self.ctlU:
            self.undef = eval(re.search("UNDEF (%s)" % _NUMBER, self.ctlU).group(1))
        self.big_endian = bool(re.search("OPTIONS.*BIG_ENDIAN", self.ctlU))
        self.yrev = bool(re.search("OPTIONS.*YREV", self.ctlU))
        self.dset = re.search("DSET (.*)", self.ctlU).group(1)
        if self.dset.startswith('^'):
            self.dset = os.path.join(os.path.dirname(self.filename), self.dset[1:])

            
        ## process x & y dimensions
        if "XDEF" in self.ctlU:
            self.x = self._procdim("XDEF")
        if "YDEF" in self.ctlU:
            self.y = self._procdim("YDEF")
        if "ZDEF" in self.ctlU:
            self.z = self._procdim("ZDEF")
        # if "TDEF" in self.ctlU:
        #     self.t = self._procdim("TDEF")

        ## process variables
        self._procvars()

            
    ## ------------------------------
    def _procdim(self, dim):

        ## if linear array is defined
        p = re.compile("%s\s+(\d+)\s+LINEAR\s+(%s)\s+(%s)" % (dim, _NUMBER, _NUMBER))
        m = p.search(self.ctlU)
        if m:
            length = int(m.group(1))
            start = float(m.group(2))
            increment = float(m.group(3))
            return np.arange(start, start+length*increment, increment)   

        ## levels are defined
        p = re.compile("%s\s+\d+\s+LEVELS((\s+%s(,)?)+)" % (dim, _NUMBER))
        m = p.search(self.ctlU)
        if m:
            return np.fromstring(m.group(1).replace(',',' '), sep=' ')

        ## and, whatever this one if for
        p = re.compile("%s\s+(\d+)\s+LINEAR\s+([:\w]+)\s+(\d{1,2})(\w{2})" % dim)
        m = p.search(self.ctlU)
        if m:
            length = int(m.group(1))
            start = parse_date(m.group(2))
            increment = parse_delta(m.group(3), m.group(4))
            return np.array([ start+i*increment for i in range(length)]) 

    def _procvars(self):
        lines = self.ctl.split('\n')
        nvars = 0
        varoffset=0
        for l in lines:
            ws = filter(None, l.split(' '))
            if len(ws) == 0:
                continue
            w0 = ws[0].lower()
            if nvars > 0:
                nvars = nvars-1
                d=GrdDataDescrip()
                d.name = ws[0].lower()
                d.levels = int(ws[1])
                if (d.levels == 0):
                    d.levels = 1
                d.offset=varoffset
                d.ctl=self
                self.data[w0]=d
                self.dataidx.append(d)
                varoffset=varoffset+len(self.x)*len(self.y)*d.levels
            elif w0 == 'vars':
                nvars = int(ws[1])
            elif w0 == 'endvars':
                nvars = 0
    
def readGrd(ctl, filename):
    grd = {}
    if ctl.big_endian:
        fmt = '>f4'
    else:
        fmt = '<f4'
    data=np.fromfile(filename, dtype=fmt)
    recs=len(data)/len(ctl.x)/len(ctl.y)
    for v in ctl.data:
        d=ctl.data[v]
        sz=d.levels*len(ctl.x)*len(ctl.y)        
        grd[v]=np.reshape(data[d.offset:(d.offset+sz)],(d.levels,len(ctl.y),len(ctl.x)))
        grd[v][grd[v]==ctl.undef] = np.nan
    return grd


def writeGrd(ctl, data, filename):
    #write out the variables in the order required by the ctl file
    first=True
    for dv in ctl.dataidx:
        if not dv.name in data:
            raise Exception("Can't find '"+dv.name+"' required by "+ctl.filename+" while writing out grd file, check the input data")
        #TODO
        #check to make sure the data is the right size
        if (first):
            dataout=data[dv.name].flatten()
            first = False
        else:
            dataout=np.append(dataout, data[dv.name].flatten())
    dataout[dataout==np.nan]=ctl.undef
            
    if ctl.big_endian:
        dataout.byteswap(True)               
    dataout.tofile(filename)
