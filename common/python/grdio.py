import  numpy as np


def _procdim(string):
    dim = []
    if (string[2].lower() == "linear"):
        start = float(string[3])
        step = float(string[4])
        num = int(string[1])
        stop = step*num+start
        dim = np.arange(start,stop,step)
    elif (string[2].lower() == 'levels'):
        dim=np.zeros(int(string[1]))
#        for i in string[3:]:
#            dim.append(float(i.replace(',','')))
#        dim=np.array(dim)
    return dim


class GrdDataDescrip:
    name=''
    levels=1
    offset=0
    description=''


class GradsCtl:
    filename=''
    undef = 9.99e33
    endian='little'


    def __init__(self, filename):
        self.data = {}
        self.dataidx=[]
        self.x=[]
        self.y=[]
        self.z=[]       
        self.filename = filename
        with open(filename,'r') as f:
            ctlstr = f.read()
            lines = ctlstr.split('\n')
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
                elif w0 == 'options':
                    for option in ws[1:]:
                        if option.lower() == 'big_endian':
                            self.endian = 'big'
                elif w0 == 'undef':
                    self.undef = float(ws[1])
                elif w0 == 'xdef':
                    self.x = _procdim(ws)
                elif w0 == 'ydef':
                    self.y = _procdim(ws)
                elif w0 == 'zdef':
                    self.z = _procdim(ws)
                elif w0 == 'vars':
                    nvars = int(ws[1])
                elif w0 == 'endvars':
                    nvars = 0
    
def readGrd(ctl, filename):
    grd = {}
    if ctl.endian=='big':
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
            
    if ctl.endian=='big':
        dataout.byteswap(True)               
    dataout.tofile(filename)
