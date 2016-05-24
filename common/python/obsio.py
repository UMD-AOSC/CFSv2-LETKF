import struct
import os

endian = None # 'big', 'little', or None
reclen = 7    # 7 for normal, 10 for post obsop

############################################################
## ob[0] : variable type
## ob[1] : longitude
## ob[2] : latitude
## ob[3] : depth/height
## ob[4] : ob value
## ob[5] : ob error
## ob[6] : ob platform (atm only)

## ob[7] : time relative to analysis
## ob[8] : y_b=H(x)
## ob[9]: 1=valid, 0=invalid
############################################################
def read(inputfile):
    e = ''
    if endian=='little':
        e='<'
    elif endian == 'big':
        e='>'
    strformat='{}4x{}f4x'.format(e,reclen)
    f=open(inputfile,'rb')
    data=f.read()
    f.close()
    parsed=[]
    j=(reclen+2)*4
    for i in range(0,len(data)/j):
        try:
            parsed.append(struct.unpack(strformat,data[i*j:(i+1)*j]))
        except:
            print 'bad end of file'
    return parsed


def write(obs, filename):
    #check to make sure the given path exists, if not create it
    pth=os.path.dirname(os.path.abspath(filename))
    if not os.path.exists(pth):
        os.makedirs(pth)

    e = ''
    if endian=='little':
        e='<'
    elif endian == 'big':
        e='>'
    strformat='{}1i{}f1i'.format(e,reclen)

    #write the file out
    f=open(filename,'wb')
    for ob in obs:
        ob = list(ob)
        args = [strformat,4*reclen]
        for i in range(reclen):
            args.append(ob[i])
        args.append(4*reclen)

        data = struct.pack(*args)
        f.write(data)
    f.close()
