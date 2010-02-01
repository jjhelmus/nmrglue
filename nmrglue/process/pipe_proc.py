"""
Function to process NMRPipe dic,data tuples similar to NMRPipe

pipe_proc
=========

Provides:
    
    1.  NMR Processing functions similar to NMRPipe.

Documentation is available in the docstrings and at http://XXX


To do:  
        test with 1D/2D and real/complex and freq/time and transposed/not
        Check for XXX in this file
        extract out more into proc_base
        unify recomputer orig functions
        Lots of documentation.

NMRPipe Functions Not Implemented (11) will raise NotImplemented exception:

    poly    - polynonial baseline correction (very little documentation)
    lp      - linear prediction
    lpc     - macro of lp
    lp2d    - 2D linear prediction
    mem     - Maximum Entropy
    ml      - maximum likelyhood freq. map (no documentation)
    ebs     - (no documentation)
    ann     - (no documentation)
    ztp     - 3D matrix transpose (no need)
    xyz2zyx - 3D matrix transpose (no need)
    mac     - macro language interpreter (no need)

Implementation Caveat:

    - hdr=True overwrites apodization values in the calling function.
    - To delete imaginaries use the di function (-di not implemented).
    - All functions require x1, xn, etc to be in units of points, use
      the spec2pts to calculate these values before calling.
    - No functions implement the dmx or nodmx flags.

History
(jjh) 2009.10.28 changed to unit_conversion object framework
(jjh) 2009.10.14 changes to pnt2spec and spec2pnt
(jjh) 2009.10.02 code refactoring
(jjh) 2009.09.18 unit conversion functions, bugfixes in ext
(jjh) 2009.09.14 sol and med functions
(jjh) 2009.09.08 base function.
(jjh) 2009.08.XX various bug fixes as making test suite    
(jjh) 2009.07.15 implemented di
(jjh) 2009.06.12 Fixed bugs in apod functions.
(jjh) 2009.06.05 finished zd, cleaned up code
(jjh) 2009.06.04 Implemented qmix,zd
(jjh) 2009.06.03 Documentation now listed unimplemented functions
(jjh) 2009.05.28 many functions implemented in last 2 weeks
(jjh) 2009.05.12 jmod function
(jjh) 2009.05.11 em, gm, gmb function

"""

# standard library modules

# external modules
import numpy as np
import scipy
import scipy.signal

# nmrglue modules
from nmrglue.fileio import pipe,fileiobase
import proc_base as p
import proc_bl


__version__ = '0.98'

pi = np.pi

###################
# Unit conversion #
###################

class unit_conversion(fileiobase.unit_conversion):
    
    # NMRPipe indexes from 1 to MAX instead on 0 to MAX-1
    # we need to modify two method to account for this off by one problem
    
    def __unit2pnt(self,val,units):
        return fileiobase.unit_conversion.__unit2pnt(self,val,units)+1

    def __pnt2unit(self,val,units):
        return fileiobase.unit_conversion.__pnt2unit(self,val-1,units)



def make_uc(dic,data,dim=-1):
    """ Make a unit conversion object"""

    if dim == -1:
        dim = data.ndim - 1 # last dimention

    fn = "FDF" + str(int(dic["FDDIMORDER"][data.ndim-1-dim]))
    size = float(data.shape[dim])

    # check for quadrature in indirect dimentions
    if (dic[fn+"QUADFLAG"] != 1) and (dim !=data.ndim-1):
        size = size/2.
        cplx = True
    else:
        cplx = False

    sw = dic[fn+"SW"]
    if sw == 0.0:
        sw = 1.0
    obs = dic[fn+"OBS"]
    if obs == 0.0:
        obs = 1.0
    car = dic[fn+"CAR"]*obs

    return unit_conversion(size,cplx,sw,obs,car)



#####################
# Dictionary Macros #
#####################


def update_minmax(dic,data):
    """
    Update the FDDISPMAX parameters returns dictionary

    """

    # maximum and minimum values
    dic["FDMAX"] = float(data.max().real)
    dic["FDDISPMAX"] = dic["FDMAX"]
    dic["FDMIN"] = float(data.min().real)
    dic["FDDISPMIN"] = dic["FDMIN"]
    dic["FDSCALEFLAG"] = 1.0    # FDMIN/MAX are valid

    return dic


def clean_minmax(dic):
    """
    clear FDDISPMAX parameters returns dictionary
    
    """

    # maximum and minimum values
    dic["FDMAX"] = 0.0
    dic["FDDISPMAX"] = 0.0
    dic["FDMIN"] = 0.0
    dic["FDDISPMIN"] = 0.0
    dic["FDSCALEFLAG"] = 0.0    # FDMIN/MAX not valid

    return dic


#########################
# Apodization functions #
#########################


def apod(dic,data,qName=None,q1=0.0,q2=0.0,q3=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """ 
    Generic Apodization
    """
    # calls the appropiate apodization function

    a_list = ['SP','EM','GM','GMB','TM','TRI','JMOD']

    if hdr:
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        qnum = dic[fn+"APODCODE"]
        qName = ["","SP","EM","GM","TM","","TRI","GMB","JMOD"][qnum]
    
    if qName not in a_list:
        raise ValueError("qName must be SP, EM, GM, GMB, TM, TRI or JMOD")

    if qName == "EM":
        return em(dic,data,q1,c,start,size,inv,one,hdr) 

    if qName == "GM":
        return gm(dic,data,q1,q2,q3,c,start,size,inv,one,hdr) 

    if qName == "GMB":
        return gmb(dic,data,q1,q2,c,start,size,inv,one,hdr)

    if qName == "JMOD":
        return jmod(dic,data,q1,q2,q3,False,False,c,start,size,inv,one,hdr)

    if qName == "SP":
        return sp(dic,data,q1,q2,q3,c,start,size,inv,one,hdr)

    if qName == "TM":
        return tm(dic,data,q1,q2,c,start,size,inv,one,hdr)

    if qName == "TRI":
        return tri(dic,data,q1,q2,q3,c,start,size,inv,one,hdr)


def em(dic,data,lb=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    Exponential Multiply Window

    """
    start = start - 1   # arrays should start at 0

    # update dictionary
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data headers
        c = dic[fn+"C1"] + 1
        lb = dic[fn+"APODQ1"]

    dic[fn+"C1"] = c - 1.0

    # set the apod flags
    dic[fn+"APODCODE"] = 2.0
    dic[fn+"APODQ1"] = lb

    sw = dic[fn+"SW"] 
    flb = lb/sw

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.em(data,lb=flb,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        data[...,start:stop] = p.em(data[...,start:stop],lb=flb,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


def gm(dic,data,g1=0.0,g2=0.0,g3=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    Lorentz-to-Gauss apodization

    """
    start = start - 1 # arrays should start at 0

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data header
        g1 = dic[fn+"APODQ1"]
        g2 = dic[fn+"APODQ2"] 
        g3 = dic[fn+"APODQ3"] 
        c = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 3.0
    dic[fn+"APODQ1"] = g1
    dic[fn+"APODQ2"] = g2
    dic[fn+"APODQ3"] = g3

    # calculate native parameters
    sw = dic[fn+"SW"]
    g1p = g1/sw
    g2p = g2/sw
    g3p = g3

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.gm(data,g1p,g2p,g3p,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        
        # pipe sets the maximum to the actual data maximum not
        # the maximum of the windowed region, so adj. g3p as necessary
        g3p = g3p*data.shape[-1]/(stop-start)
        #print start,stop,g1p,g2p,g3p
        data[...,start:stop] = p.gm(data[...,start:stop],g1p,g2p,g3p,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


def gmb(dic,data,lb=0.0,gb=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    Gaussian apodization

    """
    start = start - 1 # arrays should start at 0

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data header
        lb = dic[fn+"APODQ1"]
        gb = dic[fn+"APODQ2"] 
        c = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 7.0
    dic[fn+"APODQ1"] = lb
    dic[fn+"APODQ2"] = gb

    # calculate native parameters
    sw = dic[fn+"SW"]
    a = pi*lb/sw
    b = -a/ (2.0*gb*data.shape[-1])

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.gmb(data,a,b,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        data[...,start:stop] = p.gmb(data[...,start:stop],a,b,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


def jmod(dic,data,off=0.0,j=0.0,lb=0.0,sin=False,cos=False,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    exponentially damped j-modulation

    sin or cos set True overwrite off value.
    """
    start = start - 1 # arrays should start at 0

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc


    if sin:
        off = 0.0
    if cos:
        off = 0.5

    if hdr: # read apod values from data header
        off = dic[fn+"APODQ1"]
        j   = dic[fn+"APODQ2"]
        lb   = dic[fn+"APODQ3"]
        c  = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 8.0
    dic[fn+"APODQ1"] = off
    dic[fn+"APODQ2"] = j
    dic[fn+"APODQ3"] = lb

    # calculate native parameters
    sw = dic[fn+"SW"]
    e   = pi*lb/sw 
    end = off+j*(data.shape[-1]-1)/sw

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.jmod(data,e,off,end,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        #print start,stop,e,off,end
        end = off+j*(stop-start-1)/sw
        data[...,start:stop] = p.jmod(data[...,start:stop],e,off,end,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


def sp(dic,data,off=0.0,end=1.0,pow=1.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    sine bell
    """
    start = start - 1 # arrays should start at 0

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data header
        off = dic[fn+"APODQ1"]
        end = dic[fn+"APODQ2"]
        pow = dic[fn+"APODQ3"]
        c   = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 1.0
    dic[fn+"APODQ1"] = off
    dic[fn+"APODQ2"] = end
    dic[fn+"APODQ3"] = pow

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.sp(data,off,end,pow,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        data[...,start:stop] = p.sp(data[...,start:stop],off,end,pow,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


sine = sp   # wrapper for sine functions 


def tm(dic,data,t1=0.0,t2=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    trapezoid apodization
    """
    start = start - 1 # arrays should start at 0

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data header
        t1 = dic[fn+"APODQ1"]
        t2 = dic[fn+"APODQ2"]
        c   = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 4.0
    dic[fn+"APODQ1"] = t1
    dic[fn+"APODQ2"] = t2

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.tm(data,t1=t1,t2=t2,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        data[...,start:stop] = p.tm(data[...,start:stop],t1=t1,t2=t2,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0


    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    # check for NaN in array (when div by 0)
    if np.isnan(data).any():    
        data = np.array(np.nan_to_num(data),dtype=data.dtype)

    dic = update_minmax(dic,data)

    return dic,data


def tri(dic,data,loc="auto",lHi=0.0,rHi=0.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """
    triangular apodization

    Differences from NMRPipe.
    Right side of spectrum slightly different, but errors are small,

    """
    start = start - 1 # arrays should start at 0

    if loc == "auto":
        loc = data.shape[-1]/2

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if hdr: # read apod values from data header
        loc = dic[fn+"APODQ1"]
        lHi = dic[fn+"APODQ2"]
        rHi = dic[fn+"APODQ3"]
        c   = dic[fn+"C1"] + 1

    # update the dictionary
    dic[fn+"C1"] = c - 1.0
    
    # set the apod flags
    dic[fn+"APODCODE"] = 6.0
    dic[fn+"APODQ1"] = loc
    dic[fn+"APODQ2"] = lHi
    dic[fn+"APODQ3"] = rHi

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.tri(data,loc=loc,lHi=lHi,rHi=rHi,inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start+size
        data[...,start:stop] = p.tri(data[...,start:stop],loc,lHi,rHi,inv=inv)
        if one == False:
            data[...,:start] = 0.0
            data[...,stop:] = 0.0

    # first point scaling
    if inv:
        data[...,0] = data[...,0]/c
    else:
        data[...,0] = data[...,0]*c

    dic = update_minmax(dic,data)

    return dic,data


###################
# Shift functions #
###################

def rs(dic,data,rs=0.0,sw=False):
    """
    Right shift and zero pad
    """

    if rs < 0:  # negative right shifts are left shifts
        return ls(dic,data,ls=-rs,sw=sw)

    data = p.rs(data,pts=rs)
    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if sw and dic[fn+"FTFLAG"] == 1:    
        # we are in freq domain and must update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]+rs
        
        # recalc orig
        s = data.shape[-1]
        s2 = dic[fn+"CENTER"]
        sw = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
        dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs
            
    return dic,data

def ls(dic,data,ls=0.0,sw=False):
    """
    Left shift and zero pad
    """
    if ls < 0:
        return rs(dic,data,rs=-ls,sw=sw)

    data = p.ls(data,ls)
    dic = update_minmax(dic,data)


    if sw:
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        if dic[fn+"FTFLAG"] == 1:   # freq domain
            # update NDORIG and NDCENTER
            dic[fn+"CENTER"] = dic[fn+"CENTER"]-ls

            # recalc orig from new center
            s = data.shape[-1]
            s2 = dic[fn+"CENTER"]
            sw = dic[fn+"SW"]
            car = dic[fn+"CAR"]
            obs = dic[fn+"OBS"]
            dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs

        else:   # time domain
            dic[fn+"APOD"] = data.shape[-1] - ls
            dic[fn+"TDSIZE"] = data.shape[-1] - ls 

    return dic,data

def cs(dic,data,dir,pts=0.0,neg=False,sw=False):
    """
    Circular Shift

    Differences from NMRPipe:
    Functions syntax is different.  Parameter dir should be 'rs' or 'ls'.
    Parameter pts is the number of points to shift.
    """
    
    if dir == "ls":
        pts = -pts
    elif dir !="rs":
        raise ValueError("dir must be ls or rs")

    data = p.roll(data,pts,neg=neg)
    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if sw and dic[fn+"FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]+pts

        # recalc orig
        s = data.shape[-1]
        s2 = dic[fn+"CENTER"]
        sw = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
        dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs

    return dic,data

def fsh(dic,data,dir,pts,sw=True):
    """
    frequency shift via FT

    Differences from NMRPipe:
    - Parameter sw is on by default.
    - Does not perform Hilbert transform when data is complex (NMRPipe does) 
      rather uses provided imaginary data.
    - FDMIN and FDMAX slightly off (results of double precision FFTs? jjh)

    """

    if dir not in ["ls","rs"]:
        raise ValueError("dir must be ls or rs")

    if np.iscomplexobj(data) == False:  # real data
        null,data = _ht(dict(dic),data,zf=True)
        del_imag = True

    else:   # imaginary data
        del_imag = False
        # NMRPipe always performs a hilbert transform
        # uncommenting the next two lines will match NMRPipe's fsh real
        # channel results, no idea how to get the imaginary channel
        #null,data = _ht(dict(dic),data,zf=True)
        #data = np.array(data,dtype="complex64")

    if dir == "ls":
        pts = -pts

    data = p.fsh(data,pts)

    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    if dic[fn+"FTFLAG"] == 1 and sw: # freq domain
        
        dic[fn+"CENTER"] = dic[fn+"CENTER"] + pts

        # update NDORIG
        s = dic["FDSIZE"]
        sw  = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
        center = dic[fn+"CENTER"]
        dic[fn+"ORIG"] = -sw * (s-center)/s + car * obs


    if del_imag == False:
        return dic,data
    else:
        return dic,data.real


##############
# Transforms #
##############


def ft(dic,data,auto=False,real=False,inv=False,alt=False,neg=False,
    null=False,bruk=False):
    """
    fourier transform

    Differences from NMRPipe:
    - The dmx, and nodmx parameters are not implemented.

    """
    
    size = data.shape[-1]

    # super-flags
    if auto:

        # turn off all flags
        real = False
        inv  = False
        alt  = False
        neg  = False
        null = False
        bruk = False

        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        if dic[fn+"FTFLAG"] == 1.0:   # freq domain
            inv  = True
        else: # freq domain
           
            # Real, TPPI and Sequential data is real transform
            if dic["FD2DPHASE"] == 0 or dic["FD2DPHASE"] == 1:
                real = True

            # sign and negation in AQSIGN
            if dic[fn+"AQSIGN"] == 1 or dic[fn+"AQSIGN"] == 2:
                alt = True
            
            if dic[fn+"AQSIGN"]==16 or dic[fn+"AQSIGN"]==17 \
            or dic[fn+"AQSIGN"]==18:
                alt = True
                neg = True

        # auto debugging
        #print "real:",real
        #print "inv:",inv
        #print "alt:",alt
        #print "neg:",neg
        #print "null:",null
        #print "bruk:",bruk

    if bruk:
        real = True
        alt = True

    if real:    # keep real data
        data.imag = 0.0
    
    if alt: # sign alternate
        if inv == False:    # inv with alt, alternates the inverse
            data[...,1::2] = data[...,1::2]*-1.

    if neg: # negate the imaginary 
        data.imag = data.imag * -1.

    # update the dictionary
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    dic[fn+"FTFLAG"] = (dic[fn+"FTFLAG"]+1)%2   # troggle FT flag
    if dic[fn+"FTFLAG"]==1:
        dic[fn+"FTSIZE"] = data.shape[-1]

    if null:
        # don't perform FT, just find min/max
        dic = update_minmax(dic,data)
        return dic,data

    if inv: # inverse transform
        data = p.icomplexft(data)
        if alt:
            data[...,1::2] = data[...,1::2]*-1
    else:
        data = p.complexft(data)

    if real:
        data = data[...,size/2:]
        dic[fn+"APOD"] = dic[fn+"APOD"]/2.0
        dic[fn+"TDSIZE"] = dic[fn+"TDSIZE"]/2.0
        dic["FDSIZE"] = dic["FDSIZE"]/2.0


    dic = update_minmax(dic,data)

    return dic,data


def rft(dic,data,inv=False):
    """
    real fourier transform

    """

   
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    fn2 = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc

    # if dim is 2, direct dim is real and indirect is complex
    if data.ndim==2 and dic[fn+"QUADFLAG"] == 1 and dic[fn2+"QUADFLAG"] == 0:
        data = data[::2]

    if inv:
        data = p.irft(data.real)
    else:
        data = p.rft(data.real)

    # update the dictionary
    dic[fn+"FTFLAG"]   = (dic[fn+"FTFLAG"]+1)%2   # troggle FT flag
    dic[fn+"QUADFLAG"] = 1.0    # real data
    dic["FDQUADFLAG"]  = 1.0    # real data
    dic = update_minmax(dic,data)
        
    return dic,data


def ha(dic,data,inv=False):
    """
    hadamard transform

    This function is very slow. Implementing a FWHT will significantly improve
    speed.
    """

    data = p.ha(data)

    if inv:
        data = data/data.shape[-1]

    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    dic[fn+"FTFLAG"] = (dic[fn+"FTFLAG"]+1)%2   # troggle FT flag

    # calculation for dictionary updates
    s = data.shape[-1]
    s2 = s/2.0 + 1
    sw = dic[fn+"SW"]
    car = dic[fn+"CAR"]
    obs = dic[fn+"OBS"]
    
    # update the dictionary
    dic[fn+"CENTER"] = s2
    dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs
    dic["FDSIZE"] = s

    return dic,data


def ht(dic,data,mode="ps0-0",zf=False,td=False,auto=False):
    """
    hilbert transform 
    
    Differences from NMRPipe:
    - The ps90-180 mirror image HT mode does not work like NMRPipe's.

    """

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    if auto:
        # when no cutting and 180 P1 correction
        if dic[fn+"P1"]==180.0 and dic[fn+"X1"]==0.0 and dic[fn+"XN"]==0.0:
            zf = False
            mode = "ps90-180"
        else:
            zf = True
            mode = "ps0-0"

    dt = data.dtype

    if mode not in ["ps0-0","ps90-180"]:
        raise ValueError("mode must be ps0-0 or ps90-180")

    if mode == "ps90-180":
        # XXX this gets close but not quite right
        data = -data[::-1]

    if zf:
        N = 2**(np.ceil(np.log2(data.shape[-1]))) #not same as NMRPipe
        fac = N/data.shape[-1]
    else:
        N = data.shape[-1]
        fac = 1.0
    
    z = np.zeros(data.shape,dtype="complex64")

    if data.ndim == 1:
        z[:] = scipy.signal.hilbert(data.real,N)[:data.shape[-1]]*fac
    else:
        for i,vec in enumerate(data):
            z[i] = scipy.signal.hilbert(vec.real,N)[:data.shape[-1]]*fac
    
    # sometimes the hilbert changes the real data, we don't want this
    z.real = data.real

    dic = update_minmax(dic,data)

    # set the QUADFLAG as complex
    dic[fn+"QUADFLAG"] = 0.0
    if fn == "FDF2":
        dic["FDQUADFLAG"] = 0.0

    if td:
        dic[fn+"APOD"] = data.shape[-1]/2.
        #dic[fn+"APOD"] = dic[fn+"APOD"]/2.

    return dic,z

_ht = ht    # macro so ps can call function

####################    
# Simple functions #
####################


def null(dic,data):
    """
    No change
    """

    dic = update_minmax(dic,data)
    return dic,data


def add(dic,data,r=0.0,i=0.0,c=0.0,ri=False,x1=1.0,xn='default'):
    """ add constant to data 
   
    Differences from NMRPipe:
    - Parameter c is used even when r and i are defined.

    """

    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if ri:
        data[...,mn:mx].real = data[...,mn:mx].real+data[...,mn:mx].imag

    else:
        v = r + i*1.j + c + c*1.j
        data[...,mn:mx] = data[...,mn:mx]+v
    dic = update_minmax(dic,data)

    return dic,data

def mult(dic,data,r=1.0,i=1.0,c=1.0,inv=False,hdr=False,x1=1.0,xn='default'):
    """
    multiple data by constant

    Differences from NMRPipe:
    - Parameter c is used even when r and i are defined.

    """

    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if hdr: # read in C from header
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        c = dic[fn+"C1"] 
        r = 1.0
        i = 1.0

    rf = (r*c)  # real factor
    cf = (i*c)  # complex factor

    if inv:
        rf = 1/rf
        cf = 1/cf

    data[...,mn:mx].real = data[...,mn:mx].real*rf
    if np.iscomplex(data).any():
        data[...,mn:mx].imag = data[...,mn:mx].imag*cf

    dic = update_minmax(dic,data)

    return dic,data


def set(dic,data,r="a",i="a",c="a",x1=1.0,xn='default'):
    """
    set data to constant
    
    """

    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if r == "a" and i =="a" and c=="a":
        rc = 0
        ic = 0

    if c !="a":
        rc = c
        ic = c

    if r !="a":
        rc = r

    if i !="a":
        ic = i

    data[...,mn:mx].real = rc
    if np.iscomplex(data).any():
        data[...,mn:mx].imag = ic

    dic = update_minmax(dic,data)
    return dic,data


def rev(dic,data,sw=True):
    """
    reverse data
    """
    
    data = data[...,::-1]
    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    if sw and dic[fn+"FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]-1
        
        # recalc orig
        s = data.shape[-1]
        s2 = dic[fn+"CENTER"]
        sw = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
        dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs

    return dic,data


def dx(dic,data):
    """
    Derivative
    """

    # there is likely a way to do this with np.append 
    z = np.zeros(data.shape,dtype=data.dtype)

    z[...,0]  = data[...,1] - data[...,0]    # first point
    z[...,-1] = data[...,-1] - data[...,-2]  # last point 
    z[...,1:-1]  = data[...,2:] - data[...,:-2] # interior

    dic = update_minmax(dic,z)

    return dic,z


def integ(dic,data):
    """
    Integral (cummulative sum along axis)
    """

    data = np.cumsum(data,axis=-1)
    dic = update_minmax(dic,data)
    return dic,data

def mc(dic,data,mode="mod"):
    """ Modules/Magnitude calculation

    mode should be 'mod' or 'pow' 
    """

    if mode=="mod":
        data = np.sqrt(data.real**2+data.imag**2)
        dic["FDMCFLAG"] = 1.0
    elif mode=="pow":
        data = data.real**2+data.imag**2
        dic["FDMCFLAG"] = 2.0
    else:
        raise ValueError("mode must mod or pow")
    dic = update_minmax(dic,data)

    # change to mag. flags
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    dic[fn+"QUADFLAG"] = 1.0
    dic["FDQUADFLAG"] = 1.0

    return dic,data


def qart(dic,data,a=0.0,f=0.0,auto=False):
    """
    Scale Quad Artifacts

    R' = R
    I' = (1+a)*I + f*R

    Difference from NMRPipe:
    - auto=True performs a Gram-Schmidt orthogonalization from the real and 
      imaginary channels. Therefore none of the grid search parameters are 
      implemented.
    
    """

    if auto:
        # use Gram-Schmidt coefficents to remove quad artifacts
        A,B = calc_gram_schmidth(data)
        a = A-1
        f = B

    #n = np.empty(data.shape,dtype=data.dtype)
    #n.real = data.real
    #n.imag = (1+a)*data.imag + f*data.real
    #data = n
    data.imag = (1+a)*data.imag + f*data.real
    dic = update_minmax(dic,data)

    return dic,data


def calc_gram_schmidt(data):
    """
    calculate Gram-Schmidt orthogonalization parameters

    Calculates parameters A,B to generate orthogonal real and imag
    channel data using:

    R' = R
    I' = A*I + B*R
    """

    # similar to method in Hock and Stern "NMR Data Processing" p.61
    
    # sum of correlation between data.real and data.imag
    C = (data.real*data.imag).sum()

    # total power in real channel
    R = (data.real*data.real).sum()

    # remove correlation from imag channel
    idata = data.imag-(C/R)*data.real

    # total power in uncorrelated imag channel
    S = (idata*idata).sum()

    # imag(data'') = R/S*imag(data')
    # imag(data')  = imag(data)-C/R * real(data)
    # therefore:
    # imag(data'') = R/S*imag(data) - R*C/(S*R) * real(data) 
    # so A = R/S, B=-C/(S)

    return( R/S,-C/S)


def cbf(dic,data,last=10,reg=False,slice=slice(None)):
    """ 
    Constant Baseline correction

    Difference from NMRPipe:
    - Parameters reg and slice should be slice objects if explicit correction
      is desired (recall python arrays start at 0 not 1 so the first limit 
      should be ix1-1). 
    - noseq and nodmx are not implemented.

    """

    n = data.shape[-1]*last/100. +1
    correction = data[...,-n:].sum(axis=-1)/n
    
    if reg != False:
        n = len(range(data.shape[-1])[reg])
        correction = data[...,reg].sum(axis=-1)/n
    
    if data.ndim == 2:
        correction = np.array([correction]).transpose()
        data[slice] = data[slice] - correction[slice]
    else: 
        data= data - correction

    dic = update_minmax(dic,data)
    return dic,data


#########################
# Utility Functions     #
# slices, shuffles, etc #
#########################


def ext(dic,data,x1="default",xn="default",y1="default",yn="default",round=1,
    time=False,left=False,right=False,mid=False,pow2=False,sw=True):
    """
    extract region


    Difference from NMRPipe.
    - The sw parameter is on by default.
    - The time parameter is not implemented (poor documentation). 
    - Certain parameters take precident over others: mid > right > left
      and pow > round.

    """

    # set up axis Hz arrays for later...
    #fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc    
    #fn2 = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc   
    #
    #x_np = (dic[fn+"CENTER"]-1)*2  # number of points in X
    #x_sw = float(dic[fn+"ORIG"])     # SW in X dim
    #x_orig = dic[fn+"ORIG"]         # right edge Hz
    #x_hz = np.linspace(x_orig+x_sw*(x_np-1.)/x_np,x_orig,x_np)

    #if data.ndim == 2:
    #    fn2 = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    #    y_np = (dic[fn2+"CENTER"]-1)*2   # number of points in Y
    #    y_sw = float(dic[fn2+"ORIG"])     # SW in Y dim
    #    y_orig = dic[fn2+"ORIG"]       # right Edge Hz
    #    y_hz = np.linspace(y_orig+y_sw*(y_np-1.)/y_np,y_orig,x_np)

    # store old sizes
    old_x = float(data.shape[-1])
    if data.ndim == 2:
        old_y = float(data.shape[0])

    # slice find limits
    if x1 == "default":
        x_min = 0
    else:
        x_min = np.round(x1)-1
        
    if xn == "default":
        x_max = data.shape[-1]
    else:
        x_max = np.round(xn)

    if left:
        x_min = 0
        x_max = int(data.shape[-1]/2)
    
    if right:
        x_min = int(data.shape[-1]/2)
        x_max = data.shape[-1]

    if mid:
        x_min = int(data.shape[-1]/4)
        x_max = int(3*data.shape[-1]/4)

    r_x = round

    if pow2 and (x1 != "default" or xn != "default"):
        r_x = 2**np.ceil(np.log2(x_max-x_min)) 

    # round size to be multiple of r_x when axis is cut
    if x1 != "default" or xn != "default":
        remain_x =(x_min - x_max) % r_x     # -len_x%r_x
        x_min = x_min - np.floor(remain_x/2)
        x_max = x_max + remain_x - np.floor(remain_x/2)

    if x_min < 0:
        x_max = x_max + -x_min
        x_min = 0.0

    if x_max > data.shape[-1]:
        x_min = x_min - (x_max-data.shape[-1])
        x_max = data.shape[-1]


    if data.ndim == 2:  # 2D array so we also have to modify y
        
        if y1 == "default":
            y_min = 0
        else:
            y_min = np.round(y1)-1

        if yn == "default":
            y_max = data.shape[0]
        else:
            y_max = np.round(yn)

        r_y = round

        if pow2:
            r_y = 2**np.ceil(np.log2(y_max-y_min))

        # round only when axis is cut
        if y1 != "default" or yn !="default":
            remain_y = (y_min - y_max) % r_y
            y_min = y_min - np.floor(remain_y/2)
            y_max = y_max + remain_y - np.floor(remain_y/2)

        if y_min < 0:
            y_max = y_max + -y_min
            y_min = 0.0

        if y_max > data.shape[0]:
            y_min = y_min - (y_max-data.shape[0])
            y_min = data.shape[0]

        #print "ymin:",y_min,"ymax:",y_max
        #print "xmin:",x_min,"xmax:",x_max

        data = data[y_min:y_max,x_min:x_max]
        dic["FDSLICECOUNT"] = y_max - y_min
        dic["FDSPECNUM"] = dic["FDSLICECOUNT"]
        dic["FDSIZE"] = x_max - x_min
        
    else:       # 1D Array
        data = data[x_min:x_max]
        dic["FDSIZE"] = x_max - x_min

    # adjust sweep width and ppm calibration
    if sw:
        
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        
        s = data.shape[-1]
       

        if dic[fn+"FTFLAG"] == 0:   # time domain
            dic[fn+"CENTER"] = s/2+1
            dic[fn+"APOD"]   = s
            dic[fn+"TDSIZE"] = s
            # ORIG calculation
            sw  = dic[fn+"SW"] 
            car = dic[fn+"CAR"]
            obs = dic[fn+"OBS"]
            dic[fn+"ORIG"] = -sw * (np.ceil(s/2.)-1)/s + car*obs

        else:   # freq domain
            dic[fn+"X1"] = x_min+1
            dic[fn+"XN"] = x_max
            dic[fn+"APOD"] = np.floor(dic[fn+"APOD"] * s/old_x)

            dic[fn+"CENTER"] = dic[fn+"CENTER"] - x_min
            dic[fn+"SW"] = dic[fn+"SW"] * s/old_x
       
            # XXXORIG calculation (correct way to calc orig)
            sw  = dic[fn+"SW"]
            car = dic[fn+"CAR"]
            obs = dic[fn+"OBS"]
            s2  = dic[fn+"CENTER"]
            dic[fn+"ORIG"] = car*obs-sw*((s-s2)/s)


        if data.ndim == 2:
            
            fn = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc
            s = data.shape[0]
            if dic[fn+"QUADFLAG"] == 0:
                s = s/2

            if dic[fn+"FTFLAG"] == 0: # time domain
                dic[fn+"CENTER"] = s/2+1
                dic[fn+"APOD"]   = s
                dic[fn+"TDSIZE"] = s
                # ORIG calculation
                sw  = dic[fn+"SW"]
                car = dic[fn+"CAR"]
                obs = dic[fn+"OBS"]
                dic[fn+"ORIG"] = -sw * (np.ceil(s/2.)-1)/s + car*obs

            else:   # freq domain
                if y_min != 0:
                    dic[fn+"X1"] = y_min+1
                if y_max != data.shape[0]:
                    dic[fn+"XN"] = y_max
                if y_min != 0 or y_max != data.shape[0]:
                    dic[fn+"APOD"] = np.floor(dic[fn+"APOD"] * s/old_y)
                dic[fn+"CENTER"] = dic[fn+"CENTER"] - y_min
                dic[fn+"SW"] = dic[fn+"SW"] * s/old_y

                # XXXORIG calculation (correct way to calc orig)
                sw  = dic[fn+"SW"]
                car = dic[fn+"CAR"]
                obs = dic[fn+"OBS"]
                s2  = dic[fn+"CENTER"]
                dic[fn+"ORIG"] = car*obs-sw*((s-s2)/s)

            """
            # ORIG calculation
            sw  = dic[fn+"SW"]
            car = dic[fn+"CAR"]
            obs = dic[fn+"OBS"]
            dic[fn+"ORIG"] = -sw * (np.ceil(s/2.)-1)/s + car*obs

            dic[fn+"CENTER"] = s/2+1
            if dic[fn+"FTFLAG"] == 0:   # time domain
                dic[fn+"APOD"]   = s
                dic[fn+"TDSIZE"] = s
            else:   # freq domain
                dic[fn+"FTSIZE"] = s
            """

    dic = update_minmax(dic,data)

    return dic,data


def sign(dic,data,ri=False,r=False,i=False,left=False,right=False,alt=False,
    abs=False,sign=False):
    """
    Sign manipulation utils

    Difference from NMRPipe:
    - All sign manipulation set to true are applied in order appearing in 
      function definition, no parameters take precidence.
    """

    if ri:
        data = -data

    if r:
        data.real = -data.real

    if i:
        data.imag = -data.imag

    if left:
        data[...,:data.shape[-1]/2.] = -data[...,:data.shape[-1]/2.]
 
    if right: 
        data[...,data.shape[-1]/2.:] = -data[...,data.shape[-1]/2.:]

    if alt:
        data[...,1::2] = -data[...,1::2]

    if abs:
        data.real = np.abs(data.real)
        data.imag = np.abs(data.imag)

    if sign:
        data.real = np.sign(data.real)
        data.imag = np.sign(data.imag)

    dic = update_minmax(dic,data)
    return dic,data


def mir(dic,data,mode="left",invl=False,invr=False,sw=True):
    """
    Append mirror image

    Valid modes are: left, right, center, ps90-180, ps0-0
    invr only applied when mode=left and invl when mode=right

    Differences from NMRPipe:
    - Parameter sw set by default.
    - invr and invl are always applied

    """

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if mode not in ['left','right','center','ps90-180','ps0-0']:
        raise ValueError("invalid mode") 

    if dic[fn+"FTFLAG"] == 0: # time domain

        if mode=="left":
            data = np.append(data,data[...,::-1],axis=-1)

        if mode=="right":
            data = np.append(data[...,::-1],data,axis=-1)

        if mode=="center":
            s = data.shape[-1]
            data=np.concatenate( (data[...,s/2:],data,data[...,:s/2]),axis=-1)

        if mode=="ps90-180":
            s = data.shape[-1]
            data = np.concatenate((-data[...,s/2:],data,-data[...,:s/2]),\
                                   axis=-1)

        if invr:
            s = data.shape[-1]
            data[...,s/2:] = -data[...,s/2:]

        if invl:
            s = data.shape[-1]
            data[...,:s/2] = -data[...,:s/2]

        dic["FDSIZE"] = dic["FDSIZE"]*2

        if mode=="ps0-0":
            # no idea how this is a "center" mode
            s = int(data.shape[-1])
            data = np.concatenate( (data[...,s-1:0:-1],data),axis=-1)
            if np.iscomplexobj(data):
                data.imag[...,:s-1] = -data.imag[...,:s-1]
            dic["FDSIZE"] = dic["FDSIZE"]-1
    
    else: # freq domain

        old_size = int(dic["FDSIZE"])

        if mode=="left":
            data = np.append(data,data[...,::-1],axis=-1)
            dic[fn+"CENTER"] = old_size + dic[fn+"CENTER"]

        if mode=="right":
            data = np.append(data[...,::-1],data,axis=-1)
            dic[fn+"CENTER"] = dic[fn+"CENTER"]

        if mode=="center":
            s = data.shape[-1]
            data=np.concatenate( (data[...,s/2:],data,data[...,:s/2]),axis=-1)
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size/2.

        if mode=="ps90-180":
            s = data.shape[-1]
            data = np.concatenate((-data[...,s/2:],data,-data[...,:s/2]),\
            axis=-1)
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size/2.

        if mode=="ps0-0":
            s = int(data.shape[-1])
            data = np.concatenate( (data[...,s-1:0:-1],data),axis=-1)
            if np.iscomplexobj(data):
                data.imag[...,:s-1] = -data.imag[...,:s-1]
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size

        if invr:
            s = data.shape[-1]
            data[...,s/2:] = -data[...,s/2:]

        if invl:
            s = data.shape[-1]
            data[...,:s/2] = -data[...,:s/2]

        # dictionary updates
        dic["FDSIZE"] = data.shape[-1]
        dic[fn+"APOD"] = dic["FDSIZE"]
        dic[fn+"FTSIZE"] = dic["FDSIZE"]
        dic[fn+"TDSIZE"] = dic["FDSIZE"]
        dic[fn+"ZF"] = -dic["FDSIZE"]
        
        s = dic["FDSIZE"]
        dic[fn+"SW"] = dic[fn+"SW"]*float(s)/float(old_size)

        # this is the 'best' way to calculate ORIG
        sw  = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
        center = dic[fn+"CENTER"]
        dic[fn+"ORIG"] = -sw * (s-center)/s + car * obs


    dic = update_minmax(dic,data)
    return dic,data


def coadd(dic,data,cList=[1,1],axis='x',time=False):
    """
    Co-add Data
    """
    # make a empty blank array then add cList[i]*data[...,i:m:k] to it for 
    # each i for 0 to k=len(cList), m is needed to avoid arrays 1 element 
    # too large, when 'y' axis selected use data[i:m:k] 

    if axis not in ['x','y']:
        raise ValueError("axis must be x or y")

    if axis=='x':
        s = list(data.shape)
        s[-1] = int(s[-1]/len(cList))
        n = np.zeros(s,dtype=data.dtype)  # the new array 
        k = len(cList)
        m = s[-1] * k
        for i in range(k):
            n = n + cList[i]*data[...,i:m:k]
        data = n
        dic["FDSIZE"] = data.shape[-1]
 
    if axis=='y':
        s = list(data.shape)
        s[0] = int(s[0]/len(cList))
        n = np.zeros(s,dtype=data.dtype)  # the new array 
        k = len(cList)
        m = s[0] * k
        for i in range(k):
            n = n + cList[i]*data[i:m:k]
        data = n
        dic["FDSLICECOUNT"] = dic["FDSPECNUM"] = data.shape[0]
   

    dic = update_minmax(dic,data)

    idx = ['x','y'].index(axis)

    if time:
        fn = "FDF"+str(int(dic["FDDIMORDER"][idx])) # F1, F2, etc
        dic[fn+"APOD"] = np.floor(dic[fn+"APOD"]/len(cList))
        dic[fn+"TDSIZE"] = np.floor(dic[fn+"TDSIZE"]/len(cList))

    return dic,data


coad = coadd    # macro for coadd


def save(dic,data,name,overwrite=1):
    """
    Save Current Vector

    Note: FDPIPECOUNT is not set by this function since it has no 
    meaning in pipe_proc context
    """

    dic["FDPIPECOUNT"] = 1.0

    if dic["FDDIMCOUNT"] == 1:
        pipe.write_1D(name,data,dic,overwrite)
    else:
        pipe.write_2D(name,data,dic,overwrite)

    dic["FDPIPECOUNT"] = 0.0
    dic = update_minmax(dic,data)
    return dic,data

def shuf(dic,data,mode=None):
    """
    shuffle utilities

    Difference from NMRPipe:
    - mode="rr2ri" ignores the imaginary vector and does NOT create a 
    mis-sized matrix.  In addition true min/max values are recorded.
    - mode="bswap" updates min/max and may result in NaN in the data.
    - mode="r2i" and "i2r" not implemented as pipe does not support integer 
      NMRPipe format.

    """

    valid_modes = ["ri2c","c2ri","ri2rr","rr2ri","exlr","rolr","swap",
        "bswap","r2i","i2r","inv"]

    if mode not in valid_modes:
        raise ValueError("Invalid mode")

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if mode == "ri2c":
        # interleave real and imaginary data
        s = list(data.shape)
        s[-1] = s[-1]*2
        n = np.empty(s,dtype="float32")
        n[...,::2]  = data.real
        n[...,1::2] = data.imag
        data = n

        # update the dictionary
        dic["FDQUADFLAG"] = 1.0
        dic[fn+"QUADFLAG"] = 1.0
        dic[fn+"APOD"] = data.shape[-1]
        dic[fn+"TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]

    if mode == "c2ri":
        # seperate real and imaginary
        s = list(data.shape)
        s[-1] = int(s[-1]/2)
        n = np.empty(s,dtype="complex64")
        n.real = data.real[...,::2]
        n.imag = data.real[...,1::2]
        data = n

        # update the dictionary
        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 0.0
        dic[fn+"APOD"] = data.shape[-1]
        dic[fn+"TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]

    if mode == "ri2rr":
        # appended imaginary data (NMRPipe interleaves for user)
        s = list(data.shape)
        half = int(s[-1])
        s[-1] = half*2
        n = np.empty(s,dtype="float32")
        n[...,:half] = data.real
        n[...,half:] = data.imag
        data = n

        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0] / 2.0
            dic["FDSPECNUM"] = data.shape[0] / 2.0
        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 1.0
        dic["FDSIZE"] = data.shape[-1]

    if mode == "rr2ri":
        # unappend imaginary data (ignores current imag data)
        s = list(data.shape)
        half = int(s[-1] / 2.0)
        s[-1] = half
        n = np.empty(s,dtype="complex64")
        n.real = data[...,:half]
        n.imag = data[...,half:]
        data = n
        
        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0]
            dic["FDSPECNUM"] = data.shape[0]

        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 0.0
        dic["FDSIZE"] = data.shape[-1] 

    if mode == "exlr":
        # exchange left and right
        half = int(data.shape[-1]/2)
        n = np.empty(data.shape,data.dtype)
        n[...,:half] = data[...,half:]
        n[...,half:] = data[...,:half]
        data = n

    if mode == "rolr":
        # rotate left right halves
        half = int(data.shape[-1]/2)
        n = np.empty(data.shape,data.dtype)
        n[...,:half] = data[...,(half-1)::-1]
        n[...,half:] = data[...,:(half-1):-1]
        data = n

    if mode == "swap":
        n = np.empty(data.shape,data.dtype)
        n.real = data.imag
        n.imag = data.real
        data = n

    if mode == "bswap":
        data = data.byteswap()

    if mode == "r2i":
        raise NotImplementedError("Integer Mode not implemented")

    if mode == "i2r":
        raise NotImplementedError("Integer Mode not implemented")

    if mode == "inv":
        # These doesn't seem to do anything....
        #XXX check data with odd number of points
        pass

    # update the dictionary
    dic = update_minmax(dic,data)

    return dic,data


def qmix(dic,data,ic=1,oc=1,cList=[0],time=False):
    """
    complex mixing of input to outputs

    Difference from NMRPipe
    - ic and oc must evenly divide the number of fid (data.shape[0]) and
      cowardly refuses to make invalid length files.

    """
    ic = int(ic)
    oc = int(oc)

    if data.ndim != 2:
        raise ValueError("data must be 2D")

    if data.shape[0] % ic != 0 or data.shape[0] % oc != 0:
        raise ValueError("ic and oc must be divide the number of vectors")

    # transform matrix
    cmatrix = np.array(cList,dtype='float').reshape(ic,oc)
    cmatrix = cmatrix.transpose()

    # create a empty matrix to hold output
    n_cols = data.shape[1] 
    n_rows = data.shape[0]/float(ic)*float(oc)
    n = np.empty((n_rows,n_cols),dtype=data.dtype)

    # remix by 'block'
    for i in range(int(n_rows/oc)):
        block = data[i*ic:(i+1)*ic]
        n[i*oc:(i+1)*oc] = np.dot(cmatrix,block)


    data = n
    dic = update_minmax(dic,data)
    dic["FDSPECNUM"] = data.shape[0]
    dic["FDSLICECOUNT"] = data.shape[0]

    if time:
        fn = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc
        
        dic[fn+"APOD"]  = data.shape[0]
        dic[fn+"TDSIZE"] = data.shape[0]
        if dic[fn+"QUADFLAG"] == 0:
            dic[fn+"APOD"]  = dic[fn+"APOD"]/2.
            dic[fn+"TDSIZE"] = dic[fn+"TDSIZE"]/2.

    return dic,data


def img(dic,data,filter,dx=1.0,dy=1.0,kern=[1],conv=False,thres=None):
    """
    Image processing util

    Difference from NMRPipe:
    - This function wraps when regions extend past the edges, NMRPipe doesn't
      (even though its documentation states it does).  
    - The filter is applied to the imaginary portion of the data.  NMRPipe 
      does something different with the imaginary channel.  
      This means that only data[dy:-dy,dx:-dx].real will match NMRPipe's
      output.

    """

    # deal with thres by making a masked array
    if thres != False:
        if thres==True:
            thres = 0.0 # default value of 0.0
        data = np.ma.masked_less(data,thres)    

    # create an empty copy of the data
    n = np.empty(data.shape,dtype=data.dtype)
    w = "wrap"

    if conv:    # convolution with kernal
        kern = np.array(kern)
        n.real = scipy.ndimage.convolve(data.real,weights=kern,mode=w)
        n.imag = scipy.ndimage.convolve(data.imag,weights=kern,mode=w)
        
        data = n
        dic = update_minmax(dic,data)
        return dic,data

    # else we have a defined filter
    flts = ["median","min","max","amin","amax","range","avg","dev","thresh"]

    if filter not in flts:
        raise ValueError("filter not valid")

    s = (2*dy+1,2*dx+1) # size tuple

    if filter == "median":
        n.real = scipy.ndimage.median_filter(data.real,size=s,mode=w)
        n.imag = scipy.ndimage.median_filter(data.imag,size=s,mode=w)

    if filter == "min":
        n.real = scipy.ndimage.minimum_filter(data.real,size=s,mode=w)
        n.imag = scipy.ndimage.minimum_filter(data.imag,size=s,mode=w)
        data = n
        
    if filter == "max":
        n.real = scipy.ndimage.maximum_filter(data.real,size=s,mode=w)
        n.imag = scipy.ndimage.maximum_filter(data.imag,size=s,mode=w)
        data = n

    # these functions are much slower (rewrite like _nd_image.so for speed)
    if filter == "amin":
        n.real = scipy.ndimage.generic_filter(data.real,i_amin,size=s,mode=w)
        n.imag = scipy.ndimage.generic_filter(data.imag,i_amin,size=s,mode=w)
        data = n

    if filter == "amax":
        n.real = scipy.ndimage.generic_filter(data.real,i_amax,size=s,mode=w)
        n.imag = scipy.ndimage.generic_filter(data.imag,i_amax,size=s,mode=w)
        data = n

    if filter == "range":
        n.real = scipy.ndimage.generic_filter(data.real,i_range,size=s,mode=w)
        n.imag = scipy.ndimage.generic_filter(data.imag,i_range,size=s,mode=w)
        data = n

    if filter == "avg":
        n.real = scipy.ndimage.generic_filter(data.real,i_avg,size=s,mode=w)
        n.imag = scipy.ndimage.generic_filter(data.imag,i_avg,size=s,mode=w)
        data = n

    if filter == "dev":
        n.real = scipy.ndimage.generic_filter(data.real,i_dev,size=s,mode=w)
        n.imag = scipy.ndimage.generic_filter(data.imag,i_dev,size=s,mode=w)
        data = n

    dic = update_minmax(dic,data)

    return dic,data

# Image filter functions

def i_amin(arr):
    """ find minimum absolute value"""
    return  arr[np.abs(arr).argmin()]

def i_amax(arr):
    """ find maximum absolute value"""
    return arr[np.abs(arr).argmax()]

def i_range(arr):
    return arr.max() - arr.min()

def i_avg(arr):
    return arr.avg()

def i_dev(arr):
    return arr.std()


##################
# Misc Functions #
##################


def ps(dic,data,p0=0.0,p1=0.0,inv=False,hdr=False,noup=False,ht=False,
       zf=False,exp=False,tc=0.0):
    """
    phase shift

    Difference from NMRPipe:
    - inv=True will invert an exponential phase correction, NMRPipe doesn't.
    - FDFNP0 and FDFNP1 will be updated unless noup=True.  NMRPipe never seems
      to write to these for 2D+.  
    - Time-domain phase correction with rs, and ls are not implemented in this
      function, rather use the rs or ls functions.
    
    """

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if ht:  # Hilbert transform
      
        dic,data = _ht(dic,data,zf=zf)

    if hdr: # read from header
        p0 = dic[fn+"P0"]
        p1 = dic[fn+"P1"]
    
    if exp:
        data = p.ps_exp(data,p0=p0,tc=tc,inv=inv)
    else:
        data = p.ps(data,p0=p0,p1=p1,inv=inv)
 
    if noup == False:
        dic[fn+"P0"] = p0
        dic[fn+"P1"] = p1


    dic = update_minmax(dic,data)

    return dic,data


def zf(dic,data,zf=1,pad="auto",size="auto",
    mid=False,inter=False,auto=False,inv=False):
    """
    zero fill

    Only one of zf, pad, size should be set....

    Notes: pad over rides zf, size over rides zf and pad
           pad is ignored if inter flag set
    """
   
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc


    if inv: # recover original time domain points
        
        # calculation for dictionary updates
        s = dic[fn+"TDSIZE"]
        s2 = s/2.0 + 1
        sw = dic[fn+"SW"]
        car = dic[fn+"CAR"]
        obs = dic[fn+"OBS"]
    
        # update the dictionary
        dic[fn+"ZF"] = -1.*s
        dic[fn+"CENTER"] = s2
        dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs
        dic["FDSIZE"] = s

        return dic,data[...,:s]

    if inter:   # zero filling between points done first
        data = p.zf_inter(data,zf)
        dic[fn+"SW"] = dic[fn+"SW"]*(zf+1)
        zf = 0
        pad = 0 # NMRPipe ignores pad after a inter zf

    # set zpad, the number of zeros to be padded
    zpad = data.shape[-1]*2**zf-data.shape[-1]

    if pad != "auto":
        zpad = pad

    if size != "auto":
        zpad = size - data.shape[-1]

    # auto is applied on top of over parameters:
    if auto:
        fsize = data.shape[-1]+zpad
        fsize = 2**(np.ceil(np.log(fsize)/np.log(2)))
        zpad = fsize - data.shape[-1]

    if zpad < 0:
        zpad = 0

    data = p.zf(data,pad=zpad,mid=mid)

    # calculation for dictionary updates
    s = data.shape[-1]
    s2 = s/2.0 + 1
    sw = dic[fn+"SW"]
    car = dic[fn+"CAR"]
    obs = dic[fn+"OBS"]
    
    # update the dictionary
    dic[fn+"ZF"] = -1.*s
    dic[fn+"CENTER"] = s2
    dic[fn+"ORIG"] =  -1.*sw*((s-s2)/s)+car*obs

    dic["FDSIZE"] = s

    dic = update_minmax(dic,data)

    return dic,data


def tp(dic,data,hyper=False,nohyper=False,auto=False,nohdr=False):
    """ 
    transpose data

    XXX test if works with TPPI
    """

    dt = data.dtype
    
    if nohyper:
        hyper = False


    if auto:
        if (dic["FD2DPHASE"] == 1) or (dic["FD2DPHASE"] == 2):
            hyper = True
        else:
            hyper = False

    if hyper:   # Hypercomplex transpose
        data = p.pack_complex(p.unpack_complex(data).transpose())
        data = np.array(data,dtype=dt) 
    else:
        data = data.transpose()


    # update the dimentionality and order
    dic["FDSLICECOUNT"],dic["FDSIZE"] = data.shape[0],data.shape[1]
    dic["FDSPECNUM"] = dic["FDSLICECOUNT"]
    
    dic["FDDIMORDER1"],dic["FDDIMORDER2"]=dic["FDDIMORDER2"],dic["FDDIMORDER1"] 
    
    dic['FDDIMORDER'] = [ dic["FDDIMORDER1"], dic["FDDIMORDER2"], 
                          dic["FDDIMORDER3"], dic["FDDIMORDER4"] ]

    if nohdr != True:
        dic["FDTRANSPOSED"] = (dic["FDTRANSPOSED"]+1)%2

    dic = clean_minmax(dic)

    return dic,data


ytp = tp    # alias for tp
xy2yx = tp  # alias for tp


def smo(dic,data,n=1,center=False):
    """
    Smooth Data
    """

    a = p.smo(data,n=n)
    # NMRPipe doesn't truely smooth the left edge of the vector
    for i in range(n):
        a[...,i] = data[...,0:(n+i)].sum(axis=-1) / (n+1+i)

    if center:
        a = data - a

    dic = update_minmax(dic,a)

    return dic,a


def zd(dic,data,wide=1.0,x0=1.0,slope=0,func=0,g=1):
    """
    zero diagonal band
    """

    if x0 == 0:      # pipe takes x0=0 to be x0=1
        x0 = 1.0

    rows = data.shape[0]
    cols = data.shape[-1]

    if slope==0:    # Auto Mode
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc     
        fn2 = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc    
        sw1 = dic[fn+"SW"]
        sw2 = dic[fn2+"SW"]
        slope = cols*sw1/(rows*sw2) 

    # calculation 
    width = wide*2+1        # width of diagonal band
    c_start = x0+slope-1    # start of center diagonal band

    # maximum row is last row or 
    max_r = int(min(rows,np.floor( (cols-c_start+wide)/slope)+1))

    # window function definitions
    if func==0: # boxcar
        window = np.zeros(width)

    elif func==1: # triangle
        window = np.append(np.linspace(1,0,wide+1),np.linspace(0,1,wide+1)[1:])

    elif func==2: # sinebell
        window = 1-np.sin(np.linspace(0,pi,width))

    elif func==3: # Gaussian
        tln2 = np.sqrt(2*np.log(2))
        window = 1-scipy.signal.gaussian(width,g/tln2)
        # this can be re-written without scipy.signal if wanted
    else:
        raise ValueError("functions must be 0 to 3")

    # apply window to diagonal band row-by-row
    for r in xrange(max_r): # r from 0 to max_r-1

        w_min = 0           # window min
        w_max = len(window) # window max

        c_mid = int(r*slope+(c_start))        # middle of diagonal band
        c_min = c_mid-wide
        c_max = c_mid+wide+1

        if c_min < 0:
            w_min = int(-c_min)
            c_min = 0 
        if c_max > cols:
            w_max = int(w_max-(c_max-cols))
            c_max = cols

        w_size = c_max - c_min
        data[r,c_min:c_max] = data[r,c_min:c_max] * window[w_min:w_max]

    dic = update_minmax(dic,data)

    return dic,data


def dev(dic,data):
    """
    Development
    """
    return dic,data


def di(dic,data):
    """
    deletes imaginary data
    """

    data = data.real
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if dic[fn+"QUADFLAG"] == 0.0:
        dic[fn+"QUADFLAG"] = 1

        if fn == "FDF2":
            dic["FDSPECNUM"] = dic["FDSPECNUM"]/2.0
            dic["FDSLICECOUNT"] = dic["FDSPECNUM"]

        if dic["FDF1QUADFLAG"] == 1 and dic["FDF2QUADFLAG"]:
            dic["FDQUADFLAG"] = 1.0

    return dic,data


# baseline correction


def base(dic,data,nl=None,nw=0,first=False,last=False):
    """ linear baseline correction
    """

    if first:
        nl = [1]+nl
    if last:
        nl.append(data.shape[-1])

    # change node list into 0...spape-1 from 1...shape
    for i in xrange(len(nl)):
        nl[i] = nl[i]-1

    data = proc_bl.base(data,nl,nw)
    dic = update_minmax(dic,data)
    return dic,data


def sol(dic,data,mode="low",fl=16,fs=1,head=0):
    """
    solvent filter

    Differences from NMRPipe:
    - only Low Pass filter mode implemented (no documentation on spline
      and polynomial filters
    - Parameters associated with Spline and Polynomial filters (po, sn, sf,
      and poly) not implemented.
    - Parameter "mir" not implemented
    - Parameters "noseq" and "nodmx" not implemented

    """

    if fs not in [1,2,3]:
        raise ValueError("fs must be 1, 2 or 3")

    if fs == 1:
        sh = "boxcar"
    elif fs == 2:
        sh = "sine"
    elif fs == 3:
        sh = "sine2"

    data[...,head:] = proc_bl.sol(data[...,head:],fl*2+1,sh)
    dic = update_minmax(dic,data)
    return dic,data


def med(dic,data,nw=24,sf=16,sigma=5.0):
    """
    median baseline correction

    Difference from NMRPipe:
    - only applied Friedrichs model-free baseline flatting algorithm
      (Friedrichs JBNMR 1995 5 147-153).  NMRPipe does something else.
    - Additional parameter, sigma, adjusts spread of Gaussian which is 
      convoluted with median baseline to determind final baseline.

    """

    data = proc_bl.med(data,mw=nw,sf=sf,sigma=sigma)
    dic = update_minmax(dic,data)

    return dic,data


def poly(dic,data):
    """
    polynomial baseline correction

    Implementation: ???

    See Callaghan et al, JMR 1984 56 101-109.  Also NMRPipe describes time
    domain version.
    """
    raise NotImplementedError


###################################
# Not Implemented Function Shells #
###################################


# Linear Prediction

def lp(dic,data):
    """
    linear prediction

    Implementation: Medium

    Lots of documentation. Talkbox scikits might have some code to use

    """
    raise NotImplementedError


lpc = lp        # lpc is depreciated


def lp2d(dic,data):
    """
    2D linear prediction

    Function was NOT in original NMRPipe paper.

    Prototype in NMRPipe, not on website
    """
    raise NotImplementedError


# Maximum Entropy

def mem(dic,data):
    """
    maximum entropy reconstruction

    Implementation: Hard but well documented.

    Numerous References for implementation.  
    May be able to use scipy.maxentropy (might be Burg. MEM)
    """
    raise NotImplementedError

# Functions which will not be implemented in pipe_proc due to lack of 
# necessity or documentation.

def ml(dic,data):
    """
    maximum likelihood frequency map

    Function was NOT in original NMRPipe paper.

    This function has almost no documentation, implementation may be hard.

    """
    raise NotImplementedError


def ztp(dic,data):
    """ 3D Matrix transpose 

    There is no need in pipe_proc for this function so it will not be
    implemented.  Rather use the iter3D object from the pipe module.

    """
    raise NotImplementedError


def xyz2zyx(dic,data):
    """ 3D Matrix transpose 

    There is no need in pipe_proc for this function so it will not be
    implemented.  Rather use the iter3D object from the pipe module.

    """
    raise NotImplementedError


def ebs(dic,data):
    """
    EBS Reconstruction

    This function is not well documented in NMRPipe and therefore will not be
    implemented in pipe_proc

    """
    raise NotImplementedError


def ann(dic,data):
    """
    Fourier Analysis by Neural Net

    This function is not well documented in NMRPipe and therefore will not be
    implemented in pipe_proc

    """
    raise NotImplementedError
