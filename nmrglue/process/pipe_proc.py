"""
NMRPipe like processing functions for use with pipe dic,data pairs.

These functions try to mimic NMRPipe processing functions but many differences
exist between to two implementations.  In particular when using this module:

 * hdr=True over-rides all values in the calling function. 
 * A di flag is not used, rather the di functions should be used to delete
   the imaginary portion of a spectra.
 * x1,xn and other limits must be expressed in points.  spec2pnt or other unit 
   conversion functions should be used before calling the processing function 
   to calculate these values.
 * No functions implement the dmx or nodmx flags.

Additional differences from NMRPipe's functions are documented in the 
individual processing functions.

The following functions have not been implemented and will raise a 
NotImplemented exception:

 * ann      Fourier Analysis by Neural Net
 * ebs      EBS Reconstruction
 * lp       Linear Prediction
 * lpc      Linear Predictions
 * lp2d     2D Linear Prediction
 * mac      Macro Language Interpreter
 * mem      Maximum Entropy
 * ml       Maximum likelyhood frequency
 * poly     Polynomail baseline correction
 * xyz2zyx  3D matrix transpose
 * ztp      3D matrix transpose

"""

# external modules
import numpy as np
import scipy
import scipy.signal

# nmrglue modules
from nmrglue.fileio import pipe,fileiobase
import proc_base as p
import proc_bl

pi = np.pi

###################
# Unit conversion #
###################

class unit_conversion(fileiobase.unit_conversion):
    """ 
    Unit converter class that returns NMRPipe like index values.  Useful
    when calling pipe_proc functions

    """
    # NMRPipe indexes from 1 to MAX instead on 0 to MAX-1
    # we need to modify two method to account for this off by one problem
    def __unit2pnt(self,val,units):
        return fileiobase.unit_conversion.__unit2pnt(self,val,units)+1

    def __pnt2unit(self,val,units):
        return fileiobase.unit_conversion.__pnt2unit(self,val-1,units)

def make_uc(dic,data,dim=-1):
    """ 
    Make a unit conversion object which accepts/returns NMRPipe index values.
    """
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

def recalc_orig(dic,data,fn,axis=-1):
    """
    Recalculate origin for given axis
    """
    # ORIG calculation
    s = float(data.shape[axis])
    
    # This really should check that the axis is not the last...
    if dic[fn+"QUADFLAG"] == 0 and axis!=-1:
        s = int(s/2.)

    sw  = dic[fn+"SW"]
    car = dic[fn+"CAR"]
    obs = dic[fn+"OBS"]
    s2  = float(dic[fn+"CENTER"])
    dic[fn+"ORIG"] = car*obs-sw*((s-s2)/s)

    return dic

def update_minmax(dic,data):
    """
    Update the MAX/MIN dictionary keys
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
    Clean MAX/MIN dictionary keys
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

def apod(dic,data,qName=None,q1=1.0,q2=1.0,q3=1.0,
    c=1.0,start=1,size='default',inv=False,one=False,hdr=False):
    """ 
    Generic Apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * qName Apodization function name: SP, EM, GM, GMB, TM, TRI, JMOD.
    * q1    First parameter for apodization.
    * q2    Second parameter for apodization.
    * q3    Third parameter for apodization.
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True to set points outside of window to 1.
            False leaves points outside of window as is.
    * hdr   True to read apodization parameters from Header.

    """
    # calls the appropiate apodization function
    a_list = ['SP','EM','GM','GMB','TM','TRI','JMOD']

    if hdr:
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        qnum = dic[fn+"APODCODE"]
        qName = ["","SP","EM","GM","TM","","TRI","GMB","JMOD"][qnum]

    # Set apod codes here so that all three parameter are set
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    dic[fn+"APODQ1"] = q1
    dic[fn+"APODQ2"] = q2
    dic[fn+"APODQ3"] = q3

    

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
    Exponential apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * lp    Exponential line broadening in Hz.
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * g1    Inverse exponential width in Hz
    * g2    Gaussian broaden width in Hz
    * g3    Location of Gauss maximum (0.0 to 1.0)
    * c     First point scale value.
    * start Number of points in apodization window. Default is full size.
    * size  Starting location of apodization window.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Modified Gaussian Apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * lb    Exponential term in Hz
    * gb    Gaussian term in Hz
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Exponentially Damped J-Modulation Apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * off   Modulation start in fraction of pi radians.
    * j     J-modulation in Hz
    * lb    Exponential line Broadening in Hz
    * sin   Set True for Sin modulation, off parameter ignored.
    * cos   Set True for Cos modulation, off parameter ignored.
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Sine Bell Apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * off   Sine-bell start location as a fraction of pi radians.
    * end   Sine-bell end location as a fraction of pi radians.
    * pow   Sine-bell exponent.  
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Trapezoid Apodization

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * t1    Points in left side of trapezoid.
    * t2    Points in right side of trapezoid.
    * c     First point scale value.
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Default is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Triangular Apodization

    The right side of the apodization is slightly different than NMRPipe.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * loc   Location of Apex in points. "auto" sets to middle trace.
    * lHi   Left side starting height.
    * rHi   Right side starting height.
    * c     First point scale value. 
    * start Number of points in apodization window.
    * size  Starting location of apodization window. Dafault is full size.
    * inv   True for inverse apodization, False for normal apodization.
    * one   True leaves points outside of window as is, False zeros them.
    * hdr   True to read apodization parameters from Header.

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
    Right Shift and Zero Pad

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * rs    Number of points to right shift.  Negative values are intepreted
            as left shifting.
    * sw    True to update chemical shift calibration parameters.

    """

    if rs < 0:  # negative right shifts are left shifts
        return ls(dic,data,ls=-rs,sw=sw)

    data = p.rs(data,pts=rs)
    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if sw and dic[fn+"FTFLAG"] == 1:    
        # we are in freq domain and must update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]+rs
        dic = recalc_orig(dic,data,fn)
            
    return dic,data

def ls(dic,data,ls=0.0,sw=False):
    """
    Left Shift and Zero Pad

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * ls    Number of points to left shift.  Negative values are intepreted
            as right shifting.
    * sw    True to update chemical shift calibration parameters.

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
            dic = recalc_orig(dic,data,fn)
        else:   # time domain
            dic[fn+"APOD"] = data.shape[-1] - ls
            dic[fn+"TDSIZE"] = data.shape[-1] - ls 

    return dic,data

def cs(dic,data,dir,pts=0.0,neg=False,sw=False):
    """
    Circular Shift

    Syntax is different from NMRPipe.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * dir   Direction to shift spectra, 'rs' or 'ls'.
    * pts   Number of points to shift.
    * neg   Negate shifted points.
    * sw    True to update chemical shift calibration parameters.

    """
    if dir == "ls":
        pts = -pts
    elif dir !="rs":
        raise ValueError("dir must be ls or rs")

    data = p.cs(data,pts,neg=neg)
    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if sw and dic[fn+"FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]+pts
        dic = recalc_orig(dic,data,fn)
    
    return dic,data

def fsh(dic,data,dir,pts,sw=True):
    """
    Frequency Shift

    This function does not perfrom a Hilbert transfrom when data is complex, 
    NMRPipe seems to.  As such the results of the imaginary channel differs
    from NMRPipe. In addition MAX/MIN value are slightly different than those
    from NMRPipe.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * dir   Direction of shift, 'rs' or 'ls'.
    * pts   Number of points to shift.
    * sw    True to update chemical shift calibration parameters.

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
        # channel results, the imaginary channel is a mystery.
        #null,data = _ht(dict(dic),data,zf=True)
        #data = np.array(data,dtype="complex64")

    if dir == "ls":
        pts = -pts

    data = p.fsh(data,pts)

    dic = update_minmax(dic,data)

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    if dic[fn+"FTFLAG"] == 1 and sw: # freq domain
        dic[fn+"CENTER"] = dic[fn+"CENTER"] + pts
        dic = recalc_orig(dic,data,fn)
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
    Complex Fourier Transform

    Choosing multiply conflicting modes produces results different from
    NMRPipe.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.

    Set any of the following to True to choose a mode.
    * auto  Choose mode automatically
    * real  Transform Real-only data
    * inv   Inverse transform
    * alt   Sign Alternate
    * neg   Negate Imaginaties
    * null  Do not apply transform, only adjust headers
    * bruk  Process Redfield sequential data (same as alt=True,real=True)

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
            if dic["FDDIMCOUNT"] >= 2.:
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
    Real Fourier Transform

    Gives slightly different result in certain cases.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * inv   True to perform inverse transform.

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
    Hadamard Transform

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * inv   True to scale by 1/N where N is the length of the trace.

    This function is very slow.  Implemented a FWHT in proc_base would 
    significantly improve the speed of this functions.

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
    dic[fn+"CENTER"] = s2
    dic = recalc_orig(dic,data,fn)
    dic["FDSIZE"] = s

    return dic,data


def ht(dic,data,mode="ps0-0",zf=False,td=False,auto=False):
    """
    Hilbert Transform 

    "ps90-180" mirror image mode gives different results than NMRPipe.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * mode  "ps0-0" or "ps90-180" for mirror image mode
    * zf    Zero fill for speed
    * td    Set time-domain parameter to Size/2
    * auto  Set True to select mode and zf parameter automatically.

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
    if mode not in ["ps0-0","ps90-180"]:
        raise ValueError("mode must be ps0-0 or ps90-180")
    if mode == "ps90-180":
        # XXX this gets close but not quite right
        data = -data[::-1]
    if zf:
        N = 2**(np.ceil(np.log2(data.shape[-1]))) #not same as NMRPipe
    else:
        N = data.shape[-1]

    z = np.array(p.ht(data,N),dtype="complex64")
    dic = update_minmax(dic,data)

    # set the QUADFLAG as complex
    dic[fn+"QUADFLAG"] = 0.0
    if fn == "FDF2":
        dic["FDQUADFLAG"] = 0.0

    if td:
        dic[fn+"APOD"] = data.shape[-1]/2.

    return dic,z

_ht = ht    # macro so ps can call function

##########################
# Standard NMR Functions #
##########################

def di(dic,data):
    """
    Delete Imaginaries
    """

    data = p.di(data)
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if dic[fn+"QUADFLAG"] == 0.0:
        dic[fn+"QUADFLAG"] = 1

        if fn == "FDF2":
            dic["FDSPECNUM"] = dic["FDSPECNUM"]/2.0
            dic["FDSLICECOUNT"] = dic["FDSPECNUM"]

        if dic["FDF1QUADFLAG"] == 1 and dic["FDF2QUADFLAG"]:
            dic["FDQUADFLAG"] = 1.0

    return dic,data

def ps(dic,data,p0=0.0,p1=0.0,inv=False,hdr=False,noup=False,ht=False,
       zf=False,exp=False,tc=0.0):
    """
    Phase Shift

    inv=True will correctly invert an expoenential phase correction. FDFNP0 and
    FDFNP1 are updated unless noup=True.  rs and ls are not implemented, call
    rs or ls first.  

    Parameters:
        
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * p0    Zero order phase in degrees.
    * p1    First order phase in defrees.
    * inv   Set True to perform an inverse phase correction.
    * hdr   Use phase values in header dictionary.
    * noup  Don't update values in header dictionary.
    * ht    Set True to use Hilbert transform to reconstruct imaginaries.
    * zf    Set True to zero fill before Hilbert Transform
    * exp   Set True for Exponential correction
    * tc    Exponential decay constant.

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

def tp(dic,data,hyper=False,nohyper=False,auto=False,nohdr=False):
    """ 
    Transpose Data (2D)

    Parameters:

    * dic       Dictionary of NMRPipe parameters.
    * data      array of spectral data.
    * hyper     Hypercomplex tranpose.
    * nohyper   Supress Hypercomplex transpose.
    * auto      Choose mode automatically.
    * nohdr     Do mark the data transposed in the header dictionary.

    """
    # XXX test if works with TPPI
    dt = data.dtype
    
    if nohyper:
        hyper = False

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    fn2= "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc

    if auto:
        if (dic[fn+"QUADFLAG"]!=1) and (dic[fn2+"QUADFLAG"]!=1):
            hyper = True
        else:
            hyper = False

    if hyper:   # Hypercomplex transpose need type recast
        data = np.array(p.tp_hyper(data),dtype="complex64")
    else:
        data = p.tp(data)
        if dic[fn2+"QUADFLAG"] != 1 and nohyper!=True:
            # unpack complex as needed
            data = np.array(p.c2ri(data),dtype="complex64")
        
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

def zf(dic,data,zf=1,pad="auto",size="auto",
    mid=False,inter=False,auto=False,inv=False):
    """
    Zero Fill

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.

    One (and only one) of the following should be defined:
    * zf    Number of times to double the size.
    * pad   Number of zeros to add.
    * size  Desired final size.

    Set these to True for desired operation (some override other parameters):
    * mid   Zero fill in middle.
    * inter Zero fill between points.
    * auto  Round final size to power of 2.
    * inv   Extract time domain data points.

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
        dic = recalc_orig(dic,data,fn)
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

    data = p.zf_pad(data,pad=zpad,mid=mid)

    # calculation for dictionary updates
    s = data.shape[-1]
    s2 = s/2.0 + 1
    sw = dic[fn+"SW"]
    car = dic[fn+"CAR"]
    obs = dic[fn+"OBS"]
    
    # update the dictionary
    dic[fn+"ZF"] = -1.*s
    dic[fn+"CENTER"] = s2
    dic = recalc_orig(dic,data,fn)
    dic["FDSIZE"] = s

    dic = update_minmax(dic,data)

    return dic,data


######################
# Baseline Functions #
######################

def base(dic,data,nl=None,nw=0,first=False,last=False):
    """
    Linear Baseline Correction

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * nl    List of baseline nodes in points.
    * nw    Node width in points.
    * first Set True to include first point of data in node list.
    * last  Set True to include last point of data in node list.

    """
    if first:
        nl = [1]+nl
    if last:
        nl.append(data.shape[-1])

    # change values in node list to start at 0
    for i in xrange(len(nl)):
        nl[i] = nl[i]-1

    data = proc_bl.base(data,nl,nw)
    dic = update_minmax(dic,data)
    return dic,data

def cbf(dic,data,last=10,reg=False,slice=slice(None)):
    """ 
    Constant Baseline correction

    Parameters ref and slice should be python slice objects if explicit
    correction is desired (recall python arrays start at 0 not 1).  The
    noseq and nodmx parameters are not implemented.
    
    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * last  Percent of trace to use for calculating correction.
    * reg   Slice object describing X-axis region(s) to apply correction to.
    * slice Slice object describing Y-axis region(s) to apply correction to.

    """
    
    if reg != False:
        data = proc_bl.cbf_explicit(data,calc=reg,apply=slice)
    else:
        data = proc_bl.cbf(data,last,slice)

    dic = update_minmax(dic,data)
    return dic,data

def med(dic,data,nw=24,sf=16,sigma=5.0):
    """
    Median Baseline Correction

    This function applies Friendrich's model-free baseline flatting algorithm
    (Friendrichs JBNMR 1995 5 147-153).  NMRPipe applies a different algorithm.
    

    Parameters:
    
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * nw    Median window in points.
    * sf    Smoothing filter in points.
    * sigma Gaussian convolution width.

    """

    data = proc_bl.med(data,mw=nw,sf=sf,sigma=sigma)
    dic = update_minmax(dic,data)

    return dic,data

def sol(dic,data,mode="low",fl=16,fs=1,head=0):
    """
    Solvent Filter

    Only low pass filter implemented. mir, noseq, and nodmx parameters not 
    implemented.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * mode  Filter Mode (must be "low")
    * fl    Filter length in points
    * fs    Lowpass filter shape (1=boxcar, 2=sine , 3=Sine^2)
    * head  Number of points to skip

    Differences from NMRPipe:
    - Parameter "mir" not implemented
    - Parameters "noseq" and "nodmx" not implemented

    """

    if fs not in [1,2,3]:
        raise ValueError("fs must be 1, 2 or 3")

    if fs == 1:
        data[...,head:] = proc_bl.sol_boxcar(data[...,head:],w=fl*2+1)
    elif fs == 2:
        data[...,head:] = proc_bl.sol_sine(data[...,head:],w=fl*2+1)
    elif fs == 3:
        data[...,head:] = proc_bl.sol_sine2(data[...,head:],w=fl*2+1)

    dic = update_minmax(dic,data)
    return dic,data


###################  
# Basic Utilities #
###################

def add(dic,data,r=0.0,i=0.0,c=0.0,ri=False,x1=1.0,xn='default'):
    """ 
    Add a Constant

    Parameter c is used even when r and i are defined.  NMRPipe ignores c when
    r or i are defined.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * r     Constant to add to real data.
    * i     Constant to add to imaginary data.
    * c     Constant to add to read and imaginary data
    * ri    Add real and imaginary data into real channel.
    * x1    First point of region to add constant to.
    * xn    Last point of region to add constant to. 'default' specifies the 
            end of the vector.

    """
    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if ri:
        data[...,mn:mx].real = p.add_ri(data[...,mn:mx])

    else:
        data[...,mn:mx] = p.add(data[...,mn:mx],r,i,c)
    dic = update_minmax(dic,data)

    return dic,data

def dx(dic,data):
    """
    Derivative by central difference.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.

    """
    data = p.dx(data)
    dic = update_minmax(dic,data)
    return dic,data

def ext(dic,data,x1="default",xn="default",y1="default",yn="default",round=1,
    left=False,right=False,mid=False,pow2=False,sw=True):
    """
    Extract Region

    The time parameter is not implemented.  Using multiple conflicting 
    parameters may result in different results than NMRPipe.

    Parameters:
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * x1    X-axis extract region start
    * xn    X-axis extract region stop
    * y1    Y-axis extract region start
    * yn    Y-axis extract region stop
    * round Round extract size to nearest N points
   
    Set any of the following to True to select the desired extraction.
    * left  Extract Left Half.
    * right Extract Right Half.
    * mid   Extract Center Half.

    * pow2  Set True to round extracted size to nearest power of 2.
    * sw    Set True to update Sweep Width and ppm calibration

    """
    # this function does not weap proc_base.ext rather the slice is performed 
    # here

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
        if y_min!=1 and y_max !=data.shape[0]:  # only update when sliced
            dic["FDSLICECOUNT"] = y_max - y_min
        dic["FDSPECNUM"] = y_max - y_min
        dic["FDSIZE"] = x_max - x_min
        
    else:       # 1D Array
        data = data[x_min:x_max]
        dic["FDSIZE"] = x_max - x_min

    # adjust sweep width and ppm calibration
    if sw:
        
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
        
        s = data.shape[-1]
       

        if dic[fn+"FTFLAG"] == 0:   # time domain
            dic[fn+"CENTER"] = float(int(s/2.+1))
            dic[fn+"APOD"]   = s
            dic[fn+"TDSIZE"] = s
            dic = recalc_orig(dic,data,fn)

        else:   # freq domain
            dic[fn+"X1"] = x_min+1
            dic[fn+"XN"] = x_max
            dic[fn+"APOD"] = np.floor(dic[fn+"APOD"] * s/old_x)
            dic[fn+"CENTER"] = dic[fn+"CENTER"] - x_min
            dic[fn+"SW"] = dic[fn+"SW"] * s/old_x
            dic = recalc_orig(dic,data,fn)

        if data.ndim == 2:
            
            fn = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc
            s = data.shape[0]
            if dic[fn+"QUADFLAG"] == 0:
                s = s/2

            if dic[fn+"FTFLAG"] == 0: # time domain
                dic[fn+"CENTER"] = s/2+1
                dic[fn+"APOD"]   = s
                dic[fn+"TDSIZE"] = s
                dic = recalc_orig(dic,data,fn,-2)

            else:   # freq domain
                if y_min != 0:
                    dic[fn+"X1"] = y_min+1
                if y_max != data.shape[0]:
                    dic[fn+"XN"] = y_max
                if y_min != 0 or y_max != data.shape[0]:
                    dic[fn+"APOD"] = np.floor(dic[fn+"APOD"] * s/old_y)
                dic[fn+"CENTER"] = dic[fn+"CENTER"] - y_min
                dic[fn+"SW"] = dic[fn+"SW"] * s/old_y
                dic = recalc_orig(dic,data,fn,-2)

    dic = update_minmax(dic,data)

    return dic,data

def integ(dic,data):
    """
    Integral by Simple Sum

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.

    """
    data = p.integ(data)
    dic = update_minmax(dic,data)
    return dic,data

def mc(dic,data,mode="mod"):
    """ 
    Modules/Magnitude Calculation

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * mode  "mod" or "pow" for Modules or Square Modules

    """
    if mode=="mod":
        data = p.mc(data)
        dic["FDMCFLAG"] = 1.0
    elif mode=="pow":
        data = p.mc_pow(data)
        dic["FDMCFLAG"] = 2.0
    else:
        raise ValueError("mode must mod or pow")
    dic = update_minmax(dic,data)

    # change to mag. flags
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    dic[fn+"QUADFLAG"] = 1.0
    dic["FDQUADFLAG"] = 1.0

    return dic,data

def mir(dic,data,mode="left",invl=False,invr=False,sw=True):
    """
    Append Mirror Image

    Negations selected are applied regardless of mode selected.

    Parameters:
        
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * mode  Type of mirror image to apply valid modes are 'left', 'right',
            'center', 'ps90-180','pw0-0'
    * invl  Set True to negate left half
    * invr  Set True to negate right half
    
    """

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if mode not in ['left','right','center','ps90-180','ps0-0']:
        raise ValueError("invalid mode") 

    if dic[fn+"FTFLAG"] == 0: # time domain
        if mode=="left":    data = p.mir_left(data) 
        if mode=="right":   data = p.mir_right(data)
        if mode=="center":  data = p.mir_center(data)
        if mode=="ps90-180":    data = p.neg_edges(p.mir_center(data))
        if invr:    data = p.neg_right(data)    
        if invl:    data = p.neg_left(data)
        dic["FDSIZE"] = dic["FDSIZE"]*2
        if mode=="ps0-0":
            data = p.mir_center_onepoint(data)
            dic["FDSIZE"] = dic["FDSIZE"]-1
    
    else: # freq domain

        old_size = int(dic["FDSIZE"])

        if mode=="left":
            data = p.mir_left(data)
            dic[fn+"CENTER"] = old_size + dic[fn+"CENTER"]
        if mode=="right":
            data = p.mir_right(data)
            dic[fn+"CENTER"] = dic[fn+"CENTER"]
        if mode=="center":
            data = p.mir_center(data)
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size/2.
        if mode=="ps90-180":
            data = p.neg_edges(p.mir_center(data))
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size/2.
        if mode=="ps0-0":
            data = p.mir_center_onepoint(data)
            dic[fn+"CENTER"] = dic[fn+"CENTER"] + old_size
        if invr:
            data = p.neg_right(data)
        if invl:
            data = p.neg_left(data)

        # dictionary updates
        dic["FDSIZE"] = data.shape[-1]
        dic[fn+"APOD"] = dic["FDSIZE"]
        dic[fn+"FTSIZE"] = dic["FDSIZE"]
        dic[fn+"TDSIZE"] = dic["FDSIZE"]
        dic[fn+"ZF"] = -dic["FDSIZE"] 
        s = dic["FDSIZE"]
        dic[fn+"SW"] = dic[fn+"SW"]*float(s)/float(old_size)
        dic = recalc_orig(dic,data,fn)  # recalculate origin

    dic = update_minmax(dic,data)
    return dic,data

def mult(dic,data,r=1.0,i=1.0,c=1.0,inv=False,hdr=False,x1=1.0,xn='default'):
    """
    Multiple by a Constant

    Parameter c is used even when r and i are defined.  NMRPipe ignores c when
    r or i are defined.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * r     Constant to multply real data by.
    * i     Constant to multiply imaginary data by.
    * c     Constant to multiply both real and imaginary data by.
    * inv   Multiply by inverse of Constant (both real and imaginary)
    * hdr   Use constant value from header.
    * x1    First point of region to multiply constant by.
    * xn    Last point of region to multiply constant by. 'default' specifies 
            the end of the vector.

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

    data[...,mn:mx] = p.mult(data[...,mn:mx],r=rf,i=cf,c=1.0)
    dic = update_minmax(dic,data)
    return dic,data

def rev(dic,data,sw=True):
    """
    Reverse Data

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * sw    Adjust sweep width and ppm calibration.

    """
    data = p.rev(data)
    dic = update_minmax(dic,data)
    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc
    if sw and dic[fn+"FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn+"CENTER"] = dic[fn+"CENTER"]-1
        dic = recalc_orig(dic,data,fn)

    return dic,data

def set(dic,data,r="a",i="a",c="a",x1=1.0,xn='default'):
    """
    Set to a Constant

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * r     Constant to set real data to. "a" sets to 0 unless c is defined.
    * i     Constant to set imag data to. "a" sets to 0 unless c is defined.
    * c     Constant to set real and imaginary. "a" sets to 0 unless
            r or i is defined. 
    * x1    First point of region to set to constant.
    * xn    Last point of region to set to constant. 'default' specifies the 
            end of the vector.
    
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

    # this is so simple we do not use the proc_base functions 
    data[...,mn:mx].real = rc
    if np.iscomplex(data).any():
        data[...,mn:mx].imag = ic

    dic = update_minmax(dic,data)
    return dic,data

def shuf(dic,data,mode=None):
    """
    Shuffle Utilities

    rr2ri mode ignores any imaginary vector refusing to create a mis-sized
    vector.  bswap mode may results in NaN in the data.  r2i and i2r not 
    implemented.  All modes correctly update minimum and maximum values

    Parameters:
        
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * mode  shuffle mode to apply:

    Valid modes are:
    
    * ri2c  Interleave real/imaginary.
    * c2ri  Seperate real/imaginary.
    * ri2rr Append real/imaginary.
    * rr2ri Unappend real/imaginary.
    * exlr  Exchange left/right halfs.
    * rolr  Rotate left/right halfs.
    * swap  Swap real/imaginary.
    * bswap Byte-swap.
    * inv   Do nothing.
    
    """

    valid_modes = ["ri2c","c2ri","ri2rr","rr2ri","exlr","rolr","swap",
        "bswap","r2i","i2r","inv"]

    if mode not in valid_modes:
        raise ValueError("Invalid mode")

    fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc

    if mode == "ri2c":
        data = p.ri2c(data) # interleave real and imaginary data
        # update the dictionary
        dic["FDQUADFLAG"] = 1.0
        dic[fn+"QUADFLAG"] = 1.0
        dic[fn+"APOD"] = data.shape[-1]
        dic[fn+"TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]

    if mode == "c2ri":
        # seperate real and imaginary
        data = np.array(p.c2ri(data),dtype="complex64")
        # update the dictionary
        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 0.0
        dic[fn+"APOD"] = data.shape[-1]
        dic[fn+"TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]

    if mode == "ri2rr":
        data = p.ri2rr(data)    # appended imaginary data
        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0] / 2.0
            dic["FDSPECNUM"] = data.shape[0] / 2.0
        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 1.0
        dic["FDSIZE"] = data.shape[-1]

    if mode == "rr2ri":
        # unappend imaginary data (ignores imag data)
        data = np.array(p.rr2ri(data),dtype="complex64")
        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0]
            dic["FDSPECNUM"] = data.shape[0]

        dic["FDQUADFLAG"] = 0.0
        dic[fn+"QUADFLAG"] = 0.0
        dic["FDSIZE"] = data.shape[-1] 

    if mode == "exlr":
        data = p.exlr(data) # exchange left and right 
    if mode == "rolr":
        data = p.rolr(data) # rotate left right halves
    if mode == "swap":
        data = p.swap(data)
    if mode == "bswap":
        data = p.bswap(data)
    if mode == "r2i":
        raise NotImplementedError("Integer Mode not implemented")
    if mode == "i2r":
        raise NotImplementedError("Integer Mode not implemented")
    if mode == "inv":
        # This does not seem to do anything....
        #XXX check data with odd number of points
        pass
    # update the dictionary
    dic = update_minmax(dic,data)

    return dic,data

def sign(dic,data,ri=False,r=False,i=False,left=False,right=False,alt=False,
    abs=False,sign=False):
    """
    Sign Manipulation Utilities

    All sign manupulation modes set True are applied.

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.

    Select one or more of the desired modes by setting to True:
    * ri    Negate all data.
    * r     Negate real data.
    * i     Negate imaginary data.
    * left  Negate left half.
    * right Negate right half.
    * alt   Negate Alterate points.
    * abs   Replace data with absolute value of data.
    * sign  Replace data with sign (-1 or 1) of data.

    """
    if ri:  data = p.neg_all(data)
    if r:   data = p.neg_real(data)
    if i:   data = p.neg_imag(data)
    if left: data = p.neg_left(data)
    if right: data = p.neg_right(data)
    if alt: data = p.neg_alt(data)
    if abs: data = p.abs(data)
    if sign: data = p.sign(data)
    dic = update_minmax(dic,data)
    return dic,data

##################
# Misc Functions #
##################

def coadd(dic,data,cList=[1,1],axis='x',time=False):
    """
    Co-Addition of Data

    Parameters:
    
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * cList List of coefficients
    * axis  Axis to co-add from/to, either 'x' or 'y'
    * time  Set True to adjust time-domain headers to account for data 
            size reduction

    """
    if axis=='x': 
        data = p.coadd(data,cList,axis=-1)
        dic["FDSIZE"] = data.shape[-1]
        idx = 0
    elif axis=='y': 
        data = p.coadd(data,cList,axis=0)
        dic["FDSLICECOUNT"] = dic["FDSPECNUM"] = data.shape[0]
        idx = 1
    else:
        raise ValueError("axis must be x or y")

    dic = update_minmax(dic,data)
    if time:
        fn = "FDF"+str(int(dic["FDDIMORDER"][idx])) # F1, F2, etc
        dic[fn+"APOD"] = np.floor(dic[fn+"APOD"]/len(cList))
        dic[fn+"TDSIZE"] = np.floor(dic[fn+"TDSIZE"]/len(cList))
    return dic,data

coad = coadd    # macro for coadd

def dev(dic,data):
    """
    Development Function (Null Function)
    """
    return dic,data

def img(dic,data,filter,dx=1.0,dy=1.0,kern=[1],conv=False,thres=None):
    """
    Image Processing Utilities

    This function wraps when regions extend past the edges (NMRPipe doesn't).
    The filter is applied to both the real and imaginary channels

    Parameters:
        
    * dic    Dictionary of NMRPipe parameters.
    * data   array of spectral data.
    * filter Filter to apply
    * dx     Filter X-axis width in points.
    * dy     Filter Y-axis width in points.
    * kern   List of kernal values
    * conv   Set True to apply convolution filter
    * thres  Threshold value for use in computing filter, None turns off.

    Supported filters are:
    
    * median    Median
    * min       Minimum
    * max       Maximim
    * amin      Absolute Minimum
    * amax      Absolute Maximum
    * range     Range
    * avg       Average
    * dev       Standard Deviation

    """

    # deal with thres by making a masked array
    if thres != None:
        if thres==True:
            thres = 0.0 # default value of 0.0
        data = p.thres(data,thres) 

    if conv:    # convolution with kernal
        data = p.conv(data,kern,m="wrap")
        dic = update_minmax(dic,data)
        return dic,data

    s = (2*dy+1,2*dx+1) # size tuple
    # the various filters
    if filter == "median":  data = p.filter_median(data,s=s,m="wrap")
    elif filter == "min":   data = p.filter_min(data,s=s,m="wrap")
    elif filter == "max":   data = p.filter_max(data,s=s,m="wrap")
    elif filter == "amin":  data = p.filter_amin(data,s=s,m="wrap")
    elif filter == "amax":  data = p.filter_amax(data,s=s,m="wrap")
    elif filter == "range": data = p.filter_range(data,s=s,m="wrap")
    elif filter == "avg":   data = p.filter_avg(data,s=s,m="wrap")
    elif filter == "dev":   data = p.filter_dev(data,s=s,m="wrap")
    else:
        raise ValueError("Invalid filter")

    dic = update_minmax(dic,data)
    return dic,data

def null(dic,data):
    """
    Null Function (no change) 
    """
    dic = update_minmax(dic,data)
    return dic,data

def qart(dic,data,a=0.0,f=0.0,auto=False):
    """
    Scale Quad Artifacts

    Auto mode performs Gram-Schmidt orthogonalization.  No grid search is
    performed.

    Parameters:
    
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * a     Amplitude Adjustment
    * f     Phase Adjustment
    * auto  Perform Gram-Schmidth orthogonalization to find a and f.
    """
    if auto:
        data = p.qart_auto(data)    
    else:
        data = p.qart(data,a,f)

    dic = update_minmax(dic,data)
    return dic,data

def qmix(dic,data,ic=1,oc=1,cList=[0],time=False):
    """
    Complex Mixing of Input to Outputs

    ic and oc must evenly divide the matrix shape.  Refuses to make invalid
    length files.
 
    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * ic    Input channel count.
    * oc    Output channel count.
    * cList ic by oc list or array
    * time  Set True to adjust time-domain size headers.

    """
    ic = int(ic)
    oc = int(oc)

    if data.ndim != 2:
        raise ValueError("data must be 2D")

    if data.shape[0] % ic != 0 or data.shape[0] % oc != 0:
        raise ValueError("ic and oc must be divide the number of vectors")

    carr = np.array(cList,dtype='float').reshape(ic,oc)
    data = p.qmix(data,carr)

    #data = n
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

def save(dic,data,name,overwrite=True):
    """
    Save Current Vector

    The resulting FDPIPECOUNT header parameter does not NMRPipe's.

    Parameters:
    
    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * name  Name of file to write to.
    * overwrite Set True to overwrite existing files.
    
    """

    dic["FDPIPECOUNT"] = 1.0

    if dic["FDDIMCOUNT"] == 1:
        pipe.write_1D(name,dic,data,overwrite)
    else:
        pipe.write_2D(name,dic,data,overwrite)

    dic["FDPIPECOUNT"] = 0.0
    dic = update_minmax(dic,data)
    return dic,data

def smo(dic,data,n=1,center=False):
    """
    Smooth Data

    Parameters:

    * dic    Dictionary of NMRPipe parameters.
    * data   array of spectral data.
    * n      Smoothing window in points.
    * center Set True to perform centering (subtract smoothed data)

    """

    a = p.smo(data,n=n)
    # NMRPipe doesn't truely smooth the left edge of the vector
    for i in range(n):
        a[...,i] = data[...,0:(n+i)].sum(axis=-1) / (n+1+i)
    if center:
        # to avoid the same error center without use proc_base functions
        a = data - a

    dic = update_minmax(dic,a)

    return dic,a

def zd(dic,data,wide=1.0,x0=1.0,slope=0,func=0,g=1):
    """
    Zero Diagonal Band

    Parameters:

    * dic   Dictionary of NMRPipe parameters.
    * data  array of spectral data.
    * wide  Diagonal band width in points
    * x0    Diagonal start location in points.
    * slope Diagonal slope (X/Y ratio), 0 for auto mode.
    * func  0=boxcar, 1=Triangle, 2=Sine Bell, 3=Gaussian.
    * g     Gauss Width.

    """
    if x0 == 0:      # pipe takes x0=0 to be x0=1
        x0 = 1.0

    if slope==0:    # Auto Mode
        fn = "FDF"+str(int(dic["FDDIMORDER"][0])) # F1, F2, etc     
        fn2 = "FDF"+str(int(dic["FDDIMORDER"][1])) # F1, F2, etc    
        sw1 = dic[fn+"SW"]
        sw2 = dic[fn2+"SW"]
        slope = data.shape[-1]*sw1/(data.shape[0]*sw2) 

    if func==0:     data = p.zd_boxcar(data,wide,x0-1,slope)
    elif func==1:   data = p.zd_triangle(data,wide,x0-1,slope)
    elif func==2:   data = p.zd_sinebell(data,wide,x0-1,slope)
    elif func==3:   data = p.zd_gaussian(data,wide,x0-1,slope,g)
    else:   
        raise ValueError("func must be 0,1,2 or 3")

    dic = update_minmax(dic,data)
    return dic,data

#############################
# Not Implemented Functions #
#############################

def ann(dic,data):
    """
    Fourier Analysis by Neural Net
    """
    raise NotImplementedError

def ebs(dic,data):
    """
    EBS Reconstruction
    """
    raise NotImplementedError

def lp(dic,data):
    """
    Linear Prediction
    """
    raise NotImplementedError

lpc = lp        # lpc is depreciated

def lp2d(dic,data):
    """
    2D Linear Prediction
    """
    raise NotImplementedError

def mac(dic,data):
    """
    Macro Language Interpreter
    """
    raise NotImplementedError

def mem(dic,data):
    """
    Maximum Entropy Reconstruction
    """
    raise NotImplementedError

def ml(dic,data):
    """
    Maximum Likelihood Frequency Map
    """
    raise NotImplementedError

def poly(dic,data):
    """
    Polynomial Baseline Correction
    """
    raise NotImplementedError

def xyz2zyx(dic,data):
    """ 
    3D Matrix transpose 
    """
    raise NotImplementedError

def ztp(dic,data):
    """ 
    3D Matrix Transpose 
    """
    raise NotImplementedError
