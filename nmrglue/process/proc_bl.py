"""
NMR baseline processing functions 

proc_bl
=======

Provides:

    1. Functions which filter, smooth, and flatten NMR baselines.

calc_bl_* functions calculate baseline on 1D vector.
All other functions work on 1D or 2D data on the last (-1) axis.

Documentation is available in the docstrings and at http://XXX

Status: 

To Do: poly baseline correction
       test all functions

History:
(jjh) 2009.09.14 med and sol functions
(jjh) 2009.09.08 initial code laydown

"""

import numpy as np
import scipy
import scipy.ndimage


def base(data,nl,nw=0):
    """
    baseline correction (first-order) on nodes (nl).

    Parameters:
    data    1D or 2D data array (ndarray)
    nl      List of baseline nodes
    nw      Node width in pts

    """

    if data.ndim == 1:
        data = data-calc_bl_linear(data,nl,nw)
    else:
        out = np.array(data)
        for i,vec in enumerate(data):
            out[i] = data[i]-calc_bl_linear(vec,nl,nw)
        data = out
    return data


def calc_bl_linear(x,nl,nw=0):
    """ calculate baseline using linear approximation between nodes

    Parameters:
    x   
    nl  List of baseline nodes (pts of only noise)
    nw  +/- points to calculate node value 

    """

    bl = np.zeros_like(x)
    for i in range(len(nl)-1):
        
        # minimum and maximum index
        min = nl[i]
        max = nl[i+1]

        # linspace s1 and s2
        s1 = x[min-nw:min+nw+1].mean()
        s2 = x[max-nw:max+nw+1].mean()
    
        bl[min:max+1] = np.linspace(s1,s2,max-min+1)

    return bl

def med(data,mw=24,sf=16,sigma=5.0):
    """
    median baseline correction

    Algorith described in 
    Friedrichs, M.S. JBNMR 1995 5 147-153.

    Parameters:
    data    1D or 2D data array
    mw      Median Window size in pts.
    sf      Smooth window size in pts.
    sigma   Standard-deviation of Gaussian in convolution

    """

    if data.ndim == 1:
        data = data - calc_bl_med(data,mw,sf,sigma)
    else:
        out = np.array(data)
        for i,vec in enumerate(data):
            out[i] = vec - calc_bl_med(vec,mw,sf,sigma)
        data = out

    return data


def calc_bl_med(x,mw,sf,sigma):
    """
    1D median baseline correction

    Calculates baseline using algorithm from: 
    Friedrichs, M.S. JBNMR 1995 5 147-153

    Parameter:
    x       input 1D ndarray
    mw      Median filter width
    sf      Convolution filter size
    sigma   standard-deviation of Gaussian in convolution

    """

    # create extrema array (non extrema values are masked out)
    mask = x == scipy.ndimage.median_filter(x,size=3)
    mask[0] = False     # first pt always extrema
    mask[-1] = False    # last pt always extrema
    e = np.ma.masked_array(x,mask)

    # fill in the median vector
    half_mw = mw/2
    m = scipy.ndimage.median_filter(e,mw+1,mode="mirror")
    # using the median_filter might give slightly different than
    # algorithm but is MUCH faster

    # convolve with a gaussian
    g = scipy.signal.gaussian(sf,sigma)
    g = g/g.sum()

    return scipy.signal.convolve(m,g,mode='same')


def sol(data,w=16,shape="boxcar",filter=None):
    """ solvent filter

    Method described in Marion et al. JMR 1989 84 425-430

    Parameters:
    data    data (fid)
    w       Width (full width) of convolution shape
    shape   Lowpass shape
    filter  Custom Lowpass shape (1D array)


    """

    if filter!=None:
        s = filter
    else:
        if shape == "boxcar":
            s = scipy.signal.boxcar(w)
        elif shape == "sine":
            s = np.cos(np.pi*np.linspace(-0.5,0.5,w))
        elif shape == "sine2":
            s = np.cos(np.pi*np.linspace(-0.5,0.5,w))**2
        elif shape == "gaussian":
            s = scipy.signal.gaussian(w,w/2.)
        else:
            raise ValueError("invalid shape")

    A = s.sum()

    # convolve with shape 
    if data.ndim == 1:
        return data-scipy.signal.convolve(data,s,mode='same')/A
    else:
        return data-scipy.signal.convolve(data,np.atleast2d(s),mode='same')/A


def poly_td(data):
    """ polynomial time domain solvent subtraction

    From NMRPipe paper( appendix):

    when used with the argument -time, fits all data points to a polynomial,
    which is then subtracted from the original data.  It is intended to fit
    and subtract low-freqency solvent signal in the FID, a procedure that 
    often causes less distortions than time-domain convolution methods.
    By default, a fourth-order polynomials is used.  For speed successive
    averages of regions are usually fit, rather than fitting all of the data.



    Alg:

    1. Calculate average of blocks
    2. Fit these averages to polynomial (block parameters)
    3. Back out "real" polynomial parameters from these block parameters
    4. Subtract off the polynomial from data

    """
    # XXX 
    pass

def poly_fd(data):
    """ polynomial frequency doomain baseline correction

    From NMRPipe paper (appendix):
    
    applies a polynomial baseline correction of the order specified by 
    argument -ord via an automated base-line detection method when used
    with argument -auto.  The defauly is a forth-order polynomial. The 
    automated base-line mode works as follows: a copy of a given vector is 
    divided into a series of adjacent sections, typically eight points wide.
    The average value of each section is subtracted from all points in that 
    section, to generate a 'centered' vector.  The intensities of the entire
    centered vector are sorted, and the standard deviation of the noise is
    estimated under the assumption that a given fraction (typically about
    30%) of the smallest intensities belong to the base-line, and that the 
    noise is normally distributed.  This noise estimate is multiplied by a 
    constant, typically about 1.5, to yield a classification threshold.  
    Then, each section in the centered vector is classified as base line only
    if its standard deviation does not exceed the threshold.  These 
    classifications are used to correct the original vector.

    Alg:

    1. Divide into 'blocks'
    2. Center each block and create centered vector
    3. Calculate intensity (abs) of each centered vector
    4. Sort intensities, lower 30% belong to baseline
    5. Fit base line intensities to Normal distribution, gives estimation
       of standard deviation (std) of noise
    6. Classification threshold set to 1.5*std
    7. Qualify each block in centered vector as baseline only 
    (its std < thres) or not (std > thres)
    8. Fit baseline only points to polynomial and substract off


    """
    # XXX
    pass
