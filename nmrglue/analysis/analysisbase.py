"""
analysisbase provides general purpose analysis functions and classes used by
several nmrglue.analysis modules
"""

import numpy as np
pi = np.pi


dimension_names = ['A','Z','Y','X']

# new utilities

def recnames(dnames,ls_classes,Ms):
    """
    Determind the names of a records array provided lineshape classes, 
    
    Parameters:

    dnames      List of dimension names.
    ls_classes  List of lineshape classes.
    Ms          List of lineshape lengths.
    
    Returns list of records array names (strings).
    
    """
    names = []
    for d,M,l in zip(dnames,Ms,ls_classes):
        for p in l.pnames(M):
            names.append(d+'_'+p)
    return names

# utility functions

def find_limits(pts):
    """ 
    Find the limits which outline the provided list of points 
    
    Parameters:
    * pts   List of points [(z0,y0,x0),(z1,y1,x1),...]

    Returns (min,max) ie 
    * min   array of minimum indices array([zmin,ymin,xmin]
    * max   array of maximum indices array([zmin,ymin,xmin]
    
    """
    arr_pts = np.array(pts)
    return np.min(arr_pts,0),np.max(arr_pts,0)

def limits2slice(limits):
    """ 
    Create a set of slice objects given list of min,max limits 
    
    Parameters:
    * limits    Tuple of minimum and maximum indices

    Returns: list of slices which will return points between limits.
    
    
    """
    mins,maxs = limits
    return tuple([slice(i,j+1) for i,j in zip(mins,maxs)])

def slice2limits(slices):
    """
    Create a tuple of min,max limits from a set of slices

    Parameters:
    
    * slices: list of slices

    Returns: Tuple of minimum and maximum indices
    
    """
    mins = [s.start for s in slices]
    maxs = [s.stop-1 for s in slices]
    return mins,maxs

def squish(r,axis):
    """
    Squish array r along axis - sum along all but one axis
    """

    # put axis to be squished as the last axis
    N = int(r.ndim)
    r = r.swapaxes(axis,N-1)

    # sum along leading axis N-1 times
    for i in range(N-1):
        r = r.sum(0)
    return r

def pick2linesh(centers,linewidths,amplitudes):
    """
    Convert peakpick.pick output to linesh.fit_NDregion input

    Parameters:

    * centers       Array of estimated peak locations, shape (n_peaks,ndim).
    * linewidths    Array of estimated peak linewidths, shape (n_peaks,ndim).
    * amplitudes    Array of estimated peak amplitude, shape (n_peaks).


    Return: guesses,amp_guesses

    * guesses       P-length list (P is the number of peaks in region) of
                    N-length lists of tuples where each each tuple is the
                    optimiztion starting parameters for a given peak and
                    dimension lineshape.
    * amp_guesses   P-length list of amplitudes.


    """

    # convert center and linewidths
    guesses = []
    for cs,lws in zip(centers,linewidths):
        guesses.append([(c,lw) for c,lw in zip(cs,lws)])
    
    # convert amplitudes
    amp_guesses = list(amplitudes)

    return guesses,amp_guesses

def linesh2pick(guesses,amp_guesses):
    """
    Convert linesh.fit_NDregion input/output to peakpick.pick output

    Parameters:
    
    * guesses       P-length list (P is the number of peaks in region) of
                    N-length lists of tuples where each each tuple is the
                    optimiztion starting parameters for a given peak and
                    dimension lineshape.
    * amp_guesses   P-length list of amplitudes.

    Return: (centers,linewidths,amplitudes)


    * centers       Array of estimated peak locations, shape (n_peaks,ndim).
    * linewidths    Array of estimated peak linewidths, shape (n_peaks,ndim).
    * amplitudes    Array of estimated peak amplitude, shape (n_peaks).

    """

    # convert guesses to centers and linewidths
    linewidths = np.zeros( ( len(guesses),len(guesses[0]) ),dtype='float')
    centers = np.zeros_like(linewidths)

    for p_num,peak_params in enumerate(guesses):
        for dim_num,dim_params in enumerate(peak_params):
            centers[p_num,dim_num] = dim_params[0]
            linewidths[p_num,dim_num] = dim_params[1]

    # convert amplitudes
    amplitudes  = np.array(amp_guesses)

    return centers,linewidths,amplitudes


# Lineshape classes

# these classes are used to simulate and fit lineshapes
# they should have 6 methods:
# sim(self,M,p)     - Using parameters in p simulate a lineshape of length M.
# nparams(self,M)   - Determind the number of parameters needed for a length M 
#                     lineshape.
# guessp(self,sig)  - Estimate parameters of signal sig, these should be 
#                     parameter which might be used for initial least-squares 
#                     fitting
# pnames(self,M)    - Give names to the parameters of a lineshape of length M.
#
# add_edge(self,p,(min,max)) - take into account region limits at min,max
#                              for parameters or bounds p.
# remove_edge(self,p,(min,max)) - remove effects region limits min,max
#                                 for parameters or bounds.


# the gauss1D gives a well documented example which should be used to create
# new lineshape classes as needed

# lineshape router

def ls_str2class(l):
    """ Convert lineshape string to lineshape class """
    if l == "gauss" or l == "g":
        return gauss1D()
    elif l == "lorentz" or l == "l":
        return lorentz1D()
    elif l == "scale" or l == "s":
        return scale1D()
    elif l == "peak" or l == "p":
        return peak1D()
    else:
        raise ValueError("Unknown lineshape %s",(l))


class gauss1D():
    """
    Gaussian (normal) lineshape class
    
    Parameters (mu,sigma):

    * mu    mean (center of mean)
    * sigma variance (width of distribution)

    """
    # This class has extra documentation to give help to users wishing to
    # create their own lineshape classes

    name = "guassian"

    def sim(self,M,p):
        # simulate a lineshape of length M given parameters in p
        # unpack the 2 parameters in p, mu and sigma
        mu,sigma = p
        s2 = sigma**2   
        # simulate the lineshape
        return np.exp(-(np.arange(M)-mu)**2/(2*s2))/(np.sqrt(2*pi*s2))

    def nparam(self,M):
        # return the number of parameters needed to simulate a M length 
        # gaussian (2: mu and sigma)
        return 2

    def guessp(self,sig):
        # estimate the lineshape parameters, these are rough estimates
        # find the center and full-width at half max
        c,fwhm = center_fwhm(sig)
        # relate mu and sigma to the center and FWHM
        return (c,fwhm/2.3548)

    def pnames(self,M):
        # give names to the two parameters
        return ("mu","sigma")

    def add_edge(self,p,(min,max)):
        # return parameters corrected for region limits (edges)
        # must also be able to correct bounds for edges (so test for None)
        if p[0]==None:
            return p
        return p[0]-min,p[1]

    def remove_edge(self,p,(min,max)):
        # return parameters/bounds 'uncorrected' for region limits (edges)
        if p[0]==None:
            return p
        return p[0]+min,p[1]

class peak1D():
    """
    Peak lineshape class

    This is really a gaussian lineshape which takes a center,fwhm as parameters
    """

    name = "peak"

    def sim(self,M,p):
        mu,fwhm = p
        sigma = fwhm/2.3548
        s2 = sigma**2
        return np.exp(-(np.arange(M)-mu)**2/(2*s2))/(np.sqrt(2*pi*s2))

    def nparam(self,M):
        return 2
    
    def guessp(self,sig):
        c,fwhm=center_fwhm(sig)
        return (c,fwhm)

    def pnames(self,M):
        return ("mu","fwhm")

    def add_edge(self,p,(min,max)):
        if p[0]==None:
            return p
        return p[0]-min,p[1]

    def remove_edge(self,p,(min,max)):
        if p[0]==None:
            return p
        return p[0]+min,p[1]

class lorentz1D():
    """
    Lorentzian lineshape class

    Parameters (x0,g):

    * x0    the center of the lorentzian.
    * gamma the scale parameter .
   
    """

    name = "lorentz"

    def sim(self,M,p):
        x0,gamma = p
        return 1./pi*1./(gamma**2 + (np.arange(M)-x0)**2)

    def nparam(self,M):
        return 2

    def guessp(self,sig):
        c,fwhm = center_fwhm(sig)
        return (c,fwhm/2.)

    def pnames(self,M):
        return("x0","gamma")

    def pcorrect(self,p,(min,max)):
        return p[0]-min,p[1]

    def add_edge(self,p,(min,max)):
        if p[0]==None:
            return p
        return p[0]-min,p[1]

    def remove_edge(self,p,(min,max)):
        if p[0]==None:
            return p
        return p[0]+min,p[1]

class scale1D():
    """
    One dimensional scale class

    Simulates a lineshape with functional form:

    1.0,a0,a1,a2,....

    Where a0, a1, ... are the parameters provided.

    """

    name = "scale"
    
    def sim(self,M,p):
        l = np.empty(M,dtype='float')
        l[0] = 1
        l[1:] = p
        return l

    def nparam(self,M):
        return int(M-1)

    def guessp(self,sig):
        return sig[1:]/sig[0]

    def pnames(self,M):
        return tuple(["a%i"%i for i in range(1,M)])
    
    def add_edge(self,p,(min,max)):
        return p

    def remove_edge(self,p,(min,max)):
        return p

# basic lineshape analysis

def center_fwhm(signal):
    """
    Estimate the center and full-width half max of a signal.
    """

    # negate the signal if it appears to be a negative peak
    if -signal.min() > signal.max():
        signal = -signal

    # the center is the highest point in the signal 
    center = signal.argmax()

    # find the points that bracket the first and last crossing of the
    # half max, then use linear extrapolation to find the location of the
    # half max on either side of the maximum.  The difference between these
    # two values is a good approximation of the full width at half max
    max = signal.max()
    hmax = max/2.

    top_args = np.nonzero(signal > hmax)[0]     # all points above half-max
    l_idx = top_args[0]     # index of left hand side above half-max
    r_idx = top_args[-1]    # index of right hand side above half-max

    # solve hmax = mx+b => x = y-b/m
    # for two points x_0 and x_1 this becomes y = (hmax-x_0)/(x_1-x_0)
    # to this value add the index of x_0 to get the location of the half-max.
    
    # left side
    if l_idx==0:
        left = l_idx    # this is a bad guess but the best we can do
    else:
        x_0,x_1 = signal[l_idx-1],signal[l_idx]
        left = l_idx-1+(hmax-x_0)/(x_1-x_0)

    # right side
    if r_idx == len(signal)-1:
        right = r_idx   # this is a bad guess but the best we can do
    else:
        x_0,x_1 = signal[r_idx],signal[r_idx+1]
        right = r_idx+(hmax-x_0)/(x_1-x_0)

    return center,right-left


def center_fwhm_bymoments(signal):
    """
    Estimate the center and full-width half max of a signal using moments
    """
    
    # calculate the zeroth, first and second moment
    x = np.arange(signal.size)
    m0 = signal.sum()
    m1 = (x*signal).sum()
    m2 = (x**2.*signal).sum()

    # mu (the center) is the first over the zeroth moment
    mu = m1/m0
    # sigma (the variance) is  sqrt( abs(m2)-mu**2)
    sigma = np.sqrt(np.abs(m2)-mu**2)

    return mu,sigma*2.3548



# Windowing classes

class ndwindow(object):
    """ 
    An N-dimentional iterator to slice arrays into windows.

    Given the shape of an array and a window size, an 'ndwindow' instance
    iterators over tuples of slices which slice an the array into wsize 
    sub-arrays.  At each iteration, the index of the center of the sub-array 
    is incremented by one along the last dimension.  Array border are ignored
    so the resulting sub-array can be smaller than wsize.  If wsize contains
    even values the window is off center containing an additional point with 
    lower index.

    Parameters
    
    * size  Size of array to generate tuples of slices from.
    * wsize Size of the area to select from array (sub-array maximum size).
    
    Example

    >>> a = np.arange(12).reshape(3,4)
    >>> for s in ndwindow(a.shape,(3,3))
    ...     print a[s]
    [[0 1]
     [4 5]]
    [[0 1 2]
     [4 5 6]]
    [[1 2 3]
     [5 6 7]]
    [[2 3]
     [6 7]]
    [[0 1]
     [4 5]
     [8 9]]
    [[ 0  1  2]
     [ 4  5  6]
     [ 8  9 10]]
    etc.

    """

    def __init__(self,shape,wsize):
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(shape)
        wsize = np.array(wsize)
        self.sub = np.ceil( (wsize-1.)/2. )
        self.add = wsize-1.-self.sub

    def next(self):
        center = self.ndindex.next()
        start = [max(0,i-j) for i,j in zip(center,self.sub)]
        stop = [i+j+1 for i,j in zip(center,self.add)]
        return tuple([slice(x,y) for x,y in zip(start,stop)])
    
    def __iter__(self):
        return self

class ndwindow_index(object):
    """
    An N-dimensional interator object which returns the index of the window 
    center and a ndwindow slice array.

    Equivalent to:

    for slices,index in zip(np.ndindex(shape),ndwindow(shape,wshape)):
        return (index,slice)

    """
    def __init__(self,shape,wsize):
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(shape)
        wsize = np.array(wsize)
        self.sub = np.ceil( (wsize-1.)/2. )
        self.add = wsize-1.-self.sub

    def next(self):
        center = self.ndindex.next()
        start = [max(0,i-j) for i,j in zip(center,self.sub)]
        stop = [i+j+1 for i,j in zip(center,self.add)]
        return center,tuple([slice(x,y) for x,y in zip(start,stop)])
    
    def __iter__(self):
        return self

class ndwindow_inside(object):
    """
    An N-dimentional iterator to slice arrays into uniform size windows.

    Given the shape of an array and a window size, an 'ndwindow_inside' 
    instance iterators over tuples of slices which slice an the array into 
    uniform size wsize sub-arrays.  At each iteration, the index of the top 
    left of the sub-array is incremented by one along the last dimension utill
    the windows would extend past the array border.  All sub-arrays are
    equal sized (wsize).

    Parameters:

    * size  Size of array to generate tuples of slices from.
    * wsize Size of the area to select from array (widow size).

    Example:
    
    >>> a = np.arange(9).reshape(3,3)
    >>> for s in ndwindow_inside(a.shape,(2,2):
    ...     print a[s]
    [[0 1]
     [4 5]]
    [[1 2]
     [5 6]]
    [[2 3]
     [6 7]]
    [[4 5]
     [8 9]]
    [[ 5  6]
     [ 9 10]]
    [[ 6  7]
     [10 11]]

    """
    def __init__(self,shape,wsize):
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(tuple(np.array(shape)-np.array(wsize)+1))
        self.wsize = wsize
    
    def next(self):
        start = self.ndindex.next()
        stop = np.array(start)+np.array(self.wsize)
        return tuple([slice(x,y) for x,y in zip(start,stop)])
    
    def __iter__(self):
        return self

class ndwindow_inside_index(object):
    """
    An N-dimensional interator object which returns the index of the window 
    top-left and a ndwindow_inside slice array.

    Similar to ndwindow_index but reports top left index of window.

    """
    def __init__(self,shape,wsize):
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(tuple(np.array(shape)-np.array(wsize)+1))
        self.wsize = wsize
    
    def next(self):
        start = self.ndindex.next()
        stop = np.array(start)+np.array(self.wsize)
        return (start,tuple([slice(x,y) for x,y in zip(start,stop)]))
    
    def __iter__(self):
        return self
