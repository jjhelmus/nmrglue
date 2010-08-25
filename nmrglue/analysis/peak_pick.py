"""
Module for peak picking
"""

import numpy as np
import scipy

##########################
# Peak picking functions #
##########################

class ndwindower_inside(object):
    """
    Creates a N-dimensional window slice iterator which does not touch the
    edges.

    Example
    -------
    >>> a = np.arange(15).reshape(5,5)
    >>> for s in ndwindower_inside(a.shape,(3,3))
    ...     print a[s]
    [[0 1 2]
     [5 6 7]]

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


# ndwindower is based on numpy ndindex class

class ndwindower(object):
    """ 
    An N-dimentional window slice iterator 

    Parameters:

    * shape : Tuple 
    * window: Tuple


    Example
    -------
    >>> a = np.arange(10).reshape(5,2)
    >>> for s in ndwindower(a.shape,(2,2)
    ...     print a[s]
    [[0 1]
     [2 3]]
    [[1]
     [3]]
    [[2 3]
     [4 5]]
    [[3]
     [5]]
    [[4 5]
     [6 7]]
    [[5]
     [7]]
    [[6 7]
     [8 9]]
    [[7]
     [9]]
    [[8 9]]
    [[9]]

    """

    def __init__(self,shape,wsize):
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(shape)
        self.wsize = wsize
    
    def next(self):
        start = self.ndindex.next()
        stop = np.array(start)+np.array(self.wsize)
        return tuple([slice(x,y) for x,y in zip(start,stop)])
    
    def __iter__(self):
        return self

def peakpick_thres(data,thres,msep):
    """
    Pickpick a spectrum using a threshhold-minimum distance algorithm
    
    Find peaks (local maxima) in a arbitrary dimensional NMR spectra above a 
    set noise threshold with a minimal distance between peaks.  When the
    spectra is small and multiple copies can fit into RAM use the _fast
    version of this functions.

    Parameters:

    * data  N-dimensional 
    * thres Threshold value for minimum peak height
    * msep  N-tuple of minimum peak seperations along each axis
    
    """
    peaks = []  # create an empty list of peaks 


    middle = np.floor((np.array(msep)-1)/2.)    # index of middle of window
    ms = [slice(x,x+1) for x in middle]         # middle slice list

    # loop over the windows
    for idx,s in enumerate( ndwindower( data.shape,msep ) ):
        window = data[s]
        max = window.max()
        #print idx
        if max == window[ms] and max > thres:
            upeak = np.unravel_index(idx,data.shape)
            peaks.append(tuple( (np.array(upeak)+middle).astype('int')))
            #print peaks
    return np.array(peaks)


def peakpick_thres_fast(data,thres,msep):
    """
    Fast version of peakpick_thres functions
    """

    # find local maxima mask
    mx = scipy.ndimage.maximum_filter(data,size=msep,mode='constant')
    max_mask = mx == data
    
    # find threshold mask
    thres_mask = np.ma.greater(data,thres)

    # peaks are bitwise anding of these masks
    return np.transpose(np.nonzero(thres_mask*max_mask))


