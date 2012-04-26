"""
analysisbase provides general purpose analysis functions and classes used by
several nmrglue.analysis modules
"""

import numpy as np
pi = np.pi

# helper functions

def neighbors(pt,shape,structure):
    """
    Generate a list of all neightbors to a point
    
    Parameters:
        
    * pt        index of the point to find neighbors of.
    * shape     shape of the region.
    * structure Structure element that defines connections.

    """
    # set middle of structure to False
    s = np.copy(structure)  # copy structure
    middle = [int(np.floor(i/2.)) for i in s.shape] # find middle of structure
    s.flat[np.ravel_multi_index(middle,s.shape)] = False
    offsets = np.argwhere(s)-middle
    
    # loop over the offset adding all valid points to pts
    pts = []
    for offset in offsets:
        npt = pt-offset
        if valid_pt(npt,shape):
            pts.append(tuple(npt))
    return pts

def valid_pt(pt,shape):
    """ Determind if point is valid in a given shaped array"""
    for i,j in zip(pt,shape):
        if i < 0:   # index is not negative
            return False
        if i>=j:    # index is less than j
            return False
    return True


dimension_names = ['A','Z','Y','X']

# utility functions

def find_limits(pts):
    """ 
    Find the limits which outline the provided list of points 
    
    Parameters:

    * pts   List of points [(z0,y0,x0),(z1,y1,x1),...]

    Returns (min,max):

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
    
    Example::

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
        ...

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

    Example::
    
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
