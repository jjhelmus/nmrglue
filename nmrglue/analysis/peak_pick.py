"""
Module for peak picking
"""

import numpy as np
import numpy.ma as ma
import scipy

##########################
# Segmentation functions #
##########################

def find_connected2D(data,x,y,thres):
    """
    Find all points connected to node (x,y) above threshold.

    Parameters:
    -----------

    * data  2D array of data
    * x     X (0) axis index to starting node.
    * y     Y (1) axis index to starting node.
    * thres Threshold, below this nodes are considered noise.

    Return: tuple of x and y indexes of valid points

    """
    if data.ndim != 2:
        raise ValueError("data must be 2 dimensional")

    if data[x,y] <= thres:  # check if x,y is above threshold
        return []
    
    # initilize
    right = data.shape[0]   # index of right side of array
    top = data.shape[1]     # index of top of array
    Q = [(x,y)]             # queue
    points = [(x,y)]        # list of connected nodes

    # queue based flood fill algorithm
    while Q:    # loop until Q is empty
        x,y = Q.pop(0)  # remove first element of queue (already in points)
        # check all four directions if above thres and not in points, add to 
        # queue and add to points
        # north
        if y+1 != top and data[x,y+1] > thres and (x,y+1) not in points:
            Q.append((x,y+1))
            points.append((x,y+1))
        # south
        if y != 0 and data[x,y-1] > thres and (x,y-1) not in points:
            Q.append((x,y-1))
            points.append((x,y-1))
        # east
        if x+1 != right and data[x+1,y] > thres and (x+1,y) not in points:
            Q.append((x+1,y))
            points.append((x+1,y))
        # west
        if x != 0 and data[x-1,y] > thres and (x-1,y) not in points:
            Q.append((x-1,y))
            points.append((x-1,y))

    return (np.array(points)[:,0],np.array(points)[:,1])


def find_box_limits(x_points,y_points):
    """ 
    Find the box limits which will outline the provided points

    Parameters:
    -----------

    * x_points  Array of X (0) axis indices.
    * y_points  Array of Y (1) axis indices.

    Returns: x_min,x_max,y_min,y_max

    """
    return x_points.min(),x_points.max(),y_points.min(),y_points.max()


def map_segments2D(data,thres):
    """
    Create a map of connected segments of data above threshold

    Parameter:
    ----------

    * data  2D array of data
    * thres Threshold, below this nodes are considered noise.

    Returns: map

    map is a 2D array of integers indicating segments.  0 represents 
    noise, 1 the first segmented region, 2 the second and on.  The index
    value of these regions can be found using np.where.

    """

    # create the masked data array
    mdata = np.ma.masked_less(data,thres)
    
    # build the map to store segments IDs
    map = np.zeros(shape=data.shape,dtype='u2')
    max_segs = np.iinfo(map.dtype).max 
    mark = 1


    # loop and fill maximum until all points are masked
    while mdata.mask.all() != True:
        x,y = np.unravel_index(mdata.argmax(),mdata.shape)
        mark_map2D(mdata,map,x,y,mark)
        mark = mark+1

    return map



def mark_map2D(mdata,map,x,y,mark):
    """
    Mark connected region on segment map starting at node (x,y)

    This functions should be called from map_segments2D

    Parameters:
    -----------
    * mdata Masked 2D data array
    * map   2D integer array mapping out segments
    * x     Starting node location in 0 axis
    * y     Starting node location in 1 axis 
    * mark  Integer to mark map with.

    Returns nothing but modifies mdata mask and map.

    """
   
    # Algorithm is similar to the 2nd presented in:
    # http://en.wikipedia.org/wiki/Flood_fill

    # index limits of array
    right = mdata.shape[0]
    top = mdata.shape[1]
    
    # If node is masked return (shouldn't happen).
    if mdata.mask[x,y] == True:
        return 

    # create the queue with the node
    Q = [(x,y)]

    # loop until Q is empty
    while Q:
        # remove the first element of the queue
        wx,wy = n = Q.pop(0)

        # if working node is not masked, mark it and mask it 
        if mdata.mask[wx,wy] == False:
            #print "Marking:",wx,wy
            map[wx,wy] = mark
            mdata[wx,wy] = ma.masked
        
        # Check all four directions to see if they are not masked, if
        # so add to queue, mark, and mask them
        # north
        if wy+1 != top and mdata.mask[wx,wy+1] == False:
            #print "Marking and adding",wx,wy+1
            Q.append((wx,wy+1))
            map[wx,wy+1] = mark
            mdata[wx,wy+1] = ma.masked
        # south
        if wy != 0 and mdata.mask[wx,wy-1] == False:
            #print "Marking and adding",wx,wy-1
            Q.append((wx,wy-1))
            map[wx,wy-1] = mark
            mdata[wx,wy-1] = ma.masked
        # east
        if wx+1 != right and mdata.mask[wx+1,wy] == False:
            #print "Marking and adding",wx+1,wy
            Q.append((wx+1,wy))
            map[wx+1,wy] = mark
            mdata[wx+1,wy] = ma.masked
        # west
        if wx != 0 and mdata.mask[wx-1,wy] == False:
            #print "Marking and adding",wx-1,wy
            Q.append((wx-1,wy))
            map[wx-1,wy] = mark
            mdata[wx-1,wy] = ma.masked
    


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


