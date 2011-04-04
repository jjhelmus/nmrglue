"""
segmentation - functions to perform segmentation of NMR spectrum.
"""

import numpy as np
import numpy.ma as ma

# helper functions

def neighbors(pt,shape,diag=False):
    """
    Generate a list of all neightbors to a point
    
    Parameters:
        
    * pt        Index of the point to find neighbors of.
    * shape     Shape of the region
    * diag      True to include diagonal neightbors, False not to.

    """
    # developed from scipy.ndimage.morphology.generate_binary_structure
    rank = len(shape)
    pt = np.array(pt)

    if diag:
        conn = rank
    else:
        conn = 1
    
    # create an array of offsets
    o = np.fabs(np.indices([3]*rank) - 1).sum(0)
    offsets = np.argwhere(np.bitwise_and(o<=conn,o!=0))-1

    pts = []
    # loop over the offset adding all valid points to pts
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


# Segmentation functions

# Downward segmentation method:
# The downward segmenting method uses the flood fill algorithm to find
# all points connected to an initial node which are above a given threshold
# and to which a path exists in which each step of the path moves lower in 
# intensity.  This can be though of as all points accessible by a water drop
# following downward slopes from the initial node.

# Upward segmentation uses the same priciple except nodes must be below
# the threshold an upward path must exist.

def find_downward(data,pt,thres,diag=False):
    """
    Find points downward-connected to pt in data.

    Parameters:

    * data  2D array of data
    * pt    Starting point, top of peak, tuple.
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Return: list of indicies of downward connected nodes.

    """
    if type(pt) == int:
        pt = (pt,)
    pt = tuple(pt)
    shape = data.shape

    if data[pt] < thres:    # check that the initial point is above threshold.
        return []

    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node

        for new_pt in neighbors(pt,shape,diag):  # check all neightbors

            if thres<data[new_pt]<v and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    
    return segment

def find_all_downward(data,thres,diag=False):  
    """
    Find all downward-connected segments in data

    Parameters:

    * data  Array of data
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Returns: (centers,segments)

    * centers   List of indicies of local maximum in each segment
    * segments  List of all points in a given segment

    """

    mdata = np.ma.masked_less(data,thres)

    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all()!=True:

        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmax(),data.shape)

        # find all nodes downward connected to pt and mark them
        segment = find_downward(data,pt,thres,diag)

        for i in segment:   # probably a better way of doing this
            mdata[i] = ma.masked
        
        centers.append(pt)
        segments.append(segment)

    return centers,segments


def map_downward(data,thres,diag=False,map_dtype='u2'):
    """
    Create a map of downward-connected segments of data above threshold

    Parameters:

    * data  Array of data
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.
    * map_dtype Dtype of map.

    Returns: map

    map is an array of integers indicating segments. 0 represents noise,1 the
    first segment, 2 the second, etc.  The index value of these regions
    can be found using np.where

    """
    # create the masked data array
    mdata = np.ma.masked_less(data,thres)

    # build the map to store segment IDs
    map = np.zeros(shape=data.shape,dtype=map_dtype)
    max_segs = np.iinfo(map.dtype).max
    mark = 1

    # loop and fill the map until all points are masked
    while mdata.mask.all() != True:
        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmax(),mdata.shape)
        mark_dmap(mdata,map,pt,mark,diag)
        mark = mark+1
        if mark > max_segs:
            raise OverflowError("Overflow in map, change map_dtype")
    return map
    
def mark_dmap(mdata,map,pt,mark,diag):
    """
    Mark downward-connected region on segment map starting at node pt

    Parameters:
    
    * mdata Masked data array
    * map   Array mapping out segments
    * pt    Index of starting node
    * mark  Integer to mark map with.

    Returns nothing but modified mdata mask and map

    """
    if mdata.mask[pt] == True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked


    Q = [pt]
    
    while Q:
        pt = Q.pop(0)
        v = mdata.data[pt]

        # Check all neightbors
        for new_pt in neighbors(pt,mdata.shape,diag):

            if mdata.mask[new_pt] == False and mdata[new_pt]<v:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return


# upward versions of above functions for finding negative peaks.

def find_upward(data,pt,thres,diag=False):
    """
    Find points upward-connected to pt in data.

    Parameters:

    * data  2D array of data
    * pt    Starting point, bottom of peak, tuple.
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Return: list of indicies of upward connected nodes.

    """
    if type(pt) == int:
        pt = (pt,)
    pt = tuple(pt)
    shape = data.shape

    if data[pt] < thres:    # check that the initial point is below threshold.
        return []

    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node

        for new_pt in neighbors(pt,shape,diag):  # check all neightbors

            if thres>data[new_pt]>v and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    
    return segment

def find_all_upward(data,thres,diag=False):  
    """
    Find all upward-connected segments in data

    Parameters:

    * data  Array of data
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Returns: (centers,segments)

    * centers   List of indicies of local minima in each segment
    * segments  List of all points in a given segment

    """

    mdata = np.ma.masked_greater(data,thres)

    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all()!=True:

        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmin(),data.shape)

        # find all nodes downward connected to pt and mark them
        segment = find_upward(data,pt,thres,diag)

        for i in segment:   # probably a better way of doing this
            mdata[i] = ma.masked
        
        centers.append(pt)
        segments.append(segment)

    return centers,segments


def map_upward(data,thres,diag=False,map_dtype='u2'):
    """
    Create a map of upward-connected segments of data below threshold

    Parameters:

    * data  Array of data
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.
    * map_dtype Dtype of map.

    Returns: map

    map is an array of integers indicating segments. 0 represents noise,1 the
    first segment, 2 the second, etc.  The index value of these regions
    can be found using np.where

    """
    # create the masked data array
    mdata = np.ma.masked_greater(data,thres)

    # build the map to store segment IDs
    map = np.zeros(shape=data.shape,dtype=map_dtype)
    max_segs = np.iinfo(map.dtype).max
    mark = 1

    # loop and fill the map until all points are masked
    while mdata.mask.all() != True:
        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmin(),mdata.shape)
        mark_umap(mdata,map,pt,mark,diag)
        mark = mark+1
        if mark > max_segs:
            raise OverflowError("Overflow in map, change map_dtype")
    return map
    
def mark_umap(mdata,map,pt,mark,diag):
    """
    Mark upward-connected region on segment map starting at node pt

    Parameters:
    
    * mdata Masked data array
    * map   Array mapping out segments
    * pt    Index of starting node
    * mark  Integer to mark map with.

    Returns nothing but modified mdata mask and map

    """
    if mdata.mask[pt] == True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked

    Q = [pt]
    
    while Q:
        pt = Q.pop(0)
        v = mdata.data[pt]

        # Check all neightbors
        for new_pt in neighbors(pt,mdata.shape,diag):

            if mdata.mask[new_pt] == False and mdata[new_pt]>v:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return


# Connected segmenting method:
# The connected segmenting method uses the flood fill algorithm to find
# all points connected to an initial node which are above a given threshold.

# The nconnected method does the same except nodes must be below the
# threshold (used for finding negative peaks)

def find_connected(data,pt,thres,diag=False):
    """
    Find points connected to pt in data.

    Parameters:

    * data  2D array of data
    * pt    Starting point, (top of peak).
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Return: list of indicies of downward connected nodes.

    """
    if type(pt) == int:
        pt = (pt,)
    pt = tuple(pt)
    shape = data.shape

    if data[pt] < thres:    # check that the initial point is above threshold.
        return []

    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node

        for new_pt in neighbors(pt,shape,diag):  # check all neightbors

            if data[new_pt]>thres and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    
    return segment

def find_all_connected(data,thres,diag=False):  
    """
    Find all connected segments in data

    Parameters:

    * data  Array of data
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Returns: (centers,segments)

    * centers   List of indicies of local maximum in each segment
    * segments  List of all points in a given segment

    """
    mdata = np.ma.masked_less(data,thres)

    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all()!=True:

        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmax(),mdata.shape)

        # find all nodes downward connected to pt and mark them
        segment = find_connected(data,pt,thres,diag)

        for i in segment:   # probably a better way of doing this
            mdata[i] = ma.masked
        centers.append(pt)
        segments.append(segment)

    return centers,segments


def map_connected(data,thres,diag=False,map_dtype='u2'):
    """
    Create a map of connected segments of data above threshold

    Parameters:

    * data  Array of data
    * thres Threshold, below this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.
    * map_dtype Dtype of map.

    Returns: map

    map is an array of integers indicating segments. 0 represents noise,1 the
    first segment, 2 the second, etc.  The index value of these regions
    can be found using np.where

    """
    # create the masked data array
    mdata = np.ma.masked_less(data,thres)

    # build the map to store segment IDs
    map = np.zeros(shape=data.shape,dtype=map_dtype)
    max_segs = np.iinfo(map.dtype).max
    mark = 1

    # loop and fill the map until all points are masked
    while mdata.mask.all() != True:
        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmax(),mdata.shape)
        mark_cmap(mdata,map,pt,mark,diag)
        mark = mark+1
        if mark > max_segs:
            raise OverflowError("Overflow in map, change map_dtype")
    return map
    
def mark_cmap(mdata,map,pt,mark,diag):
    """
    Mark connected region on segment map starting at node pt

    Parameters:
    
    * mdata Masked data array
    * map   Array mapping out segments
    * pt    Index of starting node
    * mark  Integer to mark map with.

    Returns nothing but modified mdata mask and map

    """
    if mdata.mask[pt] == True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked
    Q = [pt]
    
    while Q:
        pt = Q.pop(0)

        if mdata.mask[pt]==False:
            map[pt] = mark
            mdata[pt] = ma.masked

        # Check all neightbors
        for new_pt in neighbors(pt,mdata.shape,diag):

            if mdata.mask[new_pt] == False:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return


def find_nconnected(data,pt,thres,diag=False):
    """
    Find points connected to pt in data below threshold.

    Parameters:

    * data  2D array of data
    * pt    Starting point, (top of peak).
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Return: list of indicies of connected nodes.

    """
    if type(pt) == int:
        pt = (pt,)
    pt = tuple(pt)
    shape = data.shape

    if data[pt] > thres:    # check that the initial point is above threshold.
        return []

    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node

        for new_pt in neighbors(pt,shape,diag):  # check all neightbors

            if data[new_pt]<thres and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    
    return segment

def find_all_nconnected(data,thres,diag=False):  
    """
    Find all connected segments in data below threshold.

    Parameters:

    * data  Array of data
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.

    Returns: (centers,segments)

    * centers   List of indicies of local maximum in each segment
    * segments  List of all points in a given segment

    """
    mdata = np.ma.masked_greater(data,thres)

    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all()!=True:

        # find the minimum (center of segment)
        pt = np.unravel_index(mdata.argmmin(),mdata.shape)

        # find all nodes downward connected to pt and mark them
        segment = find_nconnected(data,pt,thres,diag)

        for i in segment:   # probably a better way of doing this
            mdata[i] = ma.masked
        centers.append(pt)
        segments.append(segment)

    return centers,segments


def map_nconnected(data,thres,diag=False,map_dtype='u2'):
    """
    Create a map of connected segments of data below threshold

    Parameters:

    * data  Array of data
    * thres Threshold, above this nodes are considered noise.
    * diag  True or False to include diagonal neighbors in connection.
    * map_dtype Dtype of map.

    Returns: map

    map is an array of integers indicating segments. 0 represents noise,1 the
    first segment, 2 the second, etc.  The index value of these regions
    can be found using np.where

    """
    # create the masked data array
    mdata = np.ma.masked_greater(data,thres)

    # build the map to store segment IDs
    map = np.zeros(shape=data.shape,dtype=map_dtype)
    max_segs = np.iinfo(map.dtype).max
    mark = 1

    # loop and fill the map until all points are masked
    while mdata.mask.all() != True:
        # find the maximum (center of segment)
        pt = np.unravel_index(mdata.argmin(),mdata.shape)
        mark_nmap(mdata,map,pt,mark,diag)
        mark = mark+1
        if mark > max_segs:
            raise OverflowError("Overflow in map, change map_dtype")
    return map
    
def mark_nmap(mdata,map,pt,mark,diag):
    """
    Mark negatively connected region on segment map starting at node pt

    Parameters:
    
    * mdata Masked data array
    * map   Array mapping out segments
    * pt    Index of starting node
    * mark  Integer to mark map with.

    Returns nothing but modified mdata mask and map

    """
    if mdata.mask[pt] == True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked
    Q = [pt]
    
    while Q:
        pt = Q.pop(0)

        if mdata.mask[pt]==False:
            map[pt] = mark
            mdata[pt] = ma.masked

        # Check all neightbors
        for new_pt in neighbors(pt,mdata.shape,diag):

            if mdata.mask[new_pt] == False:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return
