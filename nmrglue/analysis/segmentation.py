"""
Functions to perform segmentation of NMR spectrum.
"""

import numpy as np
import numpy.ma as ma
import scipy.ndimage as ndimage

from .analysisbase import neighbors

# Connected segmenting method:
# The connected segmentation method finds all nodes which are above a given
# threshold and connected to the initial point.  For finding all segments
# the scipy.ndimage.label function is used for speed.


def label_connected(data, thres, structure):
    """
    Label connected features in data.  Returns labeled_array, num_features
    """
    return ndimage.label(data > thres, structure)


def find_all_connected(data, thres, find_segs=False, diag=False):
    """
    Find all connected segments.

    Parameters
    ----------
    data : ndarray
        Data to perform segmentation on.
    thres : float
        Threshold, below this nodes are considered noise.
    find_segs : bool, optional
        True to return a list of slices for the segments.
    diag : bool
        True to include diagonal neighbors in connection.

    Returns
    -------
    locations : list
        List of indices of local maximum in each segment.
    seg_slices : list, optional
        List of slices which extract a given segment from the data. Only
        returned when fig_segs is True.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    # determine labeled array of segments
    labels, num_features = label_connected(data, thres, structure)

    # determine locations of segment maxima
    locations = ndimage.maximum_position(data, labels, range(1, num_features +
                                                             1))
    # find segment slices if requested and return
    if find_segs is True:
        seg_slices = ndimage.find_objects(labels)
        return locations, seg_slices
    else:
        return locations

# nconnected method:
# The nconnected method is identical to the connected method except nodes must
# be below the threshold and local minimum are reported.  This is useful for
# finding negative peaks by setting thres to the negative of the noise level.


def label_nconnected(data, thres, structure):
    """
    Label nconnected features in data.  Returns labeled_array, num_features
    """
    return ndimage.label(data < thres, structure)


def find_all_nconnected(data, thres, find_segs=False, diag=False):
    """
    Find all negatively connected segments in data.

    Parameters
    ----------
    data : ndarray
        Data to perform segmentation on.
    thres : float
        Threshold, below this nodes are considered noise.
    find_segs : bool, optional
        True to return a list of slices for the segments.
    diag : bool
        True to include diagonal neighbors in connection.

    Returns
    -------
    locations : list
        List of indices of local maximum in each segment.
    seg_slices : list, optional
        List of slices which extract a given segment from the data. Only
        returned when fig_segs is True.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    # determine labeled array of segments
    labels, num_features = label_nconnected(data, thres, structure)

    # determine locations of segment maxima
    locations = ndimage.minimum_position(data, labels, range(1,
                                         num_features + 1))
    # find segment slices if requested and return
    if find_segs is True:
        seg_slices = ndimage.find_objects(labels)
        return locations, seg_slices
    else:
        return locations

# downward segmentation method:
# The downward segmenting method uses the flood fill algorithm to find
# all points connected to an initial node which are above a given threshold
# and to which a path exists in which each step of the path moves lower in
# intensity.  This can be though of as all points accessible by a water drop
# following downward slopes from the initial node.

# Upward segmentation uses the same principle except nodes must be below
# the threshold an upward path must exist.


def mark_dseg(mdata, map, pt, mark, structure):
    """
    Mark downward-connected region on segment map starting at node pt.

    Modifies mdata mask and map.

    Parameters
    ----------
    mdata : masked ndarray
        Masked data array.
    map :
        Array mapping out segments.
    pt : tuple of ints
        Index of starting node
    mark : int
        Integer to mark map with.

    """
    if mdata.mask[pt] is True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked

    Q = [pt]
    while Q:
        pt = Q.pop(0)
        v = mdata.data[pt]
        # Check all neighbors
        for new_pt in neighbors(pt, mdata.shape, structure):
            if mdata.mask[new_pt] == False and mdata[new_pt] < v:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return


def label_downward_seg(data, labels, seg_slice, seg_index, max_index,
                       structure):
    """ Label a segment which is downward connected """
    slabels = labels[seg_slice]
    msdata = np.ma.masked_array(data[seg_slice], mask=(slabels != seg_index))

    # mark the downward connected segment with the highest peak in the
    # selected segment with the segment index.
    argmax = np.unravel_index(msdata.argmax(), msdata.shape)
    mark_dseg(msdata, slabels, argmax, seg_index, structure)

    # mark any
    while msdata.mask.all() == False:
        argmax = np.unravel_index(msdata.argmax(), msdata.shape)
        mark_dseg(msdata, slabels, argmax, max_index, structure)
        max_index = max_index + 1
    return max_index


def label_downward(data, thres, structure):
    """
    Label connected features in data. Returns labeled_array, num_features
    """
    # find connected segments
    labels, num_features = ndimage.label(data > thres, structure)
    seg_slices = ndimage.find_objects(labels)
    max_index = int(num_features + 1)

    # loop over the segments and perform a downward segment on each
    for i, s in enumerate(seg_slices):
        max_index = label_downward_seg(data, labels, s, i + 1, max_index,
                                       structure)
    return labels, max_index - 1


def find_all_downward(data, thres, find_segs=False, diag=False):
    """
    Find all downward connected segments in data

    Parameters
    ----------
    data : ndarray
        Data to perform segmentation on.
    thres : float
        Threshold, below this nodes are considered noise.
    find_segs : bool, optional
        True to return a list of slices for the segments.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    locations : list
        List of indices of local maximum in each segment.
    seg_slices : list, optional
        List of slices which extract a given segment from the data. Only
        returned when fig_segs is True.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    # determine labeled array of segments
    labels, num_features = label_downward(data, thres, structure)

    # determine locations of segment maxima
    locations = ndimage.maximum_position(data, labels, range(1,
                                         num_features + 1))

    # find segment slices if requested and return
    if find_segs is True:
        seg_slices = ndimage.find_objects(labels)
        return locations, seg_slices
    else:
        return locations


def mark_useg(mdata, map, pt, mark, structure):
    """
    Mark upward-connected region on segment map starting at node pt

    Modifies mdata mask and map.

    Parameters
    ----------
    mdata : masked ndarray
        Masked data array.
    map :
        Array mapping out segments.
    pt : tuple of ints
        Index of starting node
    mark : int
        Integer to mark map with.

    """
    if mdata.mask[pt] is True:
        return
    else:
        map[pt] = mark
        mdata[pt] = ma.masked

    Q = [pt]
    while Q:
        pt = Q.pop(0)
        v = mdata.data[pt]
        # Check all neighbors
        for new_pt in neighbors(pt, mdata.shape, structure):
            if mdata.mask[new_pt] == False and mdata[new_pt] > v:
                Q.append(new_pt)
                map[new_pt] = mark
                mdata[new_pt] = ma.masked
    return


def label_upward_seg(data, labels, seg_slice, seg_index, max_index,
                     structure):
    """ Label a segment which is upward connected """
    slabels = labels[seg_slice]
    msdata = np.ma.masked_array(data[seg_slice],
                                mask=(slabels != seg_index))
    # mark the upward connected segment with the highest peak in the
    # selected segment with the segment index.
    argmin = np.unravel_index(msdata.argmin(), msdata.shape)
    mark_useg(msdata, slabels, argmin, seg_index, structure)

    # mark any
    while msdata.mask.all() == False:
        argmin = np.unravel_index(msdata.argmin(), msdata.shape)
        mark_useg(msdata, slabels, argmin, max_index, structure)
        max_index = max_index + 1

    return max_index


def label_upward(data, thres, structure):
    """
    Label upward connected features in data. Returns labeled_array,
    num_features
    """
    # find connected segments
    labels, num_features = ndimage.label(data < thres, structure)
    seg_slices = ndimage.find_objects(labels)
    max_index = int(num_features + 1)

    # loop over the segments and perform a downward segment on each
    for i, s in enumerate(seg_slices):
        max_index = label_upward_seg(data, labels, s, i + 1, max_index,
                                     structure)
    return labels, max_index - 1


def find_all_upward(data, thres, find_segs=False, diag=False):
    """
    Find all upward connected segments in data

    Parameters
    ----------
    data : ndarray
        Data to perform segmentation on.
    thres : float
        Threshold, below this nodes are considered noise.
    find_segs : bool, optional
        True to return a list of slices for the segments.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    locations : list
        List of indices of local maximum in each segment.
    seg_slices : list, optional
        List of slices which extract a given segment from the data. Only
        returned when fig_segs is True.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    # determine labeled array of segments
    labels, num_features = label_upward(data, thres, structure)

    # determine locations of segment maxima
    locations = ndimage.minimum_position(data, labels,
                                         range(1, num_features + 1))
    # find segment slices if requested and return
    if find_segs is True:
        seg_slices = ndimage.find_objects(labels)
        return locations, seg_slices
    else:
        return locations


##########################
# Single point functions #
##########################


def find_downward(data, pt, thres, diag=False):
    """
    Find points downward-connected to a point in data.

    Parameters
    ----------
    data : ndarray
        Array of data.
    pt : tuple of ints
        Starting point of peak.
    thres : float
        Threshold, below this nodes are considered noise.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    nodes : list
        Indices of downward-connected nodes.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    if isinstance(pt, int):
        pt = (pt, )
    pt = tuple(pt)
    shape = data.shape

    if data[pt] < thres:    # check that the initial point is above threshold.
        return []
    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node
        for new_pt in neighbors(pt, shape, structure):  # check all neighbors
            if thres < data[new_pt] < v and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    return segment


def find_connected(data, pt, thres, diag=False):
    """
    Find points connected to a point in data.

    Parameters
    ----------
    data : ndarray
        Array of data.
    pt : tuple of ints
        Starting point of peak.
    thres : float
        Threshold, below this nodes are considered noise.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    nodes : list
        Indices of connected nodes.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    if isinstance(pt, int):
        pt = (pt, )
    pt = tuple(pt)
    shape = data.shape

    if data[pt] < thres:    # check that the initial point is above threshold.
        return []
    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        for new_pt in neighbors(pt, shape, structure):  # check all neighbors
            if data[new_pt] > thres and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    return segment


def find_nconnected(data, pt, thres, diag=False):
    """
    Find points connected to pt in data below threshold.

    Parameters
    ----------
    data : ndarray
        Array of data.
    pt : tuple of ints
        Starting point of peak.
    thres : float
        Threshold, above this nodes are considered noise.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    nodes : list
        Indices of connected nodes.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    if isinstance(pt, int):
        pt = (pt, )
    pt = tuple(pt)
    shape = data.shape

    if data[pt] > thres:    # check that the initial point is above threshold.
        return []
    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        for new_pt in neighbors(pt, shape, structure):  # check all neighbors
            if data[new_pt] < thres and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    return segment


def find_upward(data, pt, thres, diag=False):
    """
    Find points upward-connected to pt in data.

    Parameters
    ----------
    data : ndarray
        Array of data.
    pt : tuple of ints
        Starting point of peak.
    thres : float
        Threshold, below this nodes are considered noise.
    diag : bool, optional
        True to include diagonal neighbors in connection.

    Returns
    -------
    nodes : list
        Indices of upward-connected nodes.

    """
    # build structure array for defining feature connections
    ndim = data.ndim
    if diag:
        structure = ndimage.generate_binary_structure(ndim, ndim)
    else:
        structure = ndimage.generate_binary_structure(ndim, 1)

    if isinstance(pt, int):
        pt = (pt, )
    pt = tuple(pt)
    shape = data.shape

    if data[pt] > thres:    # check that the initial point is below threshold.
        return []
    Q = [pt]    # queue
    segment = [pt]

    while Q:    # loop until Q is empty
        pt = Q.pop(0)   # remove first element of queue
        v = data[pt]    # value at current node
        for new_pt in neighbors(pt, shape, structure):  # check all neighbors
            if thres > data[new_pt] > v and new_pt not in segment:
                Q.append(new_pt)
                segment.append(new_pt)
    return segment
