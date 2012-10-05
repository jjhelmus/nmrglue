"""
Peak picking routines, lineshape parameter guessing, and related functions.
"""

import numpy as np
import scipy.ndimage as ndimage

from .analysisbase import ndwindow_index, valid_pt
from .lineshapes1d import gauss, ls_str2class
from .segmentation import find_all_downward, find_all_upward
from .segmentation import find_all_connected, find_all_nconnected
from ..fileio import table


def pick(data, pthres, nthres=None, msep=None, algorithm='connected',
    est_params=True, lineshapes=None, edge=None, diag=False, c_struc=None,
    c_ndil=0, cluster=True, table=True, axis_names=['A', 'Z', 'Y', 'X']):
    """
    Pick (find) peaks in a region of a NMR spectrum.

    Parameters
    ----------
    data : ndarray
        Region of NMR spectrum to pick peaks from.
    pthres : float
        Minimum peak height for positive peaks. None to not detect positive
        peaks.
    nthres : float
        Minimum peak height for negative peaks (typically a negative value).
        None to not detect negative peaks.
    msep : tuple of ints, optional
        N-tuple of minimum peak seperations along each axis. Must be provided
        if algorithm is 'thresh' or 'thresh-fast'.
    algorithm : {'thres', thresh-fast', 'downward', 'connected'}, optional
        Peak picking algorithm to use.  Default is 'connected'.
    est_params : bool, optional
        True to perform an estimate of linewidths and amplitude for all peaks
        picked.  False, the default, will return only the peak locations.
    lineshapes : list, optional
        A list of lineshape classes or string shortcuts for each dimension.
        If not specified Gaussian type lineshapes with a FWHM  linewidth
        parameter is assumed in each dimension. This parameter if only used
        if est_params is True.
    edge : tuple of ints, optional
        Tuple to add to peak locations representing the edge of the region.
        None, the default, skips this addition.
    diag : bool, optional
        True to consider diagonal points to be  touching in peak finding
        algorithm and clustering.
    c_struc : ndarray, optional
        Structure element to use when applying dilation on segments before
        applying clustering algorithm. None will apply a default square
        structure with connectivity one will be applied.
    c_ndil : int, optional
        Number of dilations to perform on segments before applying clustering
        algorithm.
    cluster : bool, optional
        True to cluster touching peaks. False does not apply clustering.
    table : bool, optional
        True to return a table. False to return lists.
    axis_names : list. optional
        List of axis names, the last n will be used for column name prefixes
        in table where n is the number of dimensions.

    Returns
    -------
    locations : list, returned when table is False
        Peak locations.
    cluster_ids : list, returned when table is False and cluster is True
        Cluster numbers for peaks.
    scales : list, returned when table is False and est_params is True
        Estimated peak scales (linewidths).
    amps : list, returned when table is False and est_params is True
        Estimated peak amplitudes.
    table : recarray, returned when table is True
        Table of request peak parameters.

    """
    ####################
    # Check parameters #
    ####################
    ndim = len(data.shape)

    # check msep
    if type(msep) == int:
        msep = (msep, )
    if algorithm in ['thres', 'thres-fast'] and len(msep) != ndim:
        raise ValueError("msep has incorrect length")

    # check algorithm
    if algorithm not in ['thres', 'thres-fast', 'downward', 'connected']:
        raise ValueError('Invalid algorithm %s' % (algorithm))

    # check  lineshapes
    if est_params:
        # expand None
        if lineshapes == None:
            lineshapes = [gauss() for i in range(ndim)]
        ls_classes = []

        # replace strings
        for l in lineshapes:
            if type(l) is str:
                ls_classes.append(ls_str2class(l))
            else:
                ls_classes.append(l)
        # check that all classes have 2 parameters
        for i, ls in enumerate(ls_classes):
            if ls.nparam(10) != 2:
                s = "Lineshape class %i does not have two parameters"
                raise ValueError(s % (i))

        if len(ls_classes) != ndim:
            raise ValueError("Incorrect number of lineshapes")

    if edge != None and len(edge) != ndim:
        raise ValueError("edge has incorrect length")

    #######################
    # find positive peaks #
    #######################
    if pthres == None:    # no locations
        ploc = []
        pseg = []

    elif est_params == True:  # find locations and segments
        if algorithm == 'thres':
            ploc, pseg = find_all_thres_fast(data, pthres, msep, True)
        elif algorithm == 'thres-fast':
            ploc, pseg = find_all_thres_fast(data, pthres, msep, True)
        elif algorithm == 'downward':
            ploc, pseg = find_all_downward(data, pthres, True, diag)
        elif algorithm == 'connected':
            ploc, pseg = find_all_connected(data, pthres, True, diag)
        else:
            raise ValueError('Invalid algorithm %s' % (algorithm))

    else:   # find only locations
        if algorithm == 'thres':
            ploc = find_all_thres_fast(data, pthres, msep, False)
        elif algorithm == 'thres-fast':
            ploc = find_all_thres_fast(data, pthres, msep, False)
        elif algorithm == 'downward':
            ploc = find_all_downward(data, pthres, False, diag)
        elif algorithm == 'connected':
            ploc = find_all_connected(data, pthres, False, diag)
        else:
            raise ValueError('Invalid algorithm %s' % (algorithm))

    #######################
    # find negative peaks #
    #######################
    if nthres == None:    # no locations
        nloc = []
        nseg = []

    elif est_params == True:  # find locations and segments
        if algorithm == 'thres':
            nloc, nseg = find_all_nthres(data, nthres, msep, True)
        elif algorithm == 'thres-fast':
            nloc, nseg = find_all_nthres_fast(data, nthres, msep, True)
        elif algorithm == 'downward':
            nloc, nseg = find_all_upward(data, nthres, True, diag)
        elif algorithm == 'connected':
            nloc, nseg = find_all_nconnected(data, nthres, True, diag)
        else:
            raise ValueError('Invalid algorithm %s' % (algorithm))

    else:   # find only locations
        if algorithm == 'thres':
            nloc = find_all_nthres(data, nthres, msep, False)
        elif algorithm == 'thres-fast':
            nloc = find_all_nthres_fast(data, nthres, msep, False)
        elif algorithm == 'downward':
            nloc = find_all_upward(data, nthres, False, diag)
        elif algorithm == 'connected':
            nloc = find_all_nconnected(data, nthres, False, diag)
        else:
            raise ValueError('Invalid algorithm %s' % (algorithm))

    # combine the positive and negative peaks
    locations = ploc + nloc

    #########################################################
    # return locations if no parameter estimation requested #
    #########################################################
    if est_params == False:
        if cluster:     # find clusters
            cluster_ids = clusters(data, locations, pthres, nthres, c_struc,
                                    None, c_ndil)
            locations = add_edge(locations, edge)
            if table:
                return pack_table(locations, cluster_ids,
                                    axis_names=axis_names)
            else:
                return locations, cluster_ids
        else:   # Do not determine clusters
            locations = add_edge(locations, edge)
            if table:
                return pack_table(locations, axis_names=axis_names)
            else:
                return locations

    ##################################
    # estimate scales and amplitudes #
    ##################################
    seg_slices = pseg + nseg
    scales = [[]] * len(locations)
    amps = [[]] * len(locations)
    #scales = np.zeros(np.array(locations).shape,dtype=float)
    #amps = np.zeros(len(locations),dtype=float)

    for i, (l, seg_slice) in enumerate(zip(locations, seg_slices)):
        null, scales[i], amps[i] = guess_params_slice(data, l, seg_slice,
                                                      ls_classes)

    ########################################################
    # return locations, scales and amplitudes as requested #
    ########################################################
    if cluster:
        cluster_ids = clusters(data, locations, pthres, nthres, c_struc, None,
                               c_ndil)
        locations = add_edge(locations, edge)
        if table:
            return pack_table(locations, cluster_ids, scales, amps, axis_names)
        else:
            return locations, cluster_ids, scales, amps
    else:
        locations = add_edge(locations, edge)
        if table:
            return pack_table(locations, scales=scales, amps=amps,
                              axis_names=axis_names)
        else:
            return locations, scales, amps


def add_edge(locations, edge):
    """
    Add edge to list of locations, returning a list of edge-added locations
    """
    if edge is None:
        return locations
    return [tuple([i + j for i, j in zip(edge, l)]) for l in locations]


def clusters(data, locations, pthres, nthres, d_struc=None, l_struc=None,
        ndil=0):
    """
    Perform cluster analysis of peak locations.

    Parameters
    ----------
    data : ndarray
        Array of data which has been peak picked.
    locations : list
        List of peak locations.
    pthres : float
        Postive peak threshold. None for no postive peaks.
    nthres : float
        Negative peak threshold. None for no negative peaks.
    d_struc : ndarray, optional
        Structure of binary dilation to apply on segments before clustering.
        None uses a square structure with connectivity of one.
    l_struc : ndarray, optional
        Structure to use for determining segment connectivity in clustering.
        None uses square structure with connectivity of one.
    dnil : int, optional
        Number of dilation to apply on segments before determining clusters.

    Returns
    -------
    cluster_ids : list
        List of cluster number corresponding to peak locations.

    """
    # make a binary array of regions above/below the noise thresholds
    if pthres == None:  # negative peaks only
        input = data < nthres
    elif nthres == None:  # postive peaks only
        input = data > pthres
    else:               # both positive and negative
        input = np.bitwise_or(data < nthres, data > pthres)

    # apply dialations to these segments
    if ndil != 0:
        input = ndimage.binary_dilation(input, d_struc, iterations=ndil)

    # label this array, these are the clusters.
    labeled_array, num_features = ndimage.label(input, l_struc)

    return [labeled_array[i] for i in locations]


def pack_table(locations, cluster_ids=None, scales=None, amps=None,
        axis_names=["A", "Z", "Y", "X"]):
    """
    Create a table from peak information.

    Parameters
    ----------
    locations : list
        List of peak locations.
    cluster_ids : list, optional
        List of cluster numbers. None will not include cluster number in the
        table.
    scales : list, optional
        List of peak scales (linewidths). None will not include peak scales in
        the table.
    amps : list, optional
        List of peak amplitudes. None will not include peak amplitudes in the
        table.
    axis_names : list, optional
        List of axis names, the last n will be used for column name prefixes
        where n is the number of dimensions.

    Returns
    -------
    table : recarray
        nmrglue table with column representing peak parameters. Peak locations
        are given column names like 'X_AXIS', 'Y_AXIS', etc. Cluster_ids are
        given a column name of 'cID'. Peak scales (linewidths) are given
        column names like 'X_LW','Y_LW'.  Peak amplitudes are given a column
        name of 'VOL'.

    """
    ndim = len(locations[0])
    anames = axis_names[-ndim:]

    dt = [(a + "_AXIS", np.float) for a in anames]
    rec = np.rec.array(locations, dtype=dt)

    if cluster_ids != None:
        rec = table.append_column(rec, cluster_ids, 'cID', 'int')
    if scales != None:
        names = [a + "_LW" for a in anames]
        for n, c in zip(names, np.array(scales).T):
            rec = table.append_column(rec, c, n, 'float')
    if amps != None:
        rec = table.append_column(rec, amps, 'VOL', 'float')

    return rec


def guess_params_slice(data, location, seg_slice, ls_classes):
    """
    Guess the parameter of a peak in a segment.

    Parameters
    ----------
    data : ndarray
        NMR data.
    location : tuple
        Peak locations.
    seg_slice : list of slices
        List slices which slice data to give the desired segment.
    lineshapes : list
        List of lineshape classes.

    Returns
    -------
    location : list
        Peak locations.
    scale : list
        Peak scales (linewidths).
    amp : list
        Peak amplitudes.

    """
    # find the rectangular region around the segment
    region = data[seg_slice]
    edge = [s.start for s in seg_slice]
    rlocation = [l - s.start for l, s in zip(location, seg_slice)]

    # amptide is estimated by the sum of all points in region
    amp = np.sum(region)

    scale = []    # list of linewidths
    nlocation = []    # list of peak centers

    # loop over the axes
    for axis, ls in enumerate(ls_classes):
        # create the 1D lineshape
        r = extract_1d(region, rlocation, axis)
        # estimate the linewidth
        loc, sc = ls.guessp(r)
        scale.append(float(sc))
        nlocation.append(float(loc))

    return tuple([l + e for l, e in zip(nlocation, edge)]), tuple(scale), amp


def extract_1d(data, location, axis):
    """
    Extract a 1D slice from data along axis at location
    """
    s = [slice(v, v + 1) for v in location]
    s[axis] = slice(None, None)
    return np.atleast_1d(np.squeeze(data[s]))


# algorithm specific peak picking routines
def find_all_thres(data, thres, msep, find_segs=False):
    """
    Peak pick a spectrum using a threshhold-minimum distance algorithm.

    Find peaks (local maxima) in a arbitrary dimensional NMR spectra above a
    set threshold with a minimal distance between peaks.  When the spectrum is
    small and multiple copies can fit into RAM use the _fast version of this
    function. Segments are found by finding the first point in each direction
    along each dimension which is below the threshold.

    Parameters
    ----------
    data : ndarray
        NMR data.
    thres : float
        Threshold value for minimum peak height
    msep : tuple
        Tuple of minimum peak seperations along each axis.
    find_segs : bool, optional
        True  to find segments and return a list of slices which select that
        segment.  False performs no segmentation discovery.

    Returns
    -------
    locations : list
        List of peak locations
    seg_slices : list, optional
        List of slices which extract a region around each peak. Only returned
        when find_segs is True.

    """
    locations = []  # create an empty list of peak locations
    wsize = tuple([2 * i + 1 for i in msep])  # window size is 2*seperation+1

    # loop over the windows
    for idx, s in ndwindow_index(data.shape, wsize):
        max = data[s].max()
        if max == data[idx] and max > thres:
            locations.append(idx)
    if find_segs:
        seg_slices = find_pseg_slice(data, locations, thres)
        return locations, seg_slices
    else:
        return locations


def find_all_nthres(data, thres, msep, find_segs=False):
    """
    Peak pick a spectrum using a threshhold-minimum distance algorithm.

    Identical to find_all_thres except local minima are found below the
    given threshold.  See :py:func:`find_all_thres` for a description of the
    algorithm and documentation.

    """
    locations = []  # create an empty list of peak locations
    wsize = tuple([2 * i + 1 for i in msep])  # window size is 2*seperation+1

    # loop over the windows
    for idx, s in ndwindow_index(data.shape, wsize):
        min = data[s].min()
        if min == data[idx] and min < thres:
            locations.append(idx)
    if find_segs:
        seg_slices = find_pseg_slice(data, locations, thres)
        return locations, seg_slices
    else:
        return locations


def find_all_thres_fast(data, thres, msep, find_segs=False):
    """
    Fast version of find_all_thres. See :py:func:`find_all_thres`.
    """
    wsize = tuple([2 * i + 1 for i in msep])  # window size is 2*seperation+1

    # find local maxima mask
    mx = ndimage.maximum_filter(data, size=wsize, mode='constant') == data

    # find positive threshold mask
    pthres = np.ma.greater(data, thres)

    # peaks are bitwise and of maximum mask and threshold mask
    locations = np.transpose(np.nonzero(np.bitwise_and(pthres, mx)))
    locations = [tuple(i) for i in locations]

    if find_segs:
        seg_slices = [find_pseg_slice(data, l, thres) for l in locations]
        return locations, seg_slices
    else:
        return locations


def find_all_nthres_fast(data, thres, msep, find_segs=False):
    """
    Fast version of find_all_nthres_fast. See :py:func:`find_all_thres`.
    """
    wsize = tuple([2 * i + 1 for i in msep])  # window size is 2*seperation+1

    # find local maxima mask
    mn = ndimage.minimum_filter(data, size=wsize, mode='constant') == data

    # find positive threshold mask
    nthres = np.ma.less(data, thres)

    # peaks are bitwise and of maximum mask and threshold mask
    locations = np.transpose(np.nonzero(np.bitwise_and(nthres, mn)))
    locations = [tuple(i) for i in locations]

    if find_segs:
        seg_slices = [find_pseg_slice(data, l, thres) for l in locations]
        return locations, seg_slices
    else:
        return locations


def find_pseg_slice(data, location, thres):
    """
    Find slices which define a segment in data above thres.
    """
    shape = data.shape
    seg_slice = []
    for dim, v in enumerate(location):
        # find start value
        al = list(location)
        start = v
        while(valid_pt(al, shape) and data[tuple(al)] > thres):
            start = start - 1
            al[dim] = start
        # find stop value
        al = list(location)
        stop = v
        while(valid_pt(al, shape) and data[tuple(al)] > thres):
            stop = stop + 1
            al[dim] = stop
        seg_slice.append(slice(start + 1, stop))
    return seg_slice


def find_nseg_slice(data, location, thres):
    """
    Find slices which define a segment in data below thres.
    """
    shape = data.shape
    seg_slice = []
    for dim, v in enumerate(location):
        # find start value
        al = list(location)
        start = v
        while(valid_pt(al, shape) and data[tuple(al)] < thres):
            start = start - 1
            al[dim] = start
        # find stop value
        al = list(location)
        stop = v
        while(valid_pt(al, shape) and data[tuple(al)] < thres):
            stop = stop + 1
            al[dim] = stop
        seg_slice.append(slice(start + 1, stop))
    return seg_slice
