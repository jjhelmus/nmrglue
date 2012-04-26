"""
Peak picking routines, lineshape parameter guessing, and related functions.

"""

# external modules
import numpy as np
import scipy.ndimage as ndimage

# analysisbase fuctions
from analysisbase import ndwindow_index,valid_pt

# lineshape classes
from lineshapes1d import gauss, ls_str2class


# spectral segmentation functions
from segmentation import find_all_downward, find_all_upward
from segmentation import find_all_connected, find_all_nconnected 

from ..fileio import table

def pick(data,pthres,nthres=None,msep=None,algorithm='connected',
            est_params=True,lineshapes=None,edge=None,diag=False,c_struc=None,
            c_ndil=0,cluster=True,table=True,axis_names=['A','Z','Y','X']):
    """
    Pick (find) peaks in a spectral region. 

    Parameters:

    * data          N-dimensional array to pick peaks in.
    * pthres        Minimum peak height for positive peaks. Set to None to not
                    detect positive peaks.
    * nthres        Minimum peak height for negative peaks (typically a 
                    negative value).  Set to None to not detect negative peaks.
    * msep          N-tuple of minimum peak seperations along each axis.
                    Must be defined if algorithm is 'thresh' or 'thresh-fast'
    * algorithm     Peak picking algorithm to use.  Options are 'thres',
                    'thres-fast', 'downward', or 'connected'
    * est_params    Set to True to perform a rough estimate of linewidths and
                    amplitude for all peaks picked.  False returns only the
                    peak locations.
    * lineshapes    A list of lineshape classes or string shortcuts for each 
                    dimension.  If not specified Gaussian type lineshapes with 
                    a FWHM  linewidth parameter is assumed in each dimension.  
                    This parameter if only used if est_params is True.
    * edge          Tuple to add to peak locations representing the edge of a
                    slices region.  None skips this addition.
    * diag          Set True to consider diagonal points to be  touching in 
                    peak finding algorithm and clustering.
    * c_struc       Structure element to use when applying dilation on segments
                    before applying clustering algorithm. None will apply 
                    default square structure with connectivity one will be 
                    used.
    * c_ndil        Number of dilations to perform on segments before applying
                    clustering algorithm.
    * cluster       Set True to cluster touching peaks.
    * table         Set True to return turn a table.
    * axis_names    List of axis names, the last n will be used for column
                    name prefixes in table where n is the number of dimensions.

    Returns:    locations,[cluster_ids,[scales,amps]] or table

    * locations
    * cluster_ids
    * scales
    * amps

    * table

    """
    ####################
    # Check parameters #
    ####################
    ndim = len(data.shape)
    
    # check msep
    if type(msep) == int:
        msep = (msep,)
    if algorithm in ['thres','thres-fast'] and len(msep) != ndim:
        raise ValueError("msep has incorrect length")
 
    # check algorithm
    if algorithm not in ['thres','thres-fast','downward','connected']:
        raise ValueError('Invalid algorithm %s'%(algorithm))
   
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
        for i,ls in enumerate(ls_classes):
            if ls.nparam(10) != 2:
                s = "Lineshape class %i does not have two parameters"
                raise ValueError(s%(i))
        
        if len(ls_classes) != ndim:
            raise ValueError("Incorrect number of lineshapes")

    if edge!=None and len(edge)!=ndim:
        raise ValueError("edge has incorrect length")


    #######################
    # find positive peaks #
    #######################
    if pthres==None:    # no locations
        ploc = []
        pseq = []

    elif est_params==True:  # find locations and segments
        if algorithm == 'thres':
            ploc,pseg = find_all_thres_fast(data,pthres,msep,True)
        elif algorithm == 'thres-fast':
            ploc,pseg = find_all_thres_fast(data,pthres,msep,True)
        elif algorithm == 'downward':
            ploc,pseg = find_all_downward(data,pthres,True,diag)
        elif algorithm == 'connected':
            ploc,pseg = find_all_connected(data,pthres,True,diag)
        else:
            raise ValueError('Invalid algorithm %s'%(algorithm))

    else:   # find only locations 
        if algorithm == 'thres':
            ploc = find_all_thres_fast(data,pthres,msep,False)
        elif algorithm == 'thres-fast':
            ploc = find_all_thres_fast(data,pthres,msep,False)
        elif algorithm == 'downward':
            ploc = find_all_downward(data,pthres,False,diag)
        elif algorithm == 'connected':
            ploc = find_all_connected(data,pthres,False,diag)
        else:
            raise ValueError('Invalid algorithm %s'%(algorithm))
    

    #######################
    # find negative peaks #
    #######################
    if nthres==None:    # no locations
        nloc = []
        nseg = []
    
    elif est_params==True:  # find locations and segments
        if algorithm == 'thres':
            nloc,nseg = find_all_nthres(data,nthres,msep,True)
        elif algorithm == 'thres-fast':
            nloc,nseg = find_all_nthres_fast(data,nthres,msep,True)
        elif algorithm == 'downward':
            nloc,nseg = find_all_upward(data,nthres,True,diag)
        elif algorithm == 'connected':
            nloc,nseg = find_all_nconnected(data,nthres,True,diag)
        else:
            raise ValueError('Invalid algorithm %s'%(algorithm))
    
    else:   # find only locations
        if algorithm == 'thres':
            nloc = find_all_nthres(data,nthres,msep,False)
        elif algorithm == 'thres-fast':
            nloc = find_all_nthres_fast(data,nthres,msep,False)
        elif algorithm == 'downward':
            nloc = find_all_upward(data,nthres,False,diag)
        elif algorithm == 'connected':
            nloc = find_all_nconnected(data,nthres,False,diag)
        else:
            raise ValueError('Invalid algorithm %s'%(algorithm))
       
    # combine the positive and negative peaks
    locations = ploc+nloc

    #########################################################
    # return locations if no parameter estimation requested #
    #########################################################
    if est_params==False:
        if cluster:     # find clusters
            cluster_ids = clusters(data,locations,pthres,nthres,c_struc,None,                              c_ndil)
            locations = add_edge(locations,edge)
            if table:
                return pack_table(locations,cluster_ids,axis_names=axis_names)
            else:
                return locations,cluster_ids
        else:   # Do not determine clusters
            locations = add_edge(locations,edge)
            if table:
                return pack_table(locations,axis_names=axis_names)
            else:
                return locations
    
    ##################################
    # estimate scales and amplitudes #
    ##################################
    seg_slices = pseg+nseg
    scales = [[]]*len(locations)
    amps = [[]] * len(locations)
    #scales = np.zeros(np.array(locations).shape,dtype=float)
    #amps = np.zeros(len(locations),dtype=float)

    for i,(l,seg_slice) in enumerate(zip(locations,seg_slices)):
        null,scales[i],amps[i]=guess_params_slice(data,l,seg_slice,ls_classes)
    
    ########################################################
    # return locations, scales and amplitudes as requested #
    ########################################################
    if cluster:
        cluster_ids = clusters(data,locations,pthres,nthres,c_struc,None,c_ndil)
        locations = add_edge(locations,edge)
        if table:
            return pack_table(locations,cluster_ids,scales,amps,axis_names)
        else:
            return locations,cluster_ids,scales,amps
    else:
        locations = add_edge(locations,edge)
        if table:
            return pack_table(locations,scales=scales,amps=amps,
                              axis_names=axis_names)
        else:
            return locations,scales,amps

def add_edge(locations,edge):
    """
    Add edge to list of locations, returning a list of edge-added locations
    """
    if edge != None:
        return [tuple([i+j for i,j in zip(edge,l)]) for l in locations]
    return locations

def clusters(data,locations,pthres,nthres,d_struc=None,l_struc=None,ndil=0):
    """
    Perform cluster analysis of peak locations

    Parameters:

    * data          Array of data which has been peak picked
    * locations     List of peak locations
    * pthres        Postive peak threshold or None for no postive peaks
    * nthres        Negative peak threshold or None for no negative peaks
    * d_struc       Structure of binary dilation to apply on segments before
                    clustering.  None uses square connectivity 1 structure.
    * l_struc       Structure to use for determining  segment connectivity
                    in clustering.  None uses square connectivity 1 structure.
    * dnil          Number of dilation to apply on segments before determining
                    clusters.

    Returns:    cluster_ids

    * cluster_ids   List of cluster_ids for each location in locations list.

    """
    # make a binary array of regions above/below the noise thresholds
    if pthres==None:    # negative peaks only
        input = data < nthres
    elif nthres==None:  # postive peaks only
        input = data > pthres
    else:               # both positive and negative
        input = np.bitwise_or(data < nthres,data > pthres)
    
    # apply dialations to these segments
    if ndil!=0:
        input = ndimage.binary_dilation(input,d_struc,iterations=ndil)

    # label this array, these are the clusters.
    labeled_array,num_features = ndimage.label(input,l_struc)
    
    return [labeled_array[i] for i in locations]


def pack_table(locations,cluster_ids=None,scales=None,amps=None,
                axis_names=["A","Z","Y","X"]):
    """
    Create a table from peak information.

    Parameters:

    * locations     List of peak locations.
    * cluster_ids   List of cluster numbers.
    * scales        List of peak scales (linewidths).
    * amps          List of peak amplitudes.
    * axis_names    List of axis names, the last n will be used for column
                    name prefixes where n is the number of dimensions.

    If any of cluster_ids, scales, or amps in None the corresponding columns
    will not be present in the table.

    Returns: table

    * table nmrglue table with column representing peak parameters.
            locations are given column names like 'X_AXIS', 'Y_AXIS', etc
            cluster_ids are given a column name of 'cID'
            scales are given column names like 'X_LW','Y_LW'
            amps are given a column name of 'VOL'

    """
    ndim = len(locations[0])
    anames = axis_names[-ndim:]
    
    dt = [(a+"_AXIS",np.float) for a in anames]
    rec = np.rec.array(locations,dtype=dt)

    if cluster_ids != None:
        rec = table.append_column(rec,cluster_ids,'cID','int')
    
    if scales != None:
        names = [a+"_LW" for a in anames]
        for n,c in zip(names,np.array(scales).T):
            rec = table.append_column(rec,c,n,'float')
    
    if amps != None:
        rec = table.append_column(rec,amps,'VOL','float')

    return rec



def guess_params_slice(data,location,seg_slice,ls_classes):
    """
    Guess the parameter of a peak in a segment given slices.

    Parameters:

    * data          Spectral data.
    * seg_slice     List slices which slice data to given the desired segment.
    * lineshapes    List of lineshape classes.

    Return: location,scale,amp

    * location  Peak location.
    * scale     Peak scale.
    * amp       Peak amplitude.

    """

    # find the rectangular region around the segment
    region= data[seg_slice]
    edge = [s.start for s in seg_slice]
    rlocation = [l-s.start for l,s in zip(location,seg_slice)]

    # amptide is estimated by the sum of all points in region
    amp = np.sum(region)

    scale  = []    # list of linewidths
    nlocation = []    # list of peak centers

    # loop over the axes
    for axis,ls in enumerate(ls_classes):
        # create the 1D lineshape
        r = extract_1d(region,rlocation,axis)
        # estimate the linewidth
        loc,sc = ls.guessp(r)
        scale.append(float(sc))
        nlocation.append(float(loc))

    return tuple([l+e for l,e in zip(nlocation,edge)]),tuple(scale),amp

def extract_1d(data,location,axis):
    """
    Extract a 1D slice from data along axis at location
    """
    s = [slice(v,v+1) for v in location]
    s[axis] = slice(None,None)
    return np.atleast_1d(np.squeeze(data[s]))


# algorithm specific peak picking routines

def find_all_thres(data,thres,msep,find_segs=False):
    """
    Peak pick a spectrum using a threshhold-minimum distance algorithm.
    
    Find peaks (local maxima) in a arbitrary dimensional NMR spectra above a 
    set threshold with a minimal distance between peaks.  When the spectrum is 
    small and multiple copies can fit into RAM use the _fast version of this 
    function. Segments are found by finding the first point in each direction
    along each dimension which is below the threshold.

    Parameters:

    * data      N-dimensional array.
    * thres     Threshold value for minimum peak height
    * msep      N-tuple of minimum peak seperations along each axis
    * find_segs True or False to return a list of slices for the segments.


    Returns: locations,seg_slices
    
    * locations     List of indicies of peak locations
    * seg_slices    List of slices which extract a region around each peak. 

    """
    locations = []  # create an empty list of peak locations
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    # loop over the windows
    for idx,s in ndwindow_index(data.shape,wsize):
        max = data[s].max()
        if max == data[idx] and max > thres:
            locations.append(idx)

    if find_segs:
        seg_slices = find_pseg_slices(data,locations,thres)
        return locations,seg_slices
    else:
        return locations


def find_all_nthres(data,thres,msep,find_segs=False):
    """
    Peak pick a spectrum using a threshhold-minimum distance algorithm.
    
    Identical to find_all_thres except local minima are found below the
    given threshold.  See find_all_thres for a description of the algorithm and
    parameters.

    """
    locations = []  # create an empty list of peak locations
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    # loop over the windows
    for idx,s in ndwindow_index(data.shape,wsize):
        min = data[s].min()
        if min == data[idx] and min < thres:
            locations.append(idx)

    if find_segs:
        seg_slices = find_pseg_slices(data,locations,thres)
        return locations,seg_slices
    else:
        return locations


def find_all_thres_fast(data,thres,msep,find_segs=False):
    """
    Fast version of find_all_thres.
    """
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    # find local maxima mask
    mx=ndimage.maximum_filter(data,size=wsize,mode='constant')==data
   
    # find positive threshold mask
    pthres = np.ma.greater(data,thres)
    
    # peaks are bitwise and of maximum mask and threshold mask
    locations = np.transpose(np.nonzero(np.bitwise_and(pthres,mx)))
    locations = [tuple(i) for i in locations]

    if find_segs:
        seg_slices = [find_pseg_slice(data,l,thres) for l in locations]
        return locations,seg_slices
    else:
        return locations


def find_all_nthres_fast(data,thres,msep,find_segs=False):
    """
    Fast version of find_all_nthres_fast.
    """
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    # find local maxima mask
    mn=ndimage.minimum_filter(data,size=wsize,mode='constant')==data
   
    # find positive threshold mask
    nthres = np.ma.less(data,thres)
    
    # peaks are bitwise and of maximum mask and threshold mask
    locations = np.transpose(np.nonzero(np.bitwise_and(nthres,mn)))
    locations = [tuple(i) for i in locations]

    if find_segs:
        seg_slices = [find_pseg_slice(data,l,thres) for l in locations]
        return locations,seg_slices
    else:
        return locations


def find_pseg_slice(data,location,thres):
    """
    Find slices which define a segment in data above thres starting a location
    """
    shape = data.shape
    seg_slice = []
    for dim,v in enumerate(location):
        
        # find start value
        al = list(location)
        start = v
        while(valid_pt(al,shape) and data[tuple(al)]>thres):
            start = start-1
            al[dim] = start

        # find stop value
        al = list(location)
        stop = v
        while(valid_pt(al,shape) and data[tuple(al)]>thres):
            stop = stop+1
            al[dim] = stop
        
        seg_slice.append(slice(start+1,stop))

    return seg_slice


def find_nseg_slice(data,location,thres):
    """
    Find slices which define a segment in data below thres starting a location
    """
    shape = data.shape
    seg_slice = []
    for dim,v in enumerate(location):
        
        # find start value
        al = list(location)
        start = v
        while(valid_pt(al,shape) and data[tuple(al)]<thres):
            start = start-1
            al[dim] = start

        # find stop value
        al = list(location)
        stop = v
        while(valid_pt(al,shape) and data[tuple(al)]<thres):
            stop = stop+1
            al[dim] = stop
        
        seg_slice.append(slice(start+1,stop))

    return seg_slice
