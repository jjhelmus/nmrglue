"""
Peak picking routines, lineshape parameter guessing, and related functions.

"""

import numpy as np
import scipy

# nmrglue analysis module imports
from analysisbase import ndwindow_index,squish,find_limits,limits2slice
from analysisbase import ls_str2class,pick2linesh,linesh2pick,slice2limits

# lineshape classes
from analysisbase import peak1D,gauss1D,lorentz1D,scale1D

from segmentation import find_downward,find_all_downward
from segmentation import find_upward,find_all_upward
from segmentation import find_connected,find_all_connected
from segmentation import find_nconnected,find_all_nconnected

# Top level peak picking routine

def pick(data,thres,msep=None,direction='both',algorithm="thres",
         est_params=True,lineshapes=None):
    """
    Pick (find) peaks in a spectrum using a given algorithm.  

    Parameters:

    * data          N-dimensional array.
    * thres         Threshold value for minimum peak height.
    * msep          N-tuple of minimum peak seperations along each axis.
                    Must be defined if algorithm is 'thresh' or 'thresh-fast'
    * direction     Direction of peaks, 'positive','negative', or 'both'
                    or short cut values 'p','n','b'.
    * algorithm     Peak picking algorithm to use.  Options are 'thres',
                    'thres-fast', 'downward', or 'connected'
    * est_params    Set to True to perform a rough estimate of linewidths and
                    amplitude for all peaks picked.  False returns only the
                    peak locations (centers)
    * lineshapes    A list of lineshape classes or string shortcuts for each 
                    dimension.  If not specified Gaussian type lineshapes with 
                    a FWHM  linewidth parameter is assumed in each dimension.  
                    This parameter if only used if est_params is True.

    Returns:    centers,[linewidths,amplitudes]

    * centers       Array of estimated peak locations, shape (n_peaks,ndim).
    * linewidths    Array of estimated peak linewidths, shape (n_peaks,ndim).
    * amplitudes    Array of estimated peak amplitude, shape (n_peaks).

    """
    # parameter checking
    ndim = len(data.shape)

    # check algorithm
    if algorithm not in ['thres','thres-fast','downward','connected']:
        raise ValueError('Invalid algorithm %s'%(algorithm))
    
    # replace direction with shortcut value
    direction = direction.lower()

    if direction == 'both':
        direction = 'b'
    if direction == 'positive':
        direction = 'p'
    if direction == 'negative':
        direction = 'n'

    # check direction
    if direction not in ['p','n','b']:
        raise ValueError("invalid dir: %s"%(direction))

    # check msep
    if type(msep) == int:
        msep = (msep,)
    if algorithm in ['thres','thres-fast'] and len(msep) != ndim:
        raise ValueError("msep has incorrect length")
    
    # check  lineshapes
    if est_params:
        # expand None
        if lineshapes == None:
            lineshapes = [peak1D() for i in range(ndim)]
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

    # find centers and segments if requested
    if algorithm == 'thres':
        centers = pick_thres(data,thres,msep,direction)
    if algorithm == 'thres-fast':
        centers = pick_thres_fast(data,thres,msep,direction)
    if algorithm == 'downward':
        if est_params:
            centers,segments=pick_downward(data,thres,direction,seg_flag=True)
        else:
            centers = pick_downward(data,thres,direction,seg_flag=False)
    if algorithm == 'connected':
        if est_params:
            centers,segments=pick_connected(data,thres,direction,seg_flag=True)
        else:
            centers = pick_connected(data,thres,direction,seg_flag=False)

    # if full parameter estimation not requests return just the centers
    if est_params==False:
        return centers
    
    # estimate the linewidths and amplitudes
    lw = np.zeros(centers.shape,dtype=float)
    amp = np.zeros(len(centers),dtype=float)

    if  algorithm in ['downward','connected']:   # we already have segments
        for i,seg in enumerate(segments):
            null,lw[i],amp[i]=guess_params_segment(data,seg,ls_classes)
    
    else:   # no segments, loop over the centers
        for i,c in enumerate(centers):
            if c.shape == ():
                center = (c,)
            else:
                center = tuple(c)

            null,lw[i],amp[i]=guess_params_center(data,center,thres,ls_classes)

    return centers,lw,amp


# region finding and formatting

def find_regions(data,centers,thres,linewidths=None,amplitudes=None):
    """
    Find regions in spectra containing peaks.

    Parameters:

    * data          N-dimensional array.
    * centers       Array of peak centers.
    * thres         Threshold value for segmenting.
    * linewidths    Array of peak linewidths, optional.
    * amplitudes    Array of peak amplitudes, optional.

    Returns: regions

    * regions   List of dictionaries defining peak containing regions in the 
                spectra.  Each region dictionary has the following keys:

                * 'min'         Tuple of the minimum corner of the region.
                * 'max'         Tuple of the maximum corner of the region.
                * 'centers'     List of peak centers in the region.
                * 'linewidths'  List of peak linewidths in the region.
                * 'amplitudes'  List of peak amplitudes in the region.
    
                The 'linewidths' and 'amplitude' keys are only created
                if linewidths and amplitudes are passed to the function.

    """
    lcen = list(centers)    # list of peak centers
    
    if linewidths!=None:
        if len(linewidths) != len(lcen):
            s="Number of linewidths does not match the number of peaks"
            raise ValueError(s)
        llwd = list(linewidths) # list of linewidths
        
    if amplitudes!=None:
        if len(amplitudes)!= len(lcen):
            s="Number of amplitudes does not match the number of peaks"
            raise ValueError(s)
        lamp = list(amplitudes) # list of amplitudes

    regions = []    # list of spectral regions

    while len(lcen)!=0: # loop until all peaks have been accounted for
        
        # pop off a peak center
        pt = tuple(lcen.pop(0))
        
        # find the connected segment containing the peak
        if data[pt] > 0:
            seg = find_connected(data,pt,thres)
        else:
            seg = find_nconnected(data,pt,-thres)

        # find limits surrounding the connected segment
        limits = find_limits(seg)
        
        # create the current region dictionary
        r = {'min':tuple(limits[0]),
             'max':tuple(limits[1]),
             'centers':[pt]}
        
        # find other peaks in region and add (from back to front)
        to_add = pts_in_limits(lcen,limits)[::-1]
        for new_pt in to_add:
            r['centers'].append( tuple(lcen.pop(new_pt)) )

        # add linewidths if needed.
        if linewidths!=None:
            lw = tuple(llwd.pop(0))
            r['linewidths'] = [lw]
        
            for new_pt in to_add:
                r['linewidths'].append( tuple(llwd.pop(new_pt)) )

        # add amplitudes if needed.
        if amplitudes!=None:
            amp = lamp.pop(0)
            r['amplitudes'] = [amp]

            for new_pt in to_add:
                r['amplitudes'].append( lamp.pop(new_pt))

        # add the current region dictionary to the list of regions
        regions.append(r)
    
    return regions


def region2linesh(region):
    """
    Convert a region dictionary to linesh input

    Parameters:

    * region    Region dictionary which has linewidth and amplitude keys

    Returns: (rslice,min,guesses,amp_guesses)

    * rslice        Slice objects which will slice the array to give a region
    * min           List of mimimum to add to the output
    * guesses       linesh.fit_NDregion guesses input
    * amp_guesses   linesh.fit_NDregion amp_guesses input

    """
    
    if 'amplitudes' not in region:
        raise ValueError("region must contain amplitude estimates")
    if 'linewidths' not in region:
        raise ValueError("region must contain linewidths estimates")

    # create the slice object
    rslice = limits2slice( (region['min'],region['max']) )

    # create the peak amplitudes
    amp_guesses = region['amplitudes']

    # build the parameters list
    min = region['min']
    guesses = []
    for cen,lws in zip(region['centers'],region['linewidths']):
        guesses.append([(c-e,l) for c,l,e in zip(cen,lws,min)])

    return rslice,guesses,amp_guesses




def linesh2region(params_best,amp_best,rslice):
    """
    Convert linesh.fit_NDregion output to a region dictionary

    Parameters:
    
    * params_best   fit_NDregion output
    * amp_best      fit_NDregion output
    * min           List of minimums to add to the output

    Returns: region dictionary

    """

    
    limits = slice2limits(rslice)
    mins = limits[0]

    # tuples for first parameters pluse minimum for each dimension 
    # for each peak in the parameter list.
    cns = [ tuple([d[0]+m for (d,m) in zip(p,mins)]) for p in params_best] 

    # similar list comprehesions 
    lws = [ tuple([d[1] for d in p]) for p in params_best ]

    # create the region dictionary
    return {'min':tuple(limits[0]),
            'max':tuple(limits[1]),
            'centers': cns,
            'linewidths': lws,
            'amplitudes':list(amp_best) }
        

def filter_by_distance(data,centers,msep,lineshapes=None,amplitudes=None):
    """
    Filter peaks which are nearby, keeping those with the largest intesity.

    Parameters:
        
        * data          N-dimensional array.     
        * centers       Array of estimated peak locations.
        * msep          N-tuple of minimum peak seperations along each axis.
        * linewidths    Array of estimated peak linewidths, optional.
        * amplitudes    Array of estimated peak amplitude, optional.

    Returns: centers,[linewidths,amplitudes]

    """

    keep = []
    lcen = list(centers)

    # loop over all peaks, this is inefficient as some points have already been
    # checked in previous iterations, but it works and keeps the indexing nice.
    for pt in centers:

        # limits for seperation region
        min = [(i-s) for i,s in zip(pt,msep)]
        max = [(i+s) for i,s in zip(pt,msep)]
    
        # find all peak centers within region
        cen_idx = pts_in_limits(centers,(min,max))
        
        # find the intensity (absolute value) at each peak center in region
        data_idx = np.take(centers,cen_idx,axis=0)
        vals = [np.abs(data[tuple(i)]) for i in data_idx]
        
        # the index of the major peak center has the largest intensity
        major_idx = cen_idx[np.argmax(vals)]

        keep.append(major_idx)   # add the major peak index to keep list
        
    # remove duplicates in keep
    keep = list(set(keep))

    # return the filtered results.
    if lineshapes==None and amplitudes==None:
        return centers.take(keep,axis=0)
    if amplitudes==None:
        return centers.take(keep,axis=0),lineshapes.take(keep,axis=0)
    if lineshapes==None:
        return centers.take(keep,axis=0),amplitudes.take(keep)
    
    return ( centers.take(keep,axis=0), lineshapes.take(keep,axis=0),
             amplitudes.take(keep) )


def pts_in_limits(pts,lms):
    """
    Find all points in pts that are within a box defined by limits lms

    Returns a list of indicies of pts which are within the box limits.
    """
    return [i for i,pt in enumerate(pts) if in_limits(pt,lms)]
    

def in_limits(pt,lms):
    """
    Return True/False depending on if point (pt) is in limits (lms).
    """
    return (False not in [min<=x<=max for x,min,max in zip(pt,lms[0],lms[1])])



# algorithm specific peak picking routines

def pick_connected(data,thres,direction='both',seg_flag=False):
    """
    Peak pick a spectrum using the connected path algorithm.

    Find peaks (local maxima/minima) in an arbitrary dimensional NMR spectra
    above/below a given threshold and not part of an already defined
    connected segment.

    Parameters:

    * data      N-dimensional array.
    * thres     Threshold value for minimum peak height.
    * direction Direction of peaks, 'positive','negative', or 'both'
                or short cut values 'p','n','b'.
    * seg_flag  Set True to return list of points in each segment, False to 
                return only peak centers.

    Return: centers,[segments]

    * centers   array of peak locations with shape (n_peaks,ndim).
    * segments  List of  all points in a given segment, optional with the 
                seg_flag.

    """
    # replace with shortcut values
    direction = direction.lower()

    if direction == 'both':
        direction = 'b'
    if direction == 'positive':
        direction = 'p'
    if direction == 'negative':
        direction = 'n'

    if direction not in ['p','n','b']:
        raise ValueError("invalid dir: %s"%(direction))

    # positive peaks
    if direction != 'n':
        pcenters,psegments = find_all_connected(data,thres)
    if direction == 'p':
        if seg_flag:
            return np.array(pcenters),psegments
        else:
            return np.array(pcenters)

    # negative peaks
    if direction !='p':
        ncenters,nsegments = find_all_nconnected(data,-thres)
    if direction == 'n':
        if seg_flag:
            return np.array(ncenters),nsegments
        else:
            return np.array(ncenters)
    
    # return both peaks
    if seg_flag:
        return np.array(pcenters+ncenters),psegments+nsegments
    else:
        return np.array(pcenters+ncenters)


def pick_downward(data,thres,direction='both',seg_flag=False):
    """
    Peak pick a spectrum using a downward/upward path algorithm.

    Find peaks (local maxima/minima) in an arbitrary dimensional NMR spectra
    above/below a given threshold and not part of an already defined
    downward/upward segment.

    Parameters:

    * data      N-dimensional array.
    * thres     Threshold value for minimum peak height.
    * direction Direction of peaks, 'positive','negative', or 'both'
                or short cut values 'p','n','b'.
    * seg_flag  Set True to return list of points in each segment, False to 
                return only peak centers.

    Return: centers,[segments]

    * centers   array of peak locations with shape (n_peaks,ndim).
    * segments  List of  all points in a given segment, optional with the 
                seg_flag.

    """
    # replace with shortcut values
    direction = direction.lower()

    if direction == 'both':
        direction = 'b'
    if direction == 'positive':
        direction = 'p'
    if direction == 'negative':
        direction = 'n'

    if direction not in ['p','n','b']:
        raise ValueError("invalid dir: %s"%(direction))

    # positive peaks
    if direction != 'n':
        pcenters,psegments = find_all_downward(data,thres)
    if direction == 'p':
        if seg_flag:
            return np.array(pcenters),psegments
        else:
            return np.array(pcenters)

    # negative peaks
    if direction !='p':
        ncenters,nsegments = find_all_upward(data,-thres)
    if direction == 'n':
        if seg_flag:
            return np.array(ncenters),nsegments
        else:
            return np.array(ncenters)
    
    # return both peaks
    if seg_flag:
        return np.array(pcenters+ncenters),psegments+nsegments
    else:
        return np.array(pcenters+ncenters)


def pick_thres(data,thres,msep,direction='both'):
    """
    Peak pick a spectrum using a threshhold-minimum distance algorithm.
    
    Find peaks (local maxima/minima) in a arbitrary dimensional NMR spectra 
    above/below a set threshold with a minimal distance between peaks.  When 
    the spectrum is small and multiple copies can fit into RAM use the _fast 
    version of this function.

    Parameters:

    * data      N-dimensional array.
    * thres     Threshold value for minimum peak height
    * msep      N-tuple of minimum peak seperations along each axis
    * direction Direction of peaks, 'positive','negative', or 'both'
                or short cut values 'p','n','b'.

    Returns: array of peak locations with shape (n_peaks,ndim)

    """
    # replace with shortcut values
    direction = direction.lower()

    if direction == 'both':
        direction = 'b'
    if direction == 'positive':
        direction = 'p'
    if direction == 'negative':
        direction = 'n'

    if direction not in ['p','n','b']:
        raise ValueError("invalid dir: %s"%(direction))

    peaks = []  # create an empty list of peaks 
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    middle = np.floor((np.array(msep)-1)/2.)    # index of middle of window
    ms = [slice(x,x+1) for x in middle]         # middle slice list

    # loop over the windows
    for idx,s in ndwindow_index(data.shape,wsize):
        window = data[s]
        max = window.max()
        min = window.min()
        #print idx
        if max == data[idx] and max > thres and direction != 'n':
            peaks.append(idx)
        if min == data[idx] and min < -thres and direction != 'p':
            peaks.append(idx)

    return np.array(peaks)


def pick_thres_fast(data,thres,msep,direction="both"):
    """
    Fast version of pick_thres functions

    See peak_thres for call details

    """
    # replace with shortcut values
    direction = direction.lower()

    if direction == 'both':
        direction = 'b'
    if direction == 'positive':
        direction = 'p'
    if direction == 'negative':
        direction = 'n'

    if direction not in ['p','n','b']:
        raise ValueError("invalid dir: %s"%(direction))

    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    if direction != 'n':
        # find local maxima mask
        mx=scipy.ndimage.maximum_filter(data,size=wsize,mode='constant')==data
   
        # find positive threshold mask
        pthres = np.ma.greater(data,thres)
    
    if direction != 'p':
        # find the local minimum mask
        mn=scipy.ndimage.minimum_filter(data,size=wsize,mode='constant')==data

        # find negative threshold mask
        nthres = np.ma.less(data,-thres)

    if direction == 'p':
        # peaks are bitwise and of maximum mask and threshold mask
        return np.transpose( np.nonzero( np.bitwise_and(pthres,mx) ) )
    if direction == 'n':
        # peaks are bitwise and of minimum mask and threshold mask
        return np.transpose( np.nonzero( np.bitwise_and(nthres,mn) ) )

    # peaks are bitwise or of above lists
    return np.transpose( np.nonzero( np.bitwise_or(
        np.bitwise_and(pthres,mx),
        np.bitwise_and(nthres,mn) ) ) )


# parameter guessing functions

def guess_params_center(data,center,thres,lineshapes):
    """
    Guess the parameter of a peak centered at center.

    Parameters:

    * data          Spectral data
    * center        Location of peak center
    * thres         Noise threshold 
    * lineshapes    List of lineshape classes

    Return: centers,linewidths,amplitudes

    * centers   Array of estimated peak centers in each dimension.
    * linewidth Array of estimated linewidths in each dimension
    * amplitude Estimate of peak amplitude

    """

    # find the segment containing the peak
    if data[center] > 0:
        segment = find_downward(data,center,thres)
    else:
        segment = find_upward(data,center,-thres)

    return guess_params_segment(data,segment,lineshapes)

def guess_params_segment(data,segment,lineshapes):
    """
    Guess the parameter of a peak in a provided segment.

    Parameters:

    * data          Spectral data
    * segment       List of points in segment.
    * lineshapes    List of lineshape classes

    Return: centers,linewidths,amplitudes

    * centers   Array of estimated peak centers in each dimension.
    * linewidth Array of estimated linewidths in each dimension
    * amplitude Estimate of peak amplitude

    """

    # find the rectangular region around the segment
    limits = find_limits(segment)
    region= data[limits2slice(limits)]

    # amptide is estimated by the sum of all points in region
    amp = np.sum(region)

    lw  = []    # list of linewidths
    cen = []    # list of peak centers
    # loop over the axes
    for axis,ls in enumerate(lineshapes):
        # create the 1D lineshape 
        r = squish(np.copy(region),axis)
        # estimate the linewidth
        center,linewidth = ls.guessp(r)
        lw.append(linewidth)
        cen.append(center)

    return np.array(cen),np.array(lw),amp
