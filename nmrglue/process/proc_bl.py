"""
A collection of NMR processing functions for filtering, smoothing, and 
correcting spectral baselines.  
"""

import numpy as np
import scipy
import scipy.ndimage

# Linear (First order) baseline correction

def base(data, nl, nw=0):
    """
    Linear (first-order) baseline correction based on node list.

    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    nl : list
        List of baseline nodes.
    nw : float, optional
        Node half-width in points.

    Returns
    -------
    ndata : ndarray
        NMR data with first order baseline correction appied.  For 2D data
        baseline correction is applied for each trace along the last 
        dimension.

    """
    if data.ndim == 1:
        data = data - calc_bl_linear(data, nl, nw)
    else:   # for 2D array loop over traces
        for i, vec in enumerate(data):
            data[i] = data[i] - calc_bl_linear(vec, nl, nw)
    return data

def calc_bl_linear(x, nl, nw=0):
    """ 
    Calculate a baseline using linear approximation between nodes

    Parameters
    ----------
    x : 1D ndarray
        One-dimensional NMR data.
    nl : list
        List of baseline nodes
    nw : float   
        Node half-width in points
    
    Returns
    -------
    baseline : ndarray
        Base calculated using linear approximation between nodes.

    """
    bl = np.zeros_like(x)
    for i in range(len(nl) - 1):
        # minimum and maximum index
        min = nl[i]
        max = nl[i + 1]
        # linspace s1 and s2
        s1 = x[min - nw:min + nw + 1].mean()
        s2 = x[max - nw:max + nw + 1].mean()
        bl[min:max + 1] = np.linspace(s1, s2, max - min + 1)
    return bl

# Constant baseline correction

def cbf(data, last=10, apply=slice(None)):
    """
    Constant baseline correction from percentage of data.

    Parameters
    ----------
    data : 1D or 2D ndarray 
        Array of 1D or 2D NMR data.
    last : float, optional
        Percent (0-100) of last axis used to calculate the baseline 
        correction.
    apply : slice, optional
        Slice describing first-axis region(s) to which the baseline correction
        should be applied.  Parameter is ignored for 1D data.

    Returns
    -------
    ndata : 1D or 2D ndarray
        NMR data with a constant baseline subtracted.

    """
    # calculate the correction
    n = data.shape[-1] * last / 100. + 1.  
    corr = data[..., -n:].sum(axis=-1) / n
    
    # apply correction
    if data.ndim == 2:
        data[apply] =  data[apply] - np.array([corr]).transpose()[apply]
        return data
    else:
        return data - corr

def cbf_explicit(data, calc=slice(None), apply=slice(None)):
    """
    Constant Baseline correction from an explicit region of the data.

    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    calc : slice, optional
        Slice describing region to use for calculating the baseline correction.
    apply : slice, optional
        Slice describing first-axis region(s) to which the baseline correction
        should be applied.  Parameter is ignored for 1D data.

    Returns
    -------
    ndata : 1D or 2D ndarray
        NMR data with a constant baseline subtracted.

    """
    # calculate correction
    n = len(range(data.shape[-1])[calc])
    corr = data[..., calc].sum(axis=-1) / n

    # apply correction
    if data.ndim == 2:
        data[apply] = data[apply] - np.array([corr]).transpose()[apply]
        return data
    else:
        return data - corr

# Median baseline correction

def med(data, mw=24, sf=16, sigma=5.0):
    """
    Median baseline correction.

    Algorith described in: Friedrichs, M.S. JBNMR 1995 5 147-153.

    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    mw : float
        Median window size in pts.
    sf : float
        Smooth window size in pts.
    sigma : float
        Standard-deviation of convoluted Gaussian window.

    Returns
    -------
    ndata : 1D or 2D ndarray
        NMR data with the median baseline subtracted.

    """
    if data.ndim == 1:
        data = data - calc_bl_med(data, mw, sf, sigma)
    else:
        for i, vec in enumerate(data):
            data[i] = vec - calc_bl_med(vec, mw, sf, sigma)
    return data

def calc_bl_med(x, mw, sf, sigma):
    """
    Calculate a baseline using median baseline correction.

    Algorithm described in: Friedrichs, M.S. JBNMR 1995 5 147-153

    Parameters
    ----------
    x : 1D ndarray
        One dimensional NMR data
    mw : float
        Median window size in pts.
    sf : float
        Smooth window size in pts.
    sigma : float
        Standard-deviation of convoluted Gaussian window.

    Returns
    -------
    baseline : 1D ndarray
        Baseline calculated using median baseline correction

    """
    # create extrema array (non extrema values are masked out)
    mask = x == scipy.ndimage.median_filter(x, size=3)
    mask[0] = False     # first pt always extrema
    mask[-1] = False    # last pt always extrema
    e = np.ma.masked_array(x, mask)

    # fill in the median vector
    m = scipy.ndimage.median_filter(e, mw + 1, mode="mirror")
    # using the median_filter might give slightly different results than
    # described algorithm but is MUCH faster

    # convolve with a gaussian
    g = scipy.signal.gaussian(sf, sigma)
    g = g / g.sum()
    return scipy.signal.convolve(m, g, mode='same')

# Solvent Filter

def sol_general(data, filter, w=16, mode='same'):
    """
    Solvent filter with generic filter.

    Algorithm described in: Marion et al. JMR 1989 84 425-430
    
    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    filter : ndarray
        Filter to convolve with data.  Not used in solvent filter functions
        which specific the filter, e.g. sol_boxcar.
    w : int, optional
        Filter length.  Not used here but is used in solent filter functions 
        which specificy the filter, e.g. sol_boxcar.
    mode : {'valid', 'same', 'full'}, optional
        Convolution mode, 'same' should be used.

    Returns
    -------
    ndata : 1D or 2D ndarray
        NMR data with solvent filter applied

    """
    A = filter.sum()
    if data.ndim == 2:
        filter = np.atleast_2d(filter)
    return data - scipy.signal.convolve(data, filter, mode=mode) / A

def sol_boxcar(data, w=16, mode='same'):
    """ 
    Solvent filter with boxcar filter.  See :py:func:`sol_general`.
    """
    filter = scipy.signal.boxcar(w)
    return sol_general(data, filter, w=w, mode=mode)

def sol_sine(data, w=16, mode='same'):
    """ 
    Solvent filter with sine-bell filter. See :py:func:`sol_general`.
    """
    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, w))
    return sol_general(data, filter, w=w, mode=mode)

def sol_sine2(data, w=16, mode='same'):
    """ 
    Solvent filter with square sine-bell filter. See :py:func:`sol_general`.
    """
    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, w)) ** 2
    return sol_general(data, filter, w=w, mode=mode)

def sol_gaussian(data, w=16, mode='same'):
    """ 
    Solvent filter with square gaussian filter. See :py:func:`sol_general`.
    """
    filter = scipy.signal.gaussian(w, w / 2.)
    return sol_general(data, filter, w=w, mode=mode)

# Polynomial Solvent Subtraction

def poly_td(data):
    """ 
    Polynomial time domain solvent subtraction - Not Implemented
    """
    # From NMRPipe paper(appendix):
    #
    # when used with the argument -time, fits all data points to a polynomial,
    # which is then subtracted from the original data.  It is intended to fit
    # and subtract low-freqency solvent signal in the FID, a procedure that 
    # often causes less distortions than time-domain convolution methods.
    # By default, a fourth-order polynomials is used.  For speed successive
    # averages of regions are usually fit, rather than fitting all of the data.
    #
    # Algorithm
    #
    # 1. Calculate average of blocks
    # 2. Fit these averages to polynomial (block parameters)
    # 3. Back out "real" polynomial parameters from these block parameters
    # 4. Subtract off the polynomial from data
    raise NotImplementedError

def poly_fd(data):
    """ 
    Polynomial frequency domain baseline correction - Not Implemented
    """
    # From NMRPipe paper (appendix):
    #
    # applies a polynomial baseline correction of the order specified by 
    # argument -ord via an automated base-line detection method when used
    # with argument -auto.  The defauly is a forth-order polynomial. The 
    # automated base-line mode works as follows: a copy of a given vector is 
    # divided into a series of adjacent sections, typically eight points wide.
    # The average value of each section is subtracted from all points in that 
    # section, to generate a 'centered' vector.  The intensities of the entire
    # centered vector are sorted, and the standard deviation of the noise is
    # estimated under the assumption that a given fraction (typically about
    # 30%) of the smallest intensities belong to the base-line, and that the 
    # noise is normally distributed.  This noise estimate is multiplied by a 
    # constant, typically about 1.5, to yield a classification threshold.  
    # Then, each section in the centered vector is classified as base line only
    # if its standard deviation does not exceed the threshold.  These 
    # classifications are used to correct the original vector.
    #
    # Algorithm
    #
    # 1. Divide into 'blocks'
    # 2. Center each block and create centered vector
    # 3. Calculate intensity (abs) of each centered vector
    # 4. Sort intensities, lower 30% belong to baseline
    # 5. Fit base line intensities to Normal distribution, gives estimation
    #    of standard deviation (std) of noise
    # 6. Classification threshold set to 1.5*std
    # 7. Qualify each block in centered vector as baseline only 
    #    (its std < thres) or not (std > thres)
    # 8. Fit baseline only points to polynomial and substract off
    raise NotImplementedError
