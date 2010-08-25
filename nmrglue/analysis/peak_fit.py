"""
peak_fit - Module for fitting peaks in NMR spectra
"""

import numpy as np
import scipy.optimize
pi = np.pi

# 1D Lineshapes

def gauss1D(M,mu,sigma):
    """
    Generate a 1D Gaussian window of length M (0...M-1)

    Parameters:
    -----------

    * M     Length of window (0,1,...,M-1)
    * mu    Location of mean (center)
    * sigma Varian (width) of distribution.

    Returns: 1D array of length M

    """
    return np.exp(-(np.arange(M)-mu)**2/(2*sigma**2))/(np.sqrt(2*pi*sigma**2))


def lorentz1D(M,x_0,gamma):
    """
    Generate a 1D Lorentz window of length M (0...M-1)

    Parameters:
    -----------

    * M     Length of window (0,1,...,M-1)
    * x_0   Location of center of distribution.
    * gamma Scale (width) of distribution.

    Returns: 1D array of length M

    """

    return 1./pi*1./(gamma**2 + (np.arange(M)-x_0)**2)

def gen_lineshape(M,params,shape,A=1.):
    """
    Generate a 1D linshape of length M

    Parameters:
    -----------
    
    * M         Length of lineshape (0,1,...,M-1).
    * params    Lineshape parameters (shape specific).
    * shape     Shape of lineshape.
    * A         Amplitude of lineshape.

    Returns:    1D array of length M

    shape can be one of: 'gauss','lorentz'

    """
    
    if shape == "gauss":
        mu,sigma = params
        return A*gauss1D(M,mu,sigma)
    elif shape == "lorentz":
        x_0,gamma = params
        return A*lorentz1D(M,x_0,gamma)
    else:
        raise ValueError("Unknown shape")

###################################
# Parameter guessing 1D functions #
###################################

def centerFWHM2params_1D(center,FWHM,shape):
    """
    Convert a center and FWHM into lineshape parameters

    Parameters:
    -----------

    * center    Location of center of peak.
    * FWHM      Full width half max of peak.
    * shape     Shape of lineshape

    Returns:    parameter tuple according to shape

    """

    if shape == 'gauss':
        # mu is the center, FWHM = 2*sqrt( 2*ln(2) )*sigma = 2.3548*sigma
        return (center,FWHM/2.3548)
    elif shape == 'lorentz':
        # center is x_0, FWHM is gamma
        return (center,FWHM)
    else:
        raise ValueError("Unknown shape")
    

def guess_FWHM_1D(line):
    """
    Guess FWHM for a line returning the full width half max (FWHM)
    """
    
    # find the points that bracket the first and last crossing of the 
    # half max, from each use a linear extrapolation to find the location
    # of the half max on either side of the maximum.  The difference between
    # these two values is the full width.
    max = line.max()
    hmax = max/2.
    
    top_args = np.nonzero(line > hmax)

    l_idx  = top_args[0][0]      #  
    r_idx = top_args[0][-1]  

    # solve hmax = mx+b => x = y-b/m
    # for two points x_0 and x_1 this becomes y = (hmax-x_0)/(x_1-x_0)
    # to this value add the index of x_0 to get the location of the half-max.

    # the left side of the line
    x_0,x_1 = line[l_idx-1],line[l_idx]
    left = l_idx-1+(hmax-x_0)/(x_1-x_0)

    # the right side of the line
    x_0,x_1 = line[r_idx],line[r_idx+1]
    right = r_idx+(hmax-x_0)/(x_1-x_0)

    return right-left

########################
# 1D fitting functions #
########################

def gen_peak1D(x_len,x_params,x_shape,A=1.):
    """
    Make a 1D peak:
    
    Parameter:
    ----------
    * x_len     Length of peak region along X (0) axis.
    * x_params  Lineshape parameters of X (0) axis (shape specific).
    * x_shape   Shape of lineshape along X (0) axis.
    * A         Amplitude of peak.

    Returns: 1D array of shape (x_len)

    """    
    return A*gen_lineshape(x_len,x_params,x_shape)

def err_peak1D(p,peak,x_len,x_shape):
    """
    Error function for 1D peak, to be called by fmin functions.
    
    Parameters:
    -----------
    * p         List of simulation parameters (A,x_params)
    * peak      Peak region to fit against.
    * x_len     Length of peak region along X (0) axis.
    * x_shape   Shape of lineshape along X (0) axis.

    Returns:    Sum of squares between peak and simulated peak

    """
    
    # split the parameter list
    A = p[0]
    x_params = p[1:]

    # simulate the peak
    sim_peak = gen_peak1D(x_len,x_params,x_shape,A)

    # calculate the error
    return ((peak-sim_peak)**2).sum()

def fit_region1D_simple(peak,x_shape,x_guess=None,A_guess=None):
    """
    Fit a 1D region to a single peak

    Parameters:
    -----------

    * peak      Peak region to fit against.
    * x_shape   Shape of lineshape along X (0) axis.
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).
    * A_guess   Guess of peak amplitude.

    If any of x_guess or A_guess are not defined (None) they will be determined
    automatically.

    Returns:
    --------

    * A_best    Peak amplitude of best fit to data
    * x_best    Best lineshape parameters of X (0) axis.

    """
    
    # automatically determind parameter is needed
    if x_guess == None or A_guess==None:
        A_g,x_g = guess_params_1D(peak,x_shape)
    if x_guess == None:
        x_guess = x_g 
    if A_guess == None:
        A_guess = A_g
    
    # prepare for downhill simplex minimization
    x_len = peak.shape[0]
    p0 = tuple([A_guess])+tuple(x_guess)
    args = (peak,x_len,x_shape)

    # perform minimization
    p_best = scipy.optimize.fmin(err_peak1D,p0,args=args)

    return p_best[0],p_best[1:]


def guess_params_1D(peak,x_shape):
    """
    Guess parameters of a 1D peak

    Parameters:
    -----------

    * peak      Peak region to fit against.
    * xshape    Shape of lineshape along X (0) axis.

    Returns: (A_guess,x_guess,y_guess)
    --------

    * A_guess   Amplitude guess
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).

    """
    
    # assuming the region contained in the peak array the total integrated
    # area is the amplitude of the peak.
    A_guess = peak.sum()

    # the maximum is a good guess at the center of the peak
    x_0,y_0 = peak.argmax()

    # guess FWHM from 1D trace along through maximum
    x_FWHM = guess_FWHM_1D(peak)

    # make the parameter tuples for each dimension
    x_guess = centerFWHM2params_1D(x_0,x_FWHM,x_shape)

    return A_guess,x_guess


########################
# 2D fitting functions #
########################

def gen_peak2D(x_len,y_len,x_params,y_params,x_shape,y_shape,A=1.):
    """ 
    Make a 2D peak.

    Parameter:
    ----------
    
    * x_len     Length of peak region along X (0) axis. 
    * y_len     Length of peak region along Y (1) axis.
    * x_params  Lineshape parameters of X (0) axis (shape specific).
    * y_params  Lineshape parameters of Y (1) axis (shape specific).
    * x_shape   Shape of lineshape along X (0) axis.
    * y_shape   Shape of lineshape along Y (1) axis.
    * A         Amplitude of peak.

    Returns:    2D array of shape (x_len,y_len)

    See gen_lineshape for allowed x_shape and y_shape values.

    """
    x_lineshape = gen_lineshape(x_len,x_params,x_shape)
    y_lineshape = gen_lineshape(y_len,y_params,y_shape)
    return A*np.outer(x_lineshape,y_lineshape)

def err_peak2D(p,peak,x_len,y_len,x_shape,y_shape,nx_params,ny_params):
    """
    Error function for 2D peak, to be called by fmin functions.

    Parameters:
    -----------

    * p         List of simulation parameters (A,x_params,y_params)
    * peak      Peak region to fit against.
    * x_len     Length of peak region along X (0) axis.
    * y_len     Length of peak region along Y (1) axis.
    * x_shape   Shape of lineshape along X (0) axis.
    * y_shape   Shape of lineshape along Y (1) axis.
    * nx_params Number of x parameters in p.
    * ny_params Number of y parameters in p.

    Returns:    Sum of squares between peak and simulated peak

    """
    
    # split the the parameter list to A, x and y parameters
    A = p[0]
    x_params = p[1:nx_params+1]
    y_params = p[-ny_params:]

    # simulate the peak
    sim_peak = gen_peak2D(x_len,y_len,x_params,y_params,x_shape,y_shape,A)

    # calculate the error
    return ((peak-sim_peak)**2).sum()


def fit_region2D_simple(peak,x_shape,y_shape,x_guess=None,y_guess=None,
                        A_guess=None):
    """
    Fit a 2D region to a single peak

    Parameters:
    -----------

    * peak      Peak region to fit against.
    * x_shape   Shape of lineshape along X (0) axis.
    * y_shape   Shape of lineshape along Y (1) axis.
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).
    * y_guess   Guess of lineshape parameters of Y (1) axis (shape specific).
    * A_guess   Guess of peak amplitude.

    If any of x_guess, y_guess or A_guess are not defined (None) they will be 
    determined automatically.

    Returns:
    --------

    * A_best    Peak amplitude of best fit to data
    * x_best    Best lineshape parameters of X (0) axis.
    * y_best    Best lineshape parameters of Y (1) axis.

    """
    
    # automatically determind parameter is needed
    if x_guess == None or y_guess == None or A_guess==None:
        A_g,x_g,y_g = guess_params_2D(peak,x_shape,y_shape)
    if x_guess == None:
        x_guess = x_g 
    if y_guess == None:
        y_guess = y_g
    if A_guess == None:
        A_guess = A_g
    

    # prepare for downhill simplex minimization
    x_len = peak.shape[0]
    y_len = peak.shape[1]
    p0 = tuple([A_guess])+tuple(x_guess)+tuple(y_guess)
    nx_params = len(x_guess)
    ny_params = len(y_guess)
    args = (peak,x_len,y_len,x_shape,y_shape,nx_params,ny_params)

    # perform minimization
    p_best = scipy.optimize.fmin(err_peak2D,p0,args=args)

    return p_best[0],p_best[1:nx_params+1],p_best[-ny_params:]


def guess_params_2D(peak,x_shape,y_shape):
    """
    Guess parameters of a 2D peak

    Parameters:
    -----------

    * peak      Peak region to fit against.
    * xshape    Shape of lineshape along X (0) axis.
    * yshape    Shape of lineshape along Y (1) axis.

    Returns: (A_guess,x_guess,y_guess)
    --------

    * A_guess   Amplitude guess
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).
    * y_guess   Guess of lineshape parameters of Y (1) axis (shape specific).

    """
    
    # assuming the region contained in the peak array the total integrated
    # area is the amplitude of the peak.
    A_guess = peak.sum()

    # the maximum is a good guess at the center of the peak
    x_0,y_0 = np.unravel_index(peak.argmax(),peak.shape)

    # guess FWHM from 1D trace along through maximum
    x_FWHM = guess_FWHM_1D(peak[:,y_0].flatten())
    y_FWHM = guess_FWHM_1D(peak[x_0,:].flatten())

    # make the parameter tuples for each dimension
    x_guess = centerFWHM2params_1D(x_0,x_FWHM,x_shape)
    y_guess = centerFWHM2params_1D(y_0,y_FWHM,y_shape)

    return A_guess,x_guess,y_guess


