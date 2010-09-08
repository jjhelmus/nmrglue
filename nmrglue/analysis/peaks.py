"""
Functions for peak picking, segmentation, and fitting.


Peak fitting
============

nmrglue has the ability to fit regions which contain a single peaks or multiple
peaks in both 1D and 2D data sets.  Regions with only one peak are fit using 
a fit_singleND function with correct dimensionality, Regions that do/may 
contain multiple peaks should be fit using the approiately dimensioned 
fit_multipleND function.


Spectrum and Region Segmentation
================================


Peak Picking
============


"""

import numpy as np
import numpy.ma as ma
import scipy
import scipy.optimize
pi = np.pi


################
# Peak fitting #
################


# to add new region to fit (like 3D+) the following functions should be 
#created where type is single/multiple/new_type:

# fit_typeND            High-level region fitting functions 
# unpack_typeNDp        Function to unpack fmin output
# err_typeND            fmin minimization functions
# sim_typeND            Simulate a region given parameters.
# guess_params_typeND   Determine initial optimization parameters.


# High-Level Fitting Functions (fit_typeND)

def fit_single1D(region,x_shape,x_guess=None,A_guess=None,
                 maxfun=None,maxiter=None,disp=False):
    """
    Fit a 1D region containing a single peak.  

    Parameters:
    -----------

    * region    1D spectral region to fit (array).
    * x_shape   Shape of lineshape along X (0) axis.
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).
    * A_guess   Guess of peak amplitude.

    Optional Optimization Parameters
    --------------------------------
    
    * maxfun    Maximum number of function evaluations.
    * maxiter   Maximum number of simplex iterations.
    * disp      Set to True to print convergence messages.


    If x_guess or A_guess are not defined (None) they will be determined
    automatically.

    Returns: (A_best,x_best)
    --------

    * A_best    Peak amplitude of best fit to data
    * x_best    Best lineshape parameters of X (0) axis.


    For a list of valid lineshape and corresponding parameter formats see
    gen_lineshape.  When the spectral region is known to, or may, contain 
    multiple peaks use fit_multiple1D.

    """

    # automatically determind parameter as needed
    if x_guess == None or A_guess==None:
        A_g,x_g = guess_params_single1D(region,x_shape)
    if x_guess == None:
        x_guess = x_g 
    if A_guess == None:
        A_guess = A_g
    
    # prepare for downhill simplex minimization
    x_len = region.shape[0]
    p0 = tuple([A_guess])+tuple(x_guess)
    args = (region,x_len,x_shape)

    # perform minimization
    p_best = scipy.optimize.fmin(err_single1D,p0,args=args,
                                maxfun=maxfun,maxiter=maxiter,disp=disp)

    return unpackp_single1D(p_best)


def fit_multiple1D(region,x_shapes=None,x_guesses=None,A_guesses=None,
                   thres=None,x_shape=None,n_peaks=None,
                   maxfun=None,maxiter=None,disp=False,):
    """
    Fit a 1D region which does/may contain multiple peaks.

    Parameters:
    -----------

    * region    1D spectral region to fit (array).
    * x_shapes  List of lineshapes along X(0) axis.
    * x_guesses List of lineshape parameters guesses of X (0) axis.
    * A_guesses List of peak amplitude guesses.
    * thres     Noise threshold (below this value data is considered noise)
    * x_shape   Shape of all lineshapes in region for use with thres.
    * n_peaks   Number of peaks in region used with thres (optional).

    Optional Optimization Parameters
    --------------------------------
    
    * maxfun    Maximum number of function evaluations.
    * maxiter   Maximum number of simplex iterations.
    * disp      Set to True to print convergence messages.

    Returns: (A_best,x_best)
    --------
    
    * A_best    Peak amplitude of best fit to data
    * x_best    Best lineshape parameters of X (0) axis.

    This functions operates in two modes depending on the parameters provided.

    Automatic mode:
    ---------------

    When thres and x_shape are both defined, the region will be segmented and
    starting parameters for peak fitting will be determined automatically.  In
    this mode n_peaks can optionally define the number of peak expected in the
    region.  If more peaks are found only the first n_peaks are fit, if fewer
    peaks are found an Exception is raised.

    Manual mode:
    ------------

    When thres or x_shape are not defined then initial peak parameters must
    be provided in x_shapes, x_guesses, and A_guesses.  thres, x_shape and 
    n_peaks are ignored in the fit.


    For a list of valid lineshape and corresponding parameter formats see
    gen_lineshape.  When the spectral region is known to contain a single peak
    use fit_single1D.

    """
   
    # Automatic mode
    if None not in [thres,x_shape]:
        # automatically find parameters
        t = guess_params_multiple1D(region,x_shape,thres,n_peaks)

        # unpack parameter that were found...
        A_guesses,x_guesses = t

        if n_peaks == None:
            n_peaks  = len(A_guesses)
        if len(A_guesses) != n_peaks:
            raise Exception("insufficient number of peaks found")
        x_shapes = [x_shape]*n_peaks
    
    # Manual mode
    if None in [x_shapes,x_guesses,A_guesses]:
        raise Exception("incomplete input")
    
    # check parameters
    n_peaks = len(A_guesses)
    if n_peaks!=len(x_guesses) or n_peaks!=len(x_shapes):
        raise ValueError("length of input parameters do not match")

    # prepare args for simplex minimization
    x_len = region.shape[0]
    nx_params = [len(i) for i in x_guesses]
    args = (region,x_len,x_shapes,nx_params)

    # prepare p0 for simplex
    p0 = []
    for i in A_guesses: 
        p0.append(i)
    for x_guess in x_guesses:
        for i in x_guess:
            p0.append(i)

    # performs minimization
    p_best = scipy.optimize.fmin(err_multiple1D,p0,args=args,
                    maxiter=maxiter,maxfun=maxfun,disp=disp)

    return unpackp_multiple1D(p_best,nx_params)


def fit_single2D(region,y_shape,x_shape,y_guess=None,x_guess=None,A_guess=None,
                        maxfun=None,maxiter=None,disp=False):

    """
    Fit a 2D region containing a single peak.  

    Parameters:
    -----------

    * region    2D spectral region to fit (array).
    * y_shape   Shape of lineshape along Y (0) axis.
    * x_shape   Shape of lineshape along X (1) axis.
    * y_guess   Guess of lineshape parameters of Y (0) axis (shape specific).
    * x_guess   Guess of lineshape parameters of X (1) axis (shape specific).
    * A_guess   Guess of peak amplitude.

    Optional Optimization Parameters
    --------------------------------
    
    * maxfun    Maximum number of function evaluations.
    * maxiter   Maximum number of simplex iterations.
    * disp      Set to True to print convergence messages.


    If x_guess, y_shape or A_guess are not defined (None) they will be 
    determined automatically.

    Returns: (A_best,y_best,x_best)
    --------

    * A_best    Peak amplitude of best fit to data
    * y_best    Best lineshape parameters of X (0) axis.
    * x_best    Best lineshape parameters of X (0) axis.


    For a list of valid lineshape and corresponding parameter formats see
    gen_lineshape.  When the spectral region is known to, or may, contain 
    multiple peaks use fit_multiple2D.

    """
    
    # automatically determind parameter is needed
    if x_guess == None or y_guess == None or A_guess==None:
        A_g,y_g,x_g = guess_params_single2D(region,y_shape,x_shape)
    if x_guess == None:
        x_guess = x_g 
    if y_guess == None:
        y_guess = y_g
    if A_guess == None:
        A_guess = A_g
    
    # prepare for downhill simplex minimization
    y_len = region.shape[0]
    x_len = region.shape[1]
    p0 = tuple([A_guess])+tuple(y_guess)+tuple(x_guess)
    nx_param = len(x_guess)
    ny_param = len(y_guess)
    args = (region,y_len,x_len,y_shape,x_shape,ny_param,nx_param)

    # perform minimization
    p_best = scipy.optimize.fmin(err_single2D,p0,args=args,
                                 maxfun=maxfun,maxiter=maxiter,disp=disp)
    
    return unpackp_single2D(p_best,ny_param,nx_param)


def fit_multiple2D(region,y_shapes=None,x_shapes=None,y_guesses=None,
                   x_guesses=None,A_guesses=None,thres=None,y_shape=None,
                   x_shape=None,n_peaks=None,
                   maxfun=None,maxiter=None,disp=False):
    """
    Fit a 2D region which does/may contain multiple peaks.

    Parameters:
    -----------

    * region    2D spectral region to fit (array).
    * y_shapes  List of lineshapes along Y (0) axis.
    * x_shapes  List of lineshapes along X (1) axis.
    * y_guesses List of lineshape parameters guesses of Y (0) axis.
    * x_guesses List of lineshape parameters guesses of X (0) axis.
    * A_guesses List of peak amplitude guesses.
    * thres     Noise threshold (below this value data is considered noise)
    * y_shape   Shape of all lineshapes in region for use with thres.
    * x_shape   Shape of all lineshapes in region for use with thres.
    * n_peaks   Number of peaks in region used with thres (optional).

    Optional Optimization Parameters
    --------------------------------
    
    * maxfun    Maximum number of function evaluations.
    * maxiter   Maximum number of simplex iterations.
    * disp      Set to True to print convergence messages.

    Returns: (A_best,x_best)
    --------
    
    * A_best    Peak amplitude of best fit to data
    * x_best    Best lineshape parameters of X (0) axis.

    This functions operates in two modes depending on the parameters provided.

    Automatic mode:
    ---------------

    When thres, y_shape and x_shape are defined, the region will be segmented 
    and starting parameters for peak fitting will be determined automatically.
    In this mode n_peaks can optionally define the number of peak expected in 
    the region.  If more peaks are found only the first n_peaks are fit, if 
    fewer peaks are found an Exception is raised.

    Manual mode:
    ------------

    When thres, y_shape or x_shape are not defined then initial peak 
    parameters must be provided in y_shapes, x_shapes, y_guesses, x_guesses, 
    and A_guesses.  thres, y_shape, x_shape and n_peaks are ignored in the fit.

    For a list of valid lineshape and corresponding parameter formats see
    gen_lineshape.  When the spectral region is known to contain a single peak
    use fit_single2D.

    """
    
    # Automatic mode
    if None not in [thres,y_shape,x_shape]:
        # automatically find parameters
        t = guess_params_multiple2D(region,y_shape,x_shape,thres,n_peaks)

        # unpack parameter that were found...
        A_guesses,y_guesses,x_guesses = t

        if n_peaks == None:
            n_peaks  = len(A_guesses)
        if len(A_guesses) != n_peaks:
            raise Exception("insufficient number of peaks found")
        y_shapes = [y_shape]*n_peaks
        x_shapes = [x_shape]*n_peaks
    
    # Manual mode
    if None in [y_shapes,x_shapes,y_guesses,x_guesses,A_guesses]:
        raise Exception("incomplete input")
    
    # check parameters
    n_peaks = len(A_guesses)
    if n_peaks!=len(y_guesses) or n_peaks!=len(y_shapes): 
        raise ValueError("length of input Y parameters do not match")
    if n_peaks!=len(x_guesses) or n_peaks!=len(x_shapes):
        raise ValueError("length of input X parameters do not match")
     
    # prepare args downward simplex minimization
    y_len,x_len = region.shape
    ny_params = [len(i) for i in y_guesses]
    nx_params = [len(i) for i in x_guesses]
    args = (region,y_len,x_len,y_shapes,x_shapes,ny_params,nx_params)

    # prepare p0 for simplex
    p0 = []
    for i in A_guesses:
        p0.append(i)
    for y_guess in y_guesses:
        for i in y_guess:
            p0.append(i)
    for x_guess in x_guesses:
        for i in x_guess:
            p0.append(i)

    # performs minimization
    p_best = scipy.optimize.fmin(err_multiple2D,p0,args=args,
                maxiter=maxiter,maxfun=maxfun,disp=disp)

    return unpackp_multiple2D(p_best,ny_params,nx_params)


# p unpacking functions (unpack_typeND)
# 
# These functions are used to unpack the parameter array returned from
# scipy.optimize.fmin into an appropiate tuple

def unpackp_single1D(p):
    """
    Unpack p into (A,x_param)
    """
    return p[0],p[1:]


def unpackp_multiple1D(p,nx_params):
    """
    Unpack p into (As,x_params)
    """ 
    # determind total number of x and y parameters
    total_x = np.array(nx_params).sum()
        
    # extract A and x parameter sets
    n_peaks = len(nx_params)
    As =  p[:n_peaks]
    x_p = p[n_peaks:n_peaks+total_x]

    # order y and x parameters into list of tuples
    x_params = make_tuples(x_p,nx_params)
 
    return As,x_params


def unpackp_single2D(p,ny_param,nx_param):
    """ 
    Unpack p into (A,y_param,x_param)
    """
    return p[0],p[1:ny_param+1],p[-nx_param:]


def unpackp_multiple2D(p,ny_params,nx_params):
    """
    Unpack p into (As,y_params,x_params)
    """
    # determind total number of x and y parameters
    total_x = np.array(nx_params).sum()
    total_y = np.array(ny_params).sum() 
        
    # extract A, y and x parameter sets
    n_peaks = len(ny_params)
    As = p[:n_peaks]
    y_p = p[n_peaks:n_peaks+total_y]
    x_p = p[n_peaks+total_y:n_peaks+total_x+total_y]

    # order y and x parameters into list of tuples
    y_params = make_tuples(y_p,ny_params)
    x_params = make_tuples(x_p,nx_params)
 
    return As,y_params,x_params


def make_tuples(l,t_sizes):
    """
    Create a tuple of tuples of given sizes from a list

    Parameters:
    -----------
    
    * l         List/array to packed tuples.
    * t_sizes   List of tuple sizes.

    """
    
    out = []  # output
    start = 0
    for s in t_sizes:
        out.append(l[start:start+s])
        start = start+s
    return tuple(out)


# err functions (err_typeND)


def err_single1D(p,region,x_len,x_shape):
    """
    Error function for a single 1D peak, to be called by fmin.
    
    Parameters:
    -----------
    
    * p         List of simulation parameters (A,x_param)
    * region    Spectral rgion to fit against.
    * x_len     Length of peak region along X (0) axis.
    * x_shape   Shape of lineshape along X (0) axis.

    Returns:    Sum of squares between peak and simulated peak (float)

    """
    # unpack parameter list
    A,x_param = unpackp_single1D(p)

    # simulate the region
    sim_region = sim_single1D(x_len,x_param,x_shape,A)

    # calculate the error
    return ((region-sim_region)**2).sum()


def err_multiple1D(p,region,x_len,x_shapes,nx_params):
    """
    Error functions for multiple 1D peaks, to be called by fmin.

    Parameters:
    -----------
    
    * p         List of simulation parameters (As,x_params)
    * region    Spectral region to fit against.
    * x_len     Length of peak region along X (0) axis.
    * x_shapes  List of lineshapes along X (0) axis.
    * nx_params List of lengths of peak-parameter which are stored in p.

    Returns:    Sum of squares between peak and simulated peak (float)

    """
    # unpack parameter list
    As,x_params = unpackp_multiple1D(p,nx_params)
   
    # simulate the region
    sim_region = sim_multiple1D(x_len,x_params,x_shapes,As)
    
    # caculate the error
    return ((region-sim_region)**2).sum()


def err_single2D(p,region,y_len,x_len,y_shape,x_shape,ny_param,nx_param):
    """
    Error function for a single 2D peak, to be called by fmin.

    Parameters:
    -----------

    * p         List of simulation parameters (A,y_param,x_param)
    * region    Spectral region to fit against.
    * y_len     Length of peak region along Y (0) axis.
    * x_len     Length of peak region along X (1) axis.
    * y_shape   Shape of lineshape along Y (0) axis.
    * x_shape   Shape of lineshape along X (1) axis.
    * ny_param  Number of y parameters in p.
    * nx_param  Number of x parameters in p.

    Returns:    Sum of squares between peak and simulated peak (float)

    """
    # unpack parameter list
    A,y_param,x_param = unpackp_single2D(p,ny_param,ny_param)

    # simulate the region
    sim_region = sim_single2D(y_len,x_len,y_param,x_param,y_shape,x_shape,A)

    # calculate the error
    return ((region-sim_region)**2).sum()


def err_multiple2D(p,region,y_len,x_len,y_shapes,x_shapes,ny_params,nx_params):
    """
    Error functions for multiple 2D peaks, to be called by fmin.

    Parameters:
    -----------

    * p         List of simulation parameters (A,y_params,x_params)
    * region    Spectral region to fit against.
    * y_len     Length of peak region along Y (0) axis.
    * x_len     Length of peak region along X (1) axis.
    * y_shapes  List of lineshapes along Y (0) axis.
    * x_shapes  List of lineshapes along X (1) axis.
    * ny_params List of lengths of peak-parameter which are stored in p.
    * nx_params List of lengths of peak-parameter which are stored in p.

    Returns:    Sum of squares between peak and simulated peak (float)

    """

    # unpack parameter list
    As,y_params,x_params = unpackp_multiple2D(p,ny_params,nx_params)

    # simulate the region
    sim_region = sim_multiple2D(y_len,x_len,y_params,x_params,y_shapes,
                                x_shapes,As)
    
    # caculate the error
    return ((region-sim_region)**2).sum()


# region simulation functions (sim_typeND)


def sim_single1D(x_len,x_param,x_shape,A=1.):
    """
    Simulate a single peak in 1D region.
    
    Parameter:
    ----------
    
    * x_len     Length of peak region along X (0) axis.
    * x_param   Lineshape parameters of X (0) axis (shape specific).
    * x_shape   Shape of lineshape along X (0) axis.
    * A         Amplitude of peak.

    Returns: 1D array of shape (x_len)

    """    
    return A*sim_lineshape(x_len,x_param,x_shape)


def sim_multiple1D(x_len,x_params,x_shapes,As):
    """
    Simulate multiple peaks in a 1D region.

    Parameters:
    -----------
    
    * x_len     Length of peak region along X (0) axis.
    * x_params  List of lineshape parameters for each peak.
    * x_shapes  List of lineshapes for each peak.
    * As        List of Amplitudes for each peak.

    Returns: 1D array of shape (x_len)

    """
    region = np.zeros(x_len)

    for x_p,x_sh,A in zip(x_params,x_shapes,As):
        peak = sim_single1D(x_len,x_p,x_sh,A)
        region = region + peak

    return region


def sim_single2D(y_len,x_len,y_param,x_param,y_shape,x_shape,A=1.):
    """ 
    Simulate a single peak in a 2D region.

    Parameter:
    ----------
    
    * y_len     Length of peak region along Y (0) axis. 
    * x_len     Length of peak region along X (1) axis.
    * y_param   Lineshape parameters of Y (0) axis (shape specific).
    * x_param   Lineshape parameters of X (1) axis (shape specific).
    * y_shape   Shape of lineshape along Y (0) axis.
    * x_shape   Shape of lineshape along X (1) axis.
    * A         Amplitude of peak.

    Returns:    2D array of shape (y_len,x_len)

    """
    x_lineshape = sim_lineshape(x_len,x_param,x_shape)
    y_lineshape = sim_lineshape(y_len,y_param,y_shape)
    return A*np.outer(y_lineshape,x_lineshape)


def sim_multiple2D(y_len,x_len,y_params,x_params,y_shapes,x_shapes,As):
    """
    Simulate multiple peaks in a 2D region.

    Parameters:
    -----------

    * y_len     Length of peak region along Y (0) axis.
    * x_len     Length of peak region along X (1) axis.
    * y_params  List of lineshape parameters for each peak.
    * x_params  List of lineshape parameters for each peak.
    * y_shapes  List of lineshapes for each peak.
    * x_shapes  List of lineshapes for each peak.
    * As        List of Amplitudes for each peak.

    Returns:    2D array of shape (y_len,x_len)

    """
    region = np.zeros((y_len,x_len))

    for y_p,x_p,y_sh,x_sh,A in zip(y_params,x_params,y_shapes,x_shapes,As):
        peak = sim_single2D(y_len,x_len,y_p,x_p,y_sh,x_sh,A)
        region = region + peak

    return region

# base lineshape functions

# To add a new lineshapes model:
# 1. Make a lineshape1D function
# 2. Add shape case to sim_linshape functions.
# 3. Create parameter packing case in centerFWHM2param_1D


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


def sim_lineshape(M,param,shape,A=1.):
    """
    Generate a 1D linshape of length M

    Parameters:
    -----------
    
    * M         Length of lineshape (0,1,...,M-1).
    * param     Lineshape parameters (shape specific).
    * shape     Shape of lineshape.
    * A         Amplitude of lineshape.

    Returns:    1D array of length M

    Shape       Param           Note
    ----        -----           ----
    'gauss'     (mu,sigma)      Gaussian lineshape
    'g'         (mu,sigma)      Shortcut for 'gauss'
    'lorentz'   (x_0,gamma)     Lorentz lineshape
    'l'         (x_0,gamma)     Shortcut for 'lorentz'


    """    
    if shape == "gauss" or shape == "g":
        mu,sigma = param
        return A*gauss1D(M,mu,sigma)
    elif shape == "lorentz" or shape == "l":
        x_0,gamma = param
        return A*lorentz1D(M,x_0,gamma)
    else:
        raise ValueError("Unknown shape")


# Parameter Guessing functions (guess_params_typeND)


def guess_params_single1D(region,x_shape):
    """
    Guess peak parameters of a single peak in a 1D region

    Parameters:
    -----------

    * region    Spectral region to fit against.
    * xshape    Shape of lineshape along X (0) axis.

    Returns: (A_guess,x_guess)
    --------

    * A_guess   Amplitude guess
    * x_guess   Guess of lineshape parameters of X (0) axis (shape specific).

    """
    # assuming the region containes the whole peak the total integrated area
    # is the amplitude of the peak.
    A_guess = region.sum()

    # the maximum is a good guess at the center of the peak
    x_0 = region.argmax()

    # guess FWHM from 1D trace along through maximum
    x_FWHM = guess_FWHM_1D(region)

    # make the parameter tuples for each dimension
    x_guess = centerFWHM2param_1D(x_0,x_FWHM,x_shape)

    return A_guess,x_guess


def guess_params_multiple1D(region,x_shape,thres,n_peaks=None):
    """
    Guess peak parameters of multiple peaks in a 1D region

    Parameters:
    -----------
    * region    Spectral region to fit against.
    * x_shape   Shape of lineshapes along X (0) axis.
    * thres     Noise threshold.
    * n_peaks   If not None, keep only the largest n_peaks segments.

    Returns: (A_guesses,x_guesses)
    --------
    * A_guesses Amplitude guesses
    * x_guesses Guesses of lineshape parameters of X (0) axis.


    """
    
    # Segment region by finding downward regions.
    centers,segments = find_all_downward1D(region,thres) 

    # cut down segments if needed
    centers,segments = centers[:n_peaks],segments[:n_peaks]

    # initilize _guesses lists
    A_guesses = []
    x_guesses = []

    # loop over segments finding peak parameter for each one
    for x,x_points in zip(centers,segments):
        
        # find limits which encloses segment and extract region
        x0,x1 = np.array(x_points).min(),np.array(x_points).max() 
        line = region[x0:x1+1]
        
        # the sum is a decent guess at the amplitude
        A_guess = line.sum()

        # the maximum is a good guess at the center of the peak
        x_max = line.argmax()

        # guess FWHM from 1D trace along through maximum
        x_FWHM = guess_FWHM_1D(line)

        # make the parameter tuples for each dimension
        x_guess = centerFWHM2param_1D(x,x_FWHM,x_shape)

        # append to output list
        x_guesses.append(x_guess)
        A_guesses.append(A_guess)

    return A_guesses,x_guesses


def guess_params_single2D(region,y_shape,x_shape):
    """
    Guess peak parameters of a single peak in a 2D region.

    Parameters:
    -----------

    * region    Spectral region to fit against.
    * y_shape   Shape of lineshape along Y (0) axis.
    * x_shape   Shape of lineshape along X (1) axis.

    Returns: (A_guess,y_guess,x_guess)
    --------

    * A_guess   Amplitude guess
    * y_guess   Guess of lineshape parameters of Y (0) axis (shape specific).
    * x_guess   Guess of lineshape parameters of X (1) axis (shape specific).

    """
    # assuming the region contains in the entire peak the total integrated area
    # is the amplitude of the peak.
    A_guess = region.sum()

    # the maximum is a good guess at the center of the peak
    y_0,x_0 = np.unravel_index(region.argmax(),region.shape)

    # guess FWHM from 1D trace along through maximum
    y_FWHM = guess_FWHM_1D(region[:,x_0].flatten())
    x_FWHM = guess_FWHM_1D(region[y_0,:].flatten())

    # make the parameter tuples for each dimension
    x_guess = centerFWHM2param_1D(x_0,x_FWHM,x_shape)
    y_guess = centerFWHM2param_1D(y_0,y_FWHM,y_shape)

    return A_guess,y_guess,x_guess


def guess_params_multiple2D(region,y_shape,x_shape,thres,n_peaks=None):
    """
    Guess peak parameters of multiple peaks in a 2D region

    Parameters:
    -----------
    * region    Spectral region to fit against.
    * y_shape   Shape of lineshapes along Y (0) axis.
    * x_shape   Shape of lineshapes along X (1) axis.
    * thres     Noise threshold.
    * n_peaks   If not None, keep only the largest n_peaks segments.

    Returns: (A_guesses,y_guesses,x_guesses)
    --------
    * A_guess   Amplitude guesses
    * y_guess   Guesses of lineshape parameters of Y (0) axis (shape specific).
    * x_guess   Guesses of lineshape parameters of X (1) axis (shape specific).

    """
    
    # Segment region by finding downward regions.
    centers,segments = find_all_downward2D(region,thres) 

    # cut down segments if needed
    centers,segments = centers[:n_peaks],segments[:n_peaks]
    
    # initilize _guesses lists
    A_guesses = []
    y_guesses = []
    x_guesses = []

    # loop over segments
    for ((y,x),(y_points,x_points)) in zip(centers,segments):
        
        # find limits of rectangle which encloses segment and extract region
        y0,y1,x0,x1 = find_box_limits(y_points,x_points)
        peak = region[y0:y1+1,x0:x1+1]
        
        # the sum is a decent guess at the amplitude
        A_guess = peak.sum()

        # the maximum is a good guess at the center of the peak
        y_max,x_max = np.unravel_index(peak.argmax(),peak.shape)

        # guess FWHM from 1D trace along through maximum
        y_FWHM = guess_FWHM_1D(peak[:,x_max].flatten())
        x_FWHM = guess_FWHM_1D(peak[y_max,:].flatten())

        # make the parameter tuples for each dimension
        x_guess = centerFWHM2param_1D(x,x_FWHM,x_shape)
        y_guess = centerFWHM2param_1D(y,y_FWHM,y_shape)

        # append to output list
        x_guesses.append(x_guess)
        y_guesses.append(y_guess)
        A_guesses.append(A_guess)

    return A_guesses,y_guesses,x_guesses


# type independant guessing functions


def guess_FWHM_1D(line):
    """
    Guess the full width half max (FWHM) for a line.

    The method used to find the FWHM in this functions does not perform 
    well if the line does not extend past at least the half-max on each side.

    """
    # Method used:
    # Find the points that bracket the first and last crossing of the 
    # half max, from each use a linear extrapolation to find the location
    # of the half max on either side of the maximum.  The difference between
    # these two values is the full width.
    max = line.max()
    hmax = max/2.

    top_args = np.nonzero(line > hmax)

    l_idx  = top_args[0][0]     # index of left hand side above half-max  
    r_idx = top_args[0][-1]     # index of right hand side above half-max 

    # solve hmax = mx+b => x = y-b/m
    # for two points x_0 and x_1 this becomes y = (hmax-x_0)/(x_1-x_0)
    # to this value add the index of x_0 to get the location of the half-max.

    # the left side of the line
    if l_idx == 0:
        left = l_idx    # this is a bad guess but the best we can do
    else:
        x_0,x_1 = line[l_idx-1],line[l_idx]
        left = l_idx-1+(hmax-x_0)/(x_1-x_0)

    # the right side of the line
    if r_idx == len(line)-1:
        right = r_idx   # this is a bad guess but the best we can do
    else:
        x_0,x_1 = line[r_idx],line[r_idx+1]
        right = r_idx+(hmax-x_0)/(x_1-x_0)

    return right-left


def centerFWHM2param_1D(center,FWHM,shape):
    """
    Convert a center and FWHM into lineshape parameters

    Parameters:
    -----------

    * center    Location of center of peak.
    * FWHM      Full width half max of peak.
    * shape     Shape of lineshape

    Returns:    parameter tuple according to shape

    """

    if shape == 'gauss' or shape == 'g':
        # mu is the center, FWHM = 2*sqrt( 2*ln(2) )*sigma = 2.3548*sigma
        return (center,FWHM/2.3548)
    elif shape == 'lorentz' or shape == 'l':
        # center is x_0, FWHM is gamma
        return (center,FWHM)
    else:
        raise ValueError("Unknown shape")
 

####################################
# Spectrum and region segmentation #
####################################

def find_box_limits(y_points,x_points):
    """ 
    Find the box limits which will outline the provided points

    Parameters:
    -----------

    * y_points  Array of Y (0) axis indices.
    * x_points  Array of X (1) axis indices.

    Returns: y_min,y_max,x_min,x_max

    """
    return y_points.min(),y_points.max(),x_points.min(),x_points.max()


# Downward segmenting method:
# The downward segmenting method uses the flood fill algorithm to find
# all points connected to an initial node which are above a given threshold
# and to which a path exists in which each step of the path moves lower in 
# intensity.  This can be though of as all points accessible by a water drop
# following downward slopes from the initial node.


def find_all_downward1D(data,thres):
    """
    Find all downward-connected segments in 1D data. 

    Parameter:
    ----------

    * data  1D array of data
    * thres Threshold, below this nodes are considered noise.

    Returns: (centers,segments)
    ---------------------------

    * centers   List of local maximum of segments
    * segments  List of arrays of points in each segment.
    
    """

    # create the masked data array
    mdata = np.ma.masked_less(data,thres)
   
    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all() != True:

        # find the maximum (center of segment)
        x = mdata.argmax()

        # find all downward points and mark them
        xs = find_downward1D(mdata,x,thres)
        mdata[(xs)] = ma.masked

        # add to centers and segments lists
        centers.append(x)
        segments.append(xs)

    return centers,segments


def find_downward1D(data,x,thres):
    """
    Find points downward-connected to node (x) in data.

    Parameters:
    -----------

    * data  1D array of data
    * x     X (0) axis index to starting node.
    * thres Threshold, below this nodes are considered noise.

    Return: x_array
    -------
    
    * x_array   Array of x-axis indices of connected nodes.

    """

    if data.ndim != 1:
        raise ValueError("data must be 1 dimensional")

    if data[x] <= thres:
        return []
    
    points = [x]

    # search left
    curr_x = x  # current x value
    while curr_x != 0 and thres<data[curr_x-1]<data[curr_x]:
        points.append(curr_x-1)
        curr_x = curr_x - 1
    # search right
    curr_x = x  # current x value
    while curr_x != len(data) and thres<data[curr_x+1]<data[curr_x]:
        points.append(curr_x+1)
        curr_x = curr_x + 1

    return points


def find_all_downward2D(data,thres):
    """
    Find all downward-connected segments in 2D data.

    Parameter:
    ----------

    * data  2D array of data
    * thres Threshold, below this nodes are considered noise.

    Returns: (centers,segments)
    --------

    * centers   List of (y,x) tuples containing local maximum of segments
    * segments  List of (y_array,x_array) for each segments
    
    """
    # create the masked data array
    mdata = np.ma.masked_less(data,thres)
   
    centers = []
    segments = []

    # loop and fill map until all points are masked
    while mdata.mask.all() != True:

        # find the maximum (center of segment)
        y,x = np.unravel_index(mdata.argmax(),mdata.shape)

        # find all downward points and mark them
        ys,xs = find_downward2D(mdata,y,x,thres)
        mdata[(ys,xs)] = ma.masked

        # add to centers and segments lists
        centers.append( (y,x) )
        segments.append( (ys,xs) )

    return centers,segments

def find_downward2D(data,y,x,thres):
    """
    Find point downward-connected to node (y,x) in data.

    Parameters:
    -----------

    * data  2D array of data
    * y     Y (0) axis index to starting node.
    * x     X (1) axis index to starting node.
    * thres Threshold, below this nodes are considered noise.

    Return: (y_array,x_array)
    -------

    * y_array   Array of y-axis indices of connected nodes.
    * x_array   Array of x-axis indices of connected nodes. 

    """
    if data.ndim != 2:
        raise ValueError("data must be 2 dimensional")

    if data[y,x] <= thres:  # check if x,y is above threshold
        return (np.array([]),np.array([]))
    
    # initilize
    top = data.shape[0]     # index of top of array
    right = data.shape[1]   # index of right side of array
    Q = [(y,x)]             # queue
    points = [(y,x)]        # list of connected nodes

    # queue based flood fill algorithm
    while Q:    # loop until Q is empty
        y,x = Q.pop(0)  # remove first element of queue (already in points)
        v = data[y,x]   # value at current node

        # check all four directions if above thres, less than current point, 
        # and above thres, if so add to queue and add to points list

        # north
        if y+1 != top and thres<data[y+1,x]<v and (y+1,x) not in points:
            Q.append((y+1,x))
            points.append((y+1,x))
        # south
        if y != 0 and thres<data[y-1,x]<v and (y-1,x) not in points:
            Q.append((y-1,x))
            points.append((y-1,x))
        # east
        if x+1 != right and thres<data[y,x+1]<v and (y,x+1) not in points:
            Q.append((y,x+1))
            points.append((y,x+1))
        # west
        if x != 0 and thres<data[y,x-1]<v and (y,x-1) not in points:
            Q.append((y,x-1))
            points.append((y,x-1))

    return (np.array(points)[:,0],np.array(points)[:,1])


# Connected segmenting method:
# The connected segmenting method uses the flood fill algorithm to find
# all points connected to an initial node which are above a given threshold.


def find_connected2D(data,y,x,thres):
    """
    Find all points connected to node (x,y) above threshold.

    Parameters:
    -----------

    * data  2D array of data
    * y     Y (0) axis index to starting node.
    * x     X (1) axis index to starting node.
    * thres Threshold, below this nodes are considered noise.

    Return: (y_array,x_array)
    -------

    * y_array
    * x_array

    """
    if data.ndim != 2:
        raise ValueError("data must be 2 dimensional")

    if data[y,x] <= thres:  # check if x,y is above threshold
        return []
    
    # initilize
    top = data.shape[0]     # index of top of array
    right = data.shape[1]   # index of right side of array
    Q = [(y,x)]             # queue
    points = [(y,x)]        # list of connected nodes

    # queue based flood fill algorithm
    while Q:    # loop until Q is empty
        y,x = Q.pop(0)  # remove first element of queue (already in points)
        # check all four directions if above thres and not in points, add to 
        # queue and add to points
        # north
        if y+1 != top and data[y+1,x] > thres and (y+1,x) not in points:
            Q.append((y+1,x))
            points.append((y+1,x))
        # south
        if y != 0 and data[y-1,x] > thres and (y-1,x) not in points:
            Q.append((y-1,x))
            points.append((y-1,x))
        # east
        if x+1 != right and data[y,x+1] > thres and (y,x+1) not in points:
            Q.append((y,x+1))
            points.append((y,x+1))
        # west
        if x != 0 and data[y,x-1] > thres and (y,x-1) not in points:
            Q.append((y,x-1))
            points.append((y,x-1))

    return (np.array(points)[:,0],np.array(points)[:,1])


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

    # loop and fill map until all points are masked
    while mdata.mask.all() != True:
        y,x = np.unravel_index(mdata.argmax(),mdata.shape)
        mark_map2D(mdata,map,y,x,mark)
        mark = mark+1
        print mark,y,x
    return map


def mark_map2D(mdata,map,y,x,mark):
    """
    Mark connected region on segment map starting at node (y,x)

    This functions should be called from map_segments2D

    Parameters:
    -----------
    * mdata Masked 2D data array
    * map   2D integer array mapping out segments
    * y     Starting node location in 0 axis
    * x     Starting node location in 1 axis 
    * mark  Integer to mark map with.

    Returns nothing but modifies mdata mask and map.

    """
   
    # Algorithm is similar to the 2nd presented in:
    # http://en.wikipedia.org/wiki/Flood_fill

    # index limits of array
    right = mdata.shape[1]
    top = mdata.shape[0]
    
    # If node is masked return (shouldn't happen).
    if mdata.mask[y,x] == True:
        return 

    # create the queue with the node
    Q = [(y,x)]

    # loop until Q is empty
    while Q:
        # remove the first element of the queue
        y,x = Q.pop(0)

        # if working node is not masked, mark it and mask it 
        if mdata.mask[y,x] == False:
            print "Marking:",y,x
            map[y,x] = mark
            mdata[y,x] = ma.masked
        
        # Check all four directions to see if they are not masked, if
        # so add to queue, mark, and mask them
        # north
        if y+1 != top and mdata.mask[y+1,x] == False:
            #print "Marking and adding",y+1,x
            Q.append((y+1,x))
            map[y+1,x] = mark
            mdata[y+1,x] = ma.masked
        # south
        if y != 0 and mdata.mask[y-1,x] == False:
            #print "Marking and adding",y-1,x
            Q.append((y-1,x))
            map[y-1,x] = mark
            mdata[y-1,x] = ma.masked
        # east
        if x+1 != right and mdata.mask[y,x+1] == False:
            #print "Marking and adding",y,x+1
            Q.append((y,x+1))
            map[y,x+1] = mark
            mdata[y,x+1] = ma.masked
        # west
        if x != 0 and mdata.mask[y,x-1] == False:
            #print "Marking and adding",y,x-1
            Q.append((y,x-1))
            map[y,x-1] = mark
            mdata[y,x-1] = ma.masked
    return

################
# Peak picking #
################

class ndwindow_inside(object):
    """
    An 
    Given the shape of an rray and a window size, an 'ndwindow_inside' 
    instance iterators over tuples of slices which slice an the array into 
    uniform size wsize sub-arrays.  At each iteration, the index of the top 
    left of the sub-array is incremented by one along the last dimension util
    the windows would extend past the array border.  All sub-arrays are
    equal sized (wsize).

    Parameters
    ----------

    * size  Size of array to generate tuples of slices from.
    * wsize Size of the area to select from array (widow size).

    Example
    -------
    
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


class ndwindow(object):
    """ 
    An N-dimentional iterator object to slice arrays into windows.

    Given the shape of an array and a window size, an 'ndwindow' instance
    iterators over tuples of slices which slice an the array into wsize 
    sub-arrays.  At each iteration, the index of the center of the sub-array 
    is incremented by one along the last dimension.  Array border are ignored
    so the resulting sub-array can be smaller than wsize.  If wsize contains
    even values the window in off center containing more points of lower
    index in the even value dimensions.

    Parameters
    -----------
    * size  Size of array to generate tuples of slices from.
    * wsize Size of the area to select from array (sub-array maximum size).
    
    Example
    -------

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
    etc.


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
        return (slices,index)

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
    wsize = tuple([2*i+1 for i in msep])    #window size is 2*seperation+1

    middle = np.floor((np.array(msep)-1)/2.)    # index of middle of window
    ms = [slice(x,x+1) for x in middle]         # middle slice list

    # loop over the windows
    for idx,s in ndwindow_index(data.shape,wsize):
        window = data[s]
        max = window.max()
        #print idx
        if max > thres and max == data[idx]:
            peaks.append(idx)
    return np.array(peaks)


def peakpick_thres_fast(data,thres,msep):
    """
    Fast version of peakpick_thres functions
    """
    # window size is 2*seperation+1
    size = tuple([2*i+1 for i in msep])

    # find local maxima mask
    mx = scipy.ndimage.maximum_filter(data,size=size,mode='constant')
    max_mask = mx == data
    
    # find threshold mask
    thres_mask = np.ma.greater(data,thres)

    # peaks are bitwise anding of these masks
    return np.transpose(np.nonzero(thres_mask*max_mask))
