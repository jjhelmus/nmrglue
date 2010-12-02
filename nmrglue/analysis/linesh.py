"""
Functions for fitting and simulating arbitrary dimensional lineshapes commonly
found in NMR experiments

"""

import numpy as np
from scipy.optimize import leastsq
from leastsqbound import leastsqbound

pi = np.pi

# lineshape classes translator
from analysisbase import ls_str2class
from analysisbase import pick2linesh,linesh2pick

# User facing functions
def sim_NDregion(shape,lineshapes,params,amps):
    """
    Simulate a arbitrary dimensional region with one or more peaks.

    Parameters:

    * shape         tuple of region shape
    * lineshapes    List of lineshapes by label (str) or a lineshape class.
                    See fit_NDregion for additional documentation.
    * params        P-length list (P is the number of peaks in region) of 
                    N-length lists of tuples where each each tuple is 
                    lineshape parameters for a given peak and dimension.
    * amps          P-length of peak amplitudes.

    Returns: array containing simulated region

    """
    # parse the user-friendly input into a format digestable by s_NDregion
    
    # parse the shape
    ndim = len(shape)
    
    # parse the lineshape parameters
    if len(lineshapes) != ndim:
        raise ValueError("Incorrect number of lineshapes provided")

    ls_classes = []
    for l in lineshapes:
        if type(l) is str:
            ls_classes.append(ls_str2class(l))
        else:
            ls_classes.append(l)

    # determind the number of parameters in each dimension.
    dim_nparam = [c.nparam(l) for l,c in zip(shape,ls_classes)]

    # parse the params parameter
    n_peaks = len(params)
    p = []
    for i,param in enumerate(params):
        if len(param) != ndim:
            err = "Incorrect number of parameters for peak %i"
            raise ValueError(err%(i))
        for j,dim_param in enumerate(param):
            if len(dim_param) != dim_nparam[j]:
                err = "Incorrect number of parameters in peak %i dimension %i"
                raise ValueError(err%(i,j))

            for g in dim_param:
                p.append(g)

    # parse the amps parameter
    if len(amps) != n_peaks:
        raise ValueError("Incorrect number of amplitudes provided")
    p = list(amps) + p # amplitudes appended to front of p

    # DEBUGGING
    #print "p",p
    #print "shape",shape
    #print "ls_classes",ls_classes
    #print "n_peaks",n_peaks

    return s_NDregion(p,shape,ls_classes,n_peaks)


def fit_NDregion(region,lineshapes,guesses,amp_guesses,guesses_bounds=None,
                 amp_bounds=None,error_flag=False,**kw):
    """
    Fit a N-dimensional region.

    Parameters:
    
    * region         N-dimensional region to fit.
    * lineshapes     List of lineshapes by label (str) or a lineshape class.
    * guesses        P-length list (P is the number of peaks in region) of 
                     N-length lists of tuples where each each tuple is the 
                     optimiztion starting parameters for a given peak and 
                     dimension lineshape.
    * amp_guesses    P-length list of amplitudes.
    * guesses_bounds List of bounds for parameter of same shape as guesses.  If
                     none of the parameters in a given dimension have limits 
                     None can be used, otherwise each dimension should have a 
                     list/tuple of (min,max) or None for each parameter.  
                     min or max may be None when there is no bound in a given 
                     direction.
    * amp_bounds     P-length list of bounds for the amplitude with format 
                     similar to guesses_bound.
    * error_flag     Set to True to estimate errors for each lineshape 
                     parameter and amplitude.

    * kw             Additional keywords passed to the scipy.optimize.leastsq 
                     function.

    Returns: param_best,amp_best,ier

    * params_best   Optimal values for lineshape parameters with same format
                    as guesses input parameter.
    * amp_best      List of optimal peak amplitudes.
    * ier           Interger flag from scipy.optimize.leastsq indicating if
                    the solution was found.  1,2,3,4 indicates that a solution
                    was found.  Otherwise the solution was not found.


    Note on the lineshape parameter:

    Elements of the lineshape parameter list can be string indicating the
    lineshape of given dimension or an instance of a lineshape class 
    which provide a sim method which takes two arguments, the first being the 
    length of the lineshape the second being a list of lineshape parameters, 
    and returns a simulated lineshape as well as a nparam method which when 
    given the length of lineshape returns the number of parameters needed to
    describe the lineshape. Currently the following strings are allowed:

    * 'g' or 'gauss'    Gaussian (normal) lineshape.
    * 'l' or 'lorentz'  Lorentzian lineshape.
    * 's' or 'scale'    Scaled lineshape.
    * 'p' or 'peak'     Gaussian lineshape which takes FWHM as a parameter.

    The following are all valid lineshapes parameters for a 2D Gaussian peak:

    ['g','g']
    ['gauss','gauss']
    [ng.analysisbase.gauss1D(),ng.analysisbase.gauss1D()]
    
    An simple example of a lineshape class which simulates the function y=c:

    class constant(): 
        def sim(self,M,p):
            c = p
            return c*np.ones(M)
        def nparam(self,M):
            return 1

    
    """
    # this function parses the user-friendly input into a format digestable
    # by f_NDregion, performs the fitting, then format the fitting results
    # into a user friendly format

    # parse the region parameter
    ndim = region.ndim
    shape = region.shape

    # parse the lineshape parameter
    if len(lineshapes) != ndim:
        raise ValueError("Incorrect number of lineshapes provided")
    
    ls_classes = []
    for l in lineshapes:
        if type(l) is str:
            ls_classes.append(ls_str2class(l))
        else:
            ls_classes.append(l)
    
    # determind the number of parameter in each dimension
    dim_nparam = [c.nparam(l) for l,c in zip(shape,ls_classes)] 

    # parse the guesses parameter
    n_peaks = len(guesses)
    p0 = []
    for i,guess in enumerate(guesses):  # peak loop
        if len(guess) != ndim:
            err = "Incorrect number of guesses for peak %i"
            raise ValueError(err%(i))
        
        for j,dim_guess in enumerate(guess):    # dimension loop
            if len(dim_guess) != dim_nparam[j]:
                err = "Incorrect number of parameters in peak %i dimension %i"
                raise ValueError(err%(i,j))
            
            for g in dim_guess: # parameter loop
                p0.append(g)

    
    # parse the bounds parameter
    if guesses_bounds == None:   # No bounds 
        peak_bounds = [[(None,None)]*i for i in dim_nparam]
        guesses_bounds = [peak_bounds]*n_peaks

    if len(guesses_bounds) != n_peaks:
        raise ("Incorrect number of parameter bounds provided")

    # build the parameter bound list to be passed to f_NDregion
    p_bounds = []
    for i,peak_bounds in enumerate(guesses_bounds): # peak loop
        
        if peak_bounds == None:
            peak_bounds = [[(None,None)]*i for i in dim_nparam]
        
        if len(peak_bounds) != ndim:
            err = "Incorrect number of bounds for peak %i"
            raise ValueError(err%(i))
        
        for j,dim_bounds in enumerate(peak_bounds):    # dimension loop
            
            if dim_bounds == None:
                dim_bounds = [(None,None)]*dim_nparam[j]
            
            if len(dim_bounds) != dim_nparam[j]:
                err = "Incorrect number of bounds for peak %i dimension %i"
                ValueError(err%(i,j))

            for k,b in enumerate(dim_bounds):    # parameter loop
                if b == None:
                    b = (None,None)

                if len(b) != 2:
                    err  = "No min/max for peak %i dim %i parameter %i"
                    raise ValueError(err%(i,j,k))
                
                p_bounds.append(b)
    
    # parse the amp_guesses parameter
    if len(amp_guesses) != n_peaks:
        raise ValueError("Incorrect number of amplitude guesses provided")
    p0 = list(amp_guesses) + p0 # amplitudes appended to front of p0
   
    # parse the amp_bounds parameter
    if amp_bounds  == None:
        amp_bounds = [(None,None)]*n_peaks

    if len(amp_bounds) != n_peaks:
        raise ValueError("Incorrect number of amplitude bounds")

    to_add = []
    for k,b in enumerate(amp_bounds):
        if b == None:
            b = (None,None)

        if len(b) != 2:
            err = "No min/max for amplitude bound %i"
            raise ValueError(err%(k)) 
        to_add.append(b)
    p_bounds = to_add + p_bounds    # amplitude bound at front of p_bounds

    
    # DEBUGGING
    #print "--------------------------------"
    #print region
    #print ls_classes
    #print p0
    #print p_bounds
    #print n_peaks

    # include full_output=True when errors requested 
    if error_flag:
        kw["full_output"] = True

    
    # perform fitting
    r = f_NDregion(region,ls_classes,p0,p_bounds,n_peaks,**kw)

    # DEBUGGING
    #print r

    # unpack results depending of if full output requested
    if "full_output" in kw and kw["full_output"]:
        p_best,cov_xi,infodic,mesg,ier = r
    else:
        p_best,ier = r
    
    # unpack and repack p_best
    # pull off the ampltides
    amp_best = p_best[:n_peaks]
    
    # split the remaining parameters into n_peaks equal sized lists
    p_list = split_list(list(p_best[n_peaks:]),n_peaks)

    # for each peak repack the flat parameter lists to reference by dimension
    param_best = [make_slist(l,dim_nparam) for l in p_list]

    # return as is if no errors requested
    if error_flag==False:
        return param_best,amp_best,ier

    # calculate errors
    p_err = calc_errors(region,ls_classes,p_best,cov_xi,n_peaks)

    # unpack and repack the error p_err
    # pull off the amplitude errors
    amp_err = p_err[:n_peaks]
    
    # split the remaining errors into n_peaks equal sized lists
    pe_list = split_list(list(p_err[n_peaks:]),n_peaks)
    
    # for each peak repack the flat errors list to reference by dimension
    param_err = [make_slist(l,dim_nparam) for l in pe_list]
    
    return param_best,amp_best,param_err,amp_err,ier


def make_slist(l,t_sizes):
    """
    Create a list of tuples of given sizes from a list

    Parameters:

    * l         List/array to pack into shaped list.
    * t_sizes   List of tuple sizes.

    """
    out = []  # output
    start = 0
    for s in t_sizes:
        out.append(l[start:start+s])
        start = start+s
    return out

def split_list(l,N):
    """ Split list l into N sublists of equal size """
    step = int(len(l)/N)
    div_points = range(0,len(l)+1,step)
    return [l[div_points[i]:div_points[i+1]] for i in xrange(N)]


def calc_errors(region,ls_classes,p,cov,n_peaks):
    """
    Calcuate the parameter errors from the Standard Errors of the Estimate.

    Parameters:
    
    * region        N-dimensional region to fit
    * ls_classes    List of lineshape classes
    * p             Parameters (array)
    * cov           Covariance matrix from leastsq fitting
    * n_peaks       Number of peaks in region

    Returns: array of standard errors of parameters in p

    """

    # calculate the residuals
    resid = err_NDregion(p,region,region.shape,ls_classes,n_peaks)
    
    SS_err = np.power(resid,2).sum()    # Sum of squared residuals
    n = region.size # size of sample XXX not sure if this always makes sense
    k = p.size-1    # free parameters
    st_err = np.sqrt(SS_err/(n-k-1))    # standard error of estimate
    
    return st_err*np.sqrt(np.diag(cov))

# internal functions

def s_NDregion(p,shape,ls_classes,n_peaks):
    """
    Simulate a arbitrary dimensional region with one or more peaks.

    Parameters:

    * p             List (and must be a list) of parameters
    * shape         tuple of region shape
    * ls_classes    List of lineshape classes
    * n_peaks       Number of peaks in region

    p is modified by this functions, pass a copy if p should be retained.

    """    
    # split the parameter list into a list of amplitudes and peak param lists
    As = [p.pop(0) for i in xrange(n_peaks)]   
    ps = split_list(p,n_peaks)
    
    # simulate the first region
    A,curr_p = As.pop(0),ps.pop(0)
    r = s_single_NDregion([A]+curr_p,shape,ls_classes)
    
    # simulate any additional regions
    for A,curr_p in zip(As,ps):
        r = r+s_single_NDregion([A]+curr_p,shape,ls_classes)
    
    return r

def s_single_NDregion(p,shape,ls_classes):
    """
    Simulate a arbitrary dimensional region with a single peak. Called 
    repeatly by s_NDregion to build up a full ND region.

    Parameters:
    * p             List (and must be a list) of parameters
    * shape         tuple of region shape
    * ls_classes    List of lineshape classes

    """
    
    A = p.pop(0)    # amplitude is ALWAYS the first parameter
    r = np.array(A,dtype='float')

    for length,ls_class in zip(shape,ls_classes):
        #print "Making lineshape of",ls_class.name,"with length:",length
        s_p = [p.pop(0) for i in xrange(ls_class.nparam(length))]
        ls = ls_class.sim(length,s_p)
        #print "Lineshape is:",ls
        r = np.kron(r,ls)   # vector direct product flattened

    return r.reshape(shape)

def err_NDregion(p,region,shape,ls_classes,n_peaks):
    """
    Error functions for a NDregion, called by f_NDregion function
    """
    sim_region = s_NDregion(list(p),shape,ls_classes,n_peaks)
    return (region-sim_region).flatten()

def f_NDregion(region,ls_classes,p0,p_bounds,n_peaks,**kw):
    """
    Fit an arbitrary dimensional regions  containing one or more peaks 
    using a contrained Levenberg-Marquard optmization algorithm.

    Parameters:

    * region        N-dimensional region to fit
    * ls_classes    List of lineshape classes
    * p0            Initial parameter guesses
    * p_bounds      (min,max) pairs for each element of p0
    * n_peaks       Number of peaks
    
    Additional keyword are passed directly to leastsqbound and in turn
    passed to scipy.optimize.leastsq after variable transformation.

    """
    args = (region,region.shape,ls_classes,n_peaks)
    p_best = leastsqbound(err_NDregion,p0,bounds=p_bounds,args=args,**kw)
    return p_best

