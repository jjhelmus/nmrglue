"""
Linear Prediction (LP) functions for extrapolating and modeling 1D NMR signals

"""

# Developer notes:

# This module contains functions for performing linear prediction on NMR 
# data, the algorithms used were selected for simplicity to show how linear
# prediction works not for computational speed nor stability.  Locations 
# where significant improvements can be made to improve speed or stability
# are indicated with SPEED and STABILITY with discussion following.

# The notation for the Linear Prediction equation, coefficients, roots, etc
# closely match those mentioned in "NMR Data Processing" by Hoch snd Stern.
# This book was references for many of the algorithms in this module.

# Reduced order LP-SVD and LP-TLS methods are not implemented but should
# be easy to add if desired.  See find_lproots_hsvd for an example.

# Possible additions
# mirror mode
# Cadzow Procedure
# 2d-LP

import numpy as np
import scipy
import scipy.linalg

###########################################
# Top level parametric modeling functions #
###########################################

def lp_model(trace,slice=slice(None),order=8,mode="f",method="svd",full=False):

    """
    Use Linear Prediction to model a NMR signal

    Parameter:

    trace   1D trace of data, the FID.
    slice   slice object selecting region of slice to use in LP equation.
    order   Prediction order (Number of LP coefficients)
    mode    Mode to generate LP filter (f-forward,b-backward,
    method  Method to use to calculate the LP filter, choose from
            'svd','qr','choleskey','hsvd'.
    full    Set to True to return amplitudes and phases calculated
            by performing a least squares fitting. False returns only
             the damping (relaxation) factors and signal frequencies


    Returns:    when full is False:  (damp,freq)
                when full if True:   (damp,freq,amp,phase)
        
        Where:

        damp    List of damping factors
        freq    List of frequencies
        amp     List of amplitudes
        phase   List of phases

    Notes: 

    When backward LP is used the signal roots are reflected before calculating
    model parameters.

    """
    # check for bad arguments
    if mode not in ["f","b"]:
        raise ValueError("mode must be 'f' or 'b'")
    if method not in ['svd','qr','cholesky','tls','hsvd']:
        raise ValueError("Invalid method")

    x = trace[slice]    # extract region to use for finding LP coefficients

    # calculate LP coefficient and factor to find poles
    if method in ['svd','qr','cholseky','tls']: 
        D,d = make_Dd(x,order,mode) # form the LP equation elements
        a = find_lpc(D,d,method)    # find LP coefficients
        poles = find_roots(a,mode)  # find roots

    elif method == "hsvd":
        poles = find_lproots_hsvd(x,M=order,K=order,mode=mode,zmethod='sm')
    else:
        raise ValueError("Invalid method")

    # reverse poles if we have backward poles 
    if mode == "b":
        poles = [1./pole for pole in poles]

    # determind the damping factor and frequencies from the roots
    damp = [root2damp(pole) for pole in poles]
    freq = [root2freq(pole) for pole in poles]
    
    if full == False:
        return damp,freq

    # perform Least Squares fitting to determind amplitudes and phases.

    # We need to find a least squares solutions to:
    # z_0*b_0^0+z_1*b_1^0+.... = x_0
    # z_0*b_0^1+z_1*b_1^1+.... = x_1
    # ...
    # Where b is the LP roots (poles), x is the signal, and a_0 are unknown
    # 
    # To solve this recast into Bz = x and solve for a using scipy.lstsq

    # SPEED
    # B is a Vandermonde matrix, this characteristic can be used to more
    # efficiently solve this least squares problem.

    # build the B matrix (a Vandermonde matrix) and solve for the coefficients
    poles = np.array(poles)
    B = np.row_stack([poles**(i) for i in range(len(x))])
    z,resid,rank,s = scipy.linalg.lstsq(B,np.array(x))

    # Now the z_n = amp_n*exp(phase_n*i), use this to determind the amplitudes
    # and phases
    amp   = [cof2amp(cof)   for cof in z]
    phase = [cof2phase(cof) for cof in z] 

    return damp,freq,amp,phase


###############################################################
# functions to determine signal parameters from LP parameters # 
###############################################################

def root2damp(pole):
    """ Calculate the damping factor from a LP root """
    # damping factor is the radius in the complex plane -1/pi*ln(|pole|)
    return -1./(np.pi)*np.log(np.abs(pole))

def root2freq(pole):
    """ Calculate the frequency from a LP root """
    # frequency is the angle from the x-axis to the point in the complex plane
    # arg(pole) / (2*pi) = atan2(imag,real) / (2*pi)
    return np.arctan2(pole.imag,pole.real)/(2.*np.pi)

def cof2amp(z):
    """ Calculate a signal amplitude from a model coefficient """
    # z = amp*exp(phase*i) so amp is abs(z)
    return abs(z)

def cof2phase(z):
    """ Calculate a signal phase from a model coefficient """
    # z = amp*exp(phase(i) so phase is the arg(z)
    return np.arctan2(z.imag,z.real)

#####################################
# Top Level extrapolation functions #
#####################################

def lp(trace,pred=1,slice=slice(None),order=8,mode="f",append="after",
    bad_roots="auto",fix_mode="on",method="svd"):
    """
    Linear Prediction to extrapolate data beyond original data

    Parameter:

    trace       1D trace of data, the FID.
    pred        Number of points to predict.
    slice       slice object selecting region of slice to use in LP equation.
    order       Prediction order (Number of LP coefficients)
    mode        Mode to generate LP filter (f-forward,b-backward,
                fb-"forward-backward",bf-"backward-forward")
    extend      Extend trace before or after trace, ("before" or "after").
    bad_roots   Type of roots which are bad and should be stabilized, 
                either those with "incr" or "decr" signals, set to None for 
                no root stabilizing. Default of "auto" set root fixing 
                based on LP mode. 
                (mode=="f" or "fb" -> "incr", mode=="b" or "bf" -> "decr")
    fix_mode    Method used to stabilize bad roots, "on" to move the roots 
                onto the unit circle, "reflect" to reflect bad roots across 
                the unit circle
    method      Method to use to calculate the LP filter, choose from
                'svd','qr','choleskey'.


    Notes: 

    In forward-backward or backward-forward mode root stabilizing is done 
    on both sets of signal roots as calculated in the first mode direction.  
    After averaging the coefficient the roots are again stabilized.

    When the append parameter does not match the LP mode, for example
    if a backward linear prediction (mode='b') is used to predict points
    after the trace (append='after'), any root fixing is done before reversing 
    the filter.

    """
    # check for bad arguments
    if mode not in ["f","b","fb","bf"]:
        raise ValueError("mode must be 'f','b', 'fb', or 'bf'")
    if append not in ["before","after"]:
        raise ValueError("append must be 'before' or 'after'")
    if bad_roots not in [None,"incr","decr","auto"]:
        raise ValueError("bad_roots must be None, 'auto', 'incr' or 'decr'")
    if fix_mode not in ["on","reflect"]:
        raise ValueError("fix_mode must be 'on' or 'reflect'")
    if method not in ['svd','qr','cholesky','tls']:
        raise ValueError("Invalid method")

    # bad_roots auto mode
    if bad_roots == "auto":
        if mode == "f" or mode =="fb":
            bad_roots = "incr"
        else:
            bad_roots = "decr" 

    x = trace[slice]    # extract region to use for finding LP coefficients

    if mode == "fb":
        a = find_lpc_fb(x,order,bad_roots,fix_mode,method)
        mode = "f"
    elif mode == "bf":
        a = find_lpc_bf(x,order,bad_roots,fix_mode,method)
        mode = "b"
    else:
        D,d = make_Dd(x,order,mode) # form the LP equation matrix and vector
        a = find_lpc(D,d,method)    # determind the LP prediction filter

    # stablize roots if needed
    if bad_roots != None:           # stablize roots if needed
        poles = find_roots(a,mode)  # find roots (poles)
        poles = fix_roots(poles,bad_roots,fix_mode) # fix roots
        # reverse filter when calculated filter is in wrong direction
        if (mode=="b" and append=="after")or(mode=="f" and append=="before"):
            poles = [1./pole for pole in poles]
            mode = {'f':'b','b':'f'}[mode]
        
        a = find_coeff(poles,mode)  # find LP filter from roots
    else:
        # reverse filter when calculated filter is in wrong direction
        if (mode=="b" and append=="after") or (mode=="f" and append=="before"):
            a = reverse_filter(a,mode)

    # extrapolate the trace using the prediction filter
    ntrace = extrapolate(trace,a,pred,append)

    return ntrace

def lp_svd(trace,pred=1,slice=slice(None),order=8,mode="f",append="after",
    bad_roots="auto",fix_mode="on"):
    """
    Linear Prediction using SVD decomposition to extrapolate data.
    
    See lp function for information on parameters.

    """
    return lp(trace,pred,slice,order,mode,append,bad_roots,fix_mode,"svd")

def lp_qr(trace,pred=1,slice=slice(None),order=8,mode="f",append="after",
    bad_roots="auto",fix_mode="on"):
    """
    Linear Prediction using a QR decomposition. to extrapolate data.
    
    See lp function for information on parameters.

    """
    return lp(trace,pred,slice,order,mode,append,bad_roots,fix_mode,"qr")

def lp_cholesky(trace,pred=1,slice=slice(None),order=8,mode="f",append="after",
    bad_roots="auto",fix_mode="on"):
    """
    Linear Prediction using a Cholesky decomposition to extrapolate data.
    
    See lp function for information on parameters.

    """
    return lp(trace,pred,slice,order,mode,append,bad_roots,fix_mode,"cholesky")

def lp_tls(trace,pred=1,slice=slice(None),order=8,mode="f",append="after",
    bad_roots="auto",fix_mode="on"):
    """
    Linear Prediction using the total least squares method to extrapolate data.
    
    See lp function for information on parameters.
    """
    return lp(trace,pred,slice,order,mode,append,bad_roots,fix_mode,"tls")


##############################
# data preperation functions #
##############################

def make_D(x,order,mode):
    """ 
    make the LP equation D matrix (Da = d')
    """
    L = len(x)-order
    if mode == "f":
        return scipy.linalg.hankel(x[:L],x[L-1:-1])
    elif mode == "b":
        return scipy.linalg.hankel(x[1:L+1],x[L:])
    else:
        raise ValueError("mode must be 'f' or 'b'")

def make_d(x,order,mode):
    """
    make the LP equation d' vector (Da = d')
    """
    if mode == "f":
        return x[order:].reshape(len(x)-order,1)
    elif mode == "b":
        L = len(x)-order
        return x[:L].reshape(L,1)
    else:
        raise ValueError("mode must be 'f' or 'b'")

def make_Dd(x,order,mode):
    """
    make the LP equation D matrix and d' vector (Da=d')
    """
    return make_D(x,order,mode),make_d(x,order,mode)

###########################################################
# LP prediction filter calculation functions (find_lpc_*) #
###########################################################

# the coefficients returned from these functions depend on the mode of the
# prediction.  Forward LP returns coefficients ordered m,m-1,...1
# Backward LP returns 1,2,...,m where m is the order of the prediction.

def find_lpc(D,d,method):
    """
    Find linear prediction filter using a provided method
    """
    if method == "svd":
        return find_lpc_svd(D,d)
    elif method == "qr":
        return find_lpc_qr(D,d)
    elif method == "cholesky":
        return find_lpc_cholesky(D,d)
    elif method == "tls":
        return find_lpc_tls(D,d)
    else:
        raise ValueError("invalid method")


def find_lpc_svd(D,d):
    """
    find linear prediction filter using single value decomposition
    """
    L = D.shape[0]
    m = D.shape[1]
    U,s,Vh = scipy.linalg.svd(D)    # SVD decomposition
    U,Vh = np.mat(U),np.mat(Vh)     # make U and Vh matrices
    
    Si = pinv_diagsvd(s,m,L)    # construct the pseudo-inverse sigma matrix
    return np.array(Vh.H*Si*U.H*d)


# the next 3 lines and the pinv_diagsvd function were adapted from the
# scipy.linalg.pinv2 function - jjh
eps = np.finfo('float').eps
feps = np.finfo('single').eps
_array_precision = {'f': 0, 'd': 1, 'F': 0, 'D': 1}

def pinv_diagsvd(s,m,L):
    """
    Construct the pseudo-inverse of the sigma matrix from singular values
    """
    t = s.dtype.char
    cond =  {0: feps*1e3, 1: eps*1e6}[_array_precision[t]]
    cutoff = s[0]*cond
    
    Si = np.zeros((m,L),t)
    for i in range(len(s)):
        if s[i] > cutoff:
            Si[i,i] = 1.0/np.conj(s[i])
    return Si

def find_lpc_qr(D,d):
    """
    find linear prediction filter using QR decomposition
    """
    L = D.shape[0]
    m = D.shape[1]
    q,r  = scipy.linalg.qr(D)
    q,r = np.mat(q),np.mat(r)

    # SPEED
    # the next line is slow and the use of pinv2 should be avoided as 
    # pseudo inversion of r involves a computationally expensive SVD 
    # decomposition which is not needed.  Rather r*x = q.H*d should be 
    # solved for x using LAPACK's ZTRTRS function (or similar function with
    # different prefix).  This is not currently available in scipy/numpy and 
    # therefore is not used here.
    return scipy.linalg.pinv2(r)*q.H*d

def find_lpc_cholesky(D,d):
    """
    find linear prediction filter using a Cholesky decomposition
    """
    # form the normal equation (D.H*D)*a = D.H*d
    
    # SPEED
    # this can be improved by using the Hankel nature of D
    D = np.mat(D)
    DhD = np.mat(np.dot(D.H,D))
    Dhd = np.mat(np.dot(D.H,d))

    (c,lower) = scipy.linalg.cho_factor(DhD)    # Compute Cholesky decomp.
    return scipy.linalg.cho_solve((c,lower),Dhd)    # solve normal equation


def find_lpc_tls(D,d):
    """
    find linear prediction filter using the Total Least Squares method
    """
    m = D.shape[1]  # the order of the prediction
    E = np.append(D,d,axis=1)   # form the augmented data matrix
    U,s,Vh = scipy.linalg.svd(E)    # SVD decompositon of augmented matrix
    V = np.conj(Vh.T)               # Hermetian transpose
    return (-1./V[m,m]*V[:m,m]).reshape( (m,1) ) 


def find_lpc_fb(x,order,bad_roots,fix_mode,method):

    """ 
    Determind LP coefficient using forward-backward linear prediction.

    Averages LP coefficients generated from solving the forward and backward
    linear prediction equations after reversing the roots of characteristic
    polynomial of the backward solution.  Method is described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 100, 202-207.

    Description of parameters can be found in lp function.

    """

    # find forward LP coefficients
    D,d = make_Dd(x,order,'f')
    a = find_lpc(D,d,method)

    # stabilize roots if needed
    if bad_roots != None:
        poles = find_roots(a,'f')
        poles = fix_roots(poles,bad_roots,fix_mode)
        a = find_coeff(poles,'f')
    # store the forward coefficients
    forward_a = a.copy()

    # find the backwards LP coefficients
    D,d = make_Dd(x,order,'b')
    a = find_lpc(D,d,method)

    # find poles, reverse poles
    poles = find_roots(a,'b')
    poles = [1./pole for pole in poles]
    # stabilize roots if needed
    if bad_roots != None:
        poles = fix_roots(poles,bad_roots,fix_mode)
    # find the backward predicted, forward ordered coefficients
    backward_a = find_coeff(poles,'f')

    # average the forward and backward coefficients
    return (forward_a+backward_a)/2.

def find_lpc_bf(x,order,bad_roots,fix_mode,method):
    """ 
    Determind LP coefficient using backward-forward linear prediction.

    Averages LP coefficients generated from solving the forward and backward
    linear prediction equations after reversing the roots of characteristic
    polynomial of the forward solution.  Similar to method described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 100, 202-207.

    Description of parameters can be found in lp function.

    """

    # find backward LP coefficients
    D,d = make_Dd(x,order,'b')
    a = find_lpc(D,d,method)

    # stabilize roots if needed
    if bad_roots != None:
        poles = find_roots(a,'b')
        poles = fix_roots(poles,bad_roots,fix_mode)
        a = find_coeff(poles,'b')
    # store the forward coefficients
    backward_a = a.copy()

    # find the forward LP coefficients
    D,d = make_Dd(x,order,'f')
    a = find_lpc(D,d,method)

    # find poles, reverse poles
    poles = find_roots(a,'f')
    poles = [1./pole for pole in poles]
    # stabilize roots if needed
    if bad_roots != None:
        poles = fix_roots(poles,bad_roots,fix_mode)
    # find the forward predicted, backward ordered coefficients
    forward_a = find_coeff(poles,'b')

    # average the forward and backward coefficients
    return (forward_a+backward_a)/2.

#####################################
# root finding and fixing functions #
#####################################

def find_lproots_hsvd(x,M,K,mode,zmethod='sm'):
    """
    Find LP roots (poles) using the HSVD method

    Parameters: 

    x       1D trace of data, the FID.
    M       Length (M+1) of data matrix to form
    K       Reduced prediction order (number of signal roots)
            K < min(M+1,len(x)-M)
    mode    Mode to perform LP (f-forward,b-backward)
    zmethod Method used to find Z' (lstsq-least squares, sm-Sherman-Morrison)

    Returns array of signal roots (poles)
    

    Perform a HSVD linear prediction to determind signal roots (poles) as 
    described in:
    Barkhuijsen, DeBeer, and VanOrmondt, JMR, 1987, 73, 553 
    
    Parameters x, M and K are the same as those described in the above article.
    zmethod refer to the method used to calculate Z', either a least-squares 
    method (lstsq) can be used to solve U_b*Z'=U_t or the Sherman-Morrison
    formula (sm) can be used to avoid the full matrix inversion with equation
    [12] being used to find Z'. The Sherman-Morrison method should be faster 
    with similar precision.

    """
    # check parameters
    if mode not in ['f','b']:
        raise ValueError("mode must be 'f' or 'b'")
    if zmethod not in ['lstsq','sm']:
        raise ValueError("zmethod must be 'lstsq' or 'sm'")
    if K > min(M+1,len(x)-M):
        raise ValueError("K must be less than min(M+1,len(x)-M)")
    
    # form the data matrix X
    N = len(x) 
    L = N-M-1
    
    if mode == "f":
        X = scipy.linalg.hankel(x[:L+1],x[L:])
    else:
        # for backward LP we need to make the hankel matrix:
        # x_N-1 x_N-2 ... x_N-M-1
        # x_N-2 x_N-3 ... x_N-M-2
        # ...
        # x_M   x_M-1 ... x_0
        X = scipy.linalg.hankel(x[:M-1:-1],x[M::-1])

    # SVD of data matrix and truncation of U to form Uk
    U,s,Vh = scipy.linalg.svd(X)
    Uk = np.mat(U[:,:K])    # trucated U matrix of rank K
    Ub = Uk[:-1]            # Uk with bottom row removed
    Ut = Uk[1:]             # Uk with top row removed

    # calculate the Z' matrix

    if zmethod == 'lstsq':  # solve Ub*Z' = Ut using least-squares
        Zp,resid,rank,s = scipy.linalg.lstsq(Ub,Ut)
    else:
        # solve using equation [12]:
        # Z' = (Ek + (u*uh / (1-uh*u)) ) * Ub.H*Ut
        uh = Uk[-1] # bottom row of Uk
        u = uh.H
        Zp = (np.eye(K,dtype=u.dtype) + (u*uh/(1.-uh*u)))*Ub.H*Ut

    # diagonalization (find eigenvalues) of Z' to yield roots
    return scipy.linalg.eigvals(Zp)


def find_roots(a,mode="f"):
    """
    Find LP roots (poles) from a set of LP coefficients

    Parameters:
        
    a       LP coefficients
    mode    "f" or "b" depending on the ordering of a, "f" assumes the
            coefficients are stored m,m-1,...1; "b" are stores 1,2,...,m
            where m is the order of the Linear prediction.
    """
    if mode not in ['f','b']:
        raise ValueError("mode must be 'f'or 'b'")

    # STABILITY
    # the algorithm here is that used by numpy roots, build the companion
    # matrix and find it's eigenvalues.  These values should be polished for
    # better numerical accuracy.

    # np.roots expects a array, p, with coefficients
    # p[0] * x**n + p[1] * x**(n-1] + ... + p[n-1]*x + p[n]
    # in forward mode LP the coefficients are ordered m,m-1,...1
    # in backward mode LP the coefficient are ordered is 1,2,...,m
    # To calculate the roots, create a leading 1.0+0.0j and reverse if needed.
    p = np.empty( len(a)+1, dtype=a.dtype)
    p[0] = (1.0+0.0j)
    
    if mode =="f":      # reverse for forward LP
        p[1:] = -a.flat[::-1]
    else:   # backward LP
        p[1:] = -a.flat[:]
    return np.roots(p)

def find_coeff(poles,mode="f"):
    """
    Find LP coefficients from a set of LP roots (poles)
    
    Parmeters:
    
    poles   Array of LP roots
    mode    "f" or "b"

    """

    # STABILITY
    # the algorithm used here is numpy poly function which convolves
    # the roots to find the coefficients, the accuracy of this method
    # depends on the dtype of the poles parameter.
    if mode not in ['f','b']:
        raise ValueError("mode must be 'f'or 'b'")

    if mode == 'f': # reverse resulting coefficients
        return np.squeeze(-np.poly(poles)[:0:-1])
    else:   # keep coefficients as is
        return np.squeeze(-np.poly(poles)[1:])

def reverse_filter(a,mode):
    """
    Reverse a filter (ie change forward LP to backwards LP)
    """
    nmode = {'f':'b','b':'f'}[mode]
    # find roots, replace each root with 1/root, then recalculate filter
    return find_coeff( [1./pole for pole in find_roots(a,mode)],nmode)

def fix_roots(poles,fix_roots="incr",fix_mode="reflect"):
    """
    Fix (stabilize) LP roots 

    Parameters:
       
    poles       Array of LP roots/poles.
    fix_roots   "incr" or "decr" indicating which roots to remove (increasing
                or decreasing signals.
    fix_mode    Method to use when regularizing roots, "on" moves each bad
                root onto the unit circle, "reflect" reflects the root over
                the unit circle.
    """

    if fix_roots not in ["incr","decr"]:
        raise ValueError("fix_roots must be 'incr' or 'decr'")
    if fix_mode not in ["on","reflect"]:
        raise ValueError("fix_mode must be 'on' or 'reflect'")

    if fix_roots == "incr":     # remove increasing signals
        for i,pole in enumerate(poles):
            if np.abs(pole) > 1:
                #print "Fixing root:",i
                if fix_mode == "on":
                    poles[i] = pole/np.abs(pole)
                else:
                    poles[i] = 1/np.conj(pole)
    else:   # remove decreasing signals
        for i,pole in enumerate(poles):
            if np.abs(pole) < 1:
                #print "Fixing root:",i
                if fix_mode == "on":
                    poles[i] = pole/np.abs(pole)
                else:
                    poles[i] = 1/np.conj(pole)
        
    return poles


###########################
# Extrapolation functions #
###########################

def extrapolate(trace,a,pred,append):
    """
    Extrapolate points using LP prediction filter

    Parameters:
    
    trace   1D vector to add points to
    a       LP coefficients (ordered according to direction of extrapolation)
    pred    Number of points to predict using LP prediction filter
    append  "b" or "a" indicating if extrapolation should append before
            or after known points.

    """
    m = len(a)      # LP order
    M = len(trace)  # number of points in original trace
    ntrace = np.empty((M+pred),dtype=trace.dtype)

    if append not in ["after","before"]:
        raise ValueError("append must be 'a' or 'b'")

    if append == "after":   # append after trace
        ntrace[:M] = trace
        for i in xrange(pred):
            ntrace[M+i] = np.sum(np.multiply(ntrace[M-m+i:M+i],a.flat))
        return ntrace

    if append == "before":   # append before trace
        ntrace[-M:] = trace
        for i in xrange(pred):
            ntrace[pred-i-1]=np.sum(np.multiply(ntrace[pred-i:pred+m-i],a.flat))
        return ntrace
