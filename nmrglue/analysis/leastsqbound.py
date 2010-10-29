"""
Constrained multivariate Levenberg-Marquardt optimization
"""

from scipy.optimize import leastsq
import numpy as np

def internal2external_grad(xi,bounds):
    """ 
    Calculate the internal to external gradiant
    
    Calculates the partial of external over internal
    
    """
    
    ge = np.empty_like(xi)

    for i,(v,bound) in enumerate(zip(xi,bounds)):
        
        a = bound[0]    # minimum
        b = bound[1]    # maximum

        if a == None and b == None:    # No constraints
            ge[i] = 1.0

        elif b == None:      # only min
            ge[i] = v/np.sqrt(v**2+1)

        elif a == None:      # only max
            ge[i] = -v/np.sqrt(v**2+1)

        else:       # both min and max
            ge[i] = (b-a)*np.cos(v)/2.

    return ge

def i2e_cov_x(xi,bounds,cov_x):

    grad = internal2external_grad(xi,bounds)
    grad = grad = np.atleast_2d(grad)
    return np.dot(grad.T,grad)*cov_x


def internal2external(xi,bounds):
    """ Convert a series of internal variables to external variables"""
    
    xe = np.empty_like(xi)

    for i,(v,bound) in enumerate(zip(xi,bounds)):
        
        a = bound[0]    # minimum
        b = bound[1]    # maximum

        if a == None and b == None:    # No constraints
            xe[i] = v

        elif b == None:      # only min
            xe[i] = a-1.+np.sqrt(v**2.+1.)

        elif a == None:      # only max
            xe[i] = b+1.-np.sqrt(v**2.+1.)

        else:       # both min and max
            xe[i] = a+((b-a)/2.)*( np.sin(v)+1.)

    return xe

def external2internal(xe,bounds):
    """ Convert a series of external variables to internal variables"""

    xi = np.empty_like(xe)

    for i,(v,bound) in enumerate(zip(xe,bounds)):
        
        a = bound[0]    # minimum
        b = bound[1]    # maximum

        if a == None and b == None: # No constraints
            xi[i] = v

        elif b == None:     # only min
            xi[i] = np.sqrt( (v-a+1.)**2.-1 )

        elif a == None:     # only max
            xi[i] = np.sqrt( (b-v+1.)**2.-1 )

        else:   # both min and max
            xi[i] = np.arcsin( (2.*(v-a)/(b-a))-1.)

    return xi

def err(p,bounds,efunc,args):
    
    pe = internal2external(p,bounds)    # convert to external variables
    return efunc(pe,*args)

def calc_cov_x(infodic,p):
    """
    Calculate cov_x from fjac, ipvt and p as is done in leastsq
    """

    fjac = infodic['fjac']
    ipvt = infodic['ipvt']
    n = len(p)

    # adapted from leastsq function in scipy/optimize/minpack.py
    perm = np.take(np.eye(n),ipvt-1,0)
    r = np.triu(np.transpose(fjac)[:n,:])
    R = np.dot(r,perm)
    try:
        cov_x = np.linalg.inv(np.dot(np.transpose(R),R))
    except LinAlgError:
        cov_x = None
    return cov_x


def leastsqbound(func,x0,bounds,args=(),**kw):
    """
    Constrained multivariant Levenberg-Marquard optimization

    Minimize the sum of squares of a given function using the 
    Levenberg-Marquard algorithm. Contraints on parameters are inforced using 
    variable transformations as described in the MINUIT User's Guide by
    Fred James and Matthias Winkler.

    Parameters:

    * func      functions to call for optimization.
    * x0        Starting estimate for the minimization.
    * bounds    (min,max) pair for each element of x, defining the bounds on
                that parameter.  Use None for one of min or max when there is
                no bound in that direction.
    * args      Any extra arguments to func are places in this tuple.

    Returns: (x,{cov_x,infodict,mesg},ier)

    Return is described in the scipy.optimize.leastsq function.  x and con_v  
    are corrected to take into account the parameter transformation, infodic 
    is not corrected.

    Additional keyword arguments are passed directly to the 
    scipy.optimize.leastsq algorithm. 

    """
    # check for full output
    if "full_output" in kw and kw["full_output"]:
        full=True
    else:
        full=False

    # convert x0 to internal variables
    i0 = external2internal(x0,bounds)

    # perfrom unconstrained optimization using internal variables
    r = leastsq(err,i0,args=(bounds,func,args),**kw)

    # unpack return convert to external variables and return
    if full:
        xi,cov_xi,infodic,mesg,ier = r
        xe = internal2external(xi,bounds)
        cov_xe = i2e_cov_x(xi,bounds,cov_xi)
        # XXX correct infodic 'fjac','ipvt', and 'qtf' 
        return xe,cov_xe,infodic,mesg,ier 

    else:
        xi,ier = r
        xe = internal2external(xi,bounds)
        return xe,ier


