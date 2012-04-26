"""
lineshapes1d.py One-dimensional lineshape functions and classes
"""

import numpy as np
import scipy.special    # needed for complex error function (wo
pi = np.pi

####################################
# 1D lineshape simulator functions #
####################################

# Gaussian (normal) lineshape simulator functions
def sim_gauss_sigma(x, x0, sigma):
    """
    Simulate a Gaussian (normal) lineshape with unit height at the center. 
    
    Simulate discrete points of a continuous Gaussian (normal) distribution
    with unit height at the center.  Sigma (the standard deviation of the
    distribution) is used as the distribution scale parameter. 

    Functional form:

        f(x; x0, sigma) = exp(-(x-x0)^2  / (2 * sigma^2))

    Parameters:

    * x     Array of values at which to evaluate distribution.
    * x0    Center (mean) of Gaussian distribution.
    * sigma Scale (variance) of the Gaussian distribution.

    """
    return np.exp(-(x-x0)**2 / (2.0*sigma**2) )


def sim_gauss_fwhm(x, x0, fwhm):
    """
    Simulate a Gaussian (normal) lineshape with unit height at the center. 
    
    Simulate discrete points of a continuous Gaussian (normal) distribution
    with unit height at the center.  FWHM (full-width at half-maximum ) is used
    as the distribution scale parameter.

    Functional form:

        f(x; x0, fwhm) = exp( -(x-x0)^2 * 4 * ln(2)  / (fwhm^2))

    Parameters:

    * x     Array of values at which to evaluate distribution.
    * x0    Center (mean) of Gaussian distribution.
    * fwhm  Full-width at half-maximum of distribution.
    
    """
    return np.exp(-(x-x0)**2 * 4 * np.log(2) / (fwhm**2))

# Lorentzian lineshape simulator functions
def sim_lorentz_gamma(x, x0, gamma):
    """ 
    Simulate a Lorentzian lineshape with unit height at the center.

    Simulates discrete points of the continuous Cauchy-Lorentz (Breit-Wigner)
    distribution with unit height at the center.  Gamma (the half-width at
    half-maximum, HWHM) is used as the scale parameter.

    Functional form:

        f(x; x0, gamma) = g^2 / ((x-x0)^2 + g^2) 
    
    Parameters:
        
    * x     Array of values at which to evaluate distribution.
    * x0    Center of the distribution.
    * gamma Scale parameter, half-width at half-maximum, of distribution.
    
    """
    return gamma**2 / (gamma**2 + (x-x0)**2)

def sim_lorentz_fwhm(x, x0, fwhm):
    """ 
    Simulate a Lorentzian lineshape with unit height at the center.

    Simulates discrete points of the continuous Cauchy-Lorentz (Breit-Wigner)
    distribution with unit height at the center.  FWHM (full-width at
    half-maximum) is used as the scale parameter.

    Functional form:

        f(x; x0, fwhm) = (0.5 * fwhm)^2 / ((x-x0)^2 + (0.5 * fwhm)^2) 
    
    Parameters:
        
    * x     Array of values at which to evaluate distribution.
    * x0    Center of the distribution.
    * fwhm  Full-width at half-maximum of distribution.
    
    """
    return (0.5 * fwhm)**2 / ((0.5 * fwhm)**2 + (x-x0)**2)

# Voigt lineshape simulator functions

def sim_voigt_fwhm(x, x0, fwhm_g, fwhm_l): 
    """
    Simulate a Voigt lineshape with unit height at the center.

    Simulates discrete points of the continuous Voigt profile with unit height
    at the center.  Full-width at half-maximum (FWHM) of each component are used
    as the scale parameters for the Gaussian and Lorentzian distribution.

    Functional Form:

        f(x; x0, fwhm_g, fwhm_l) = Re[w(z)] / Re[(w(z0)]

    Where:
        
        z = sqrt(ln(2)) * (2 * (x - x0) + 1j * fwhm_l) / fwhm_g
        z0 = sqrt(ln(2)) * 1j * fwhm_l / fwhm_g
        w(z) is the complex error function of z

    Parameters

    * x         Array of values at which to evalutate distribution.
    * x0        Center of the distribution.
    * fwhm_g    Full-width at half-maximum of the Gaussian component.
    * fwhm_l    Full-width at half-maximum of the Lorentzian component.

    """
    z = np.sqrt(np.log(2)) * (2.0 * (x - x0) + 1.j * fwhm_l) / fwhm_g
    z0 = np.sqrt(np.log(2)) * 1.j *  fwhm_l / fwhm_g
    return scipy.special.wofz(z).real / scipy.special.wofz(z0).real


def sim_voigt_sigmagamma(x, x0, sigma, gamma):
    """
    Simulate a Voigt lineshape with unit height at the center.

    Simulates discrete points of the continuous Voigt profile with unit height
    at the center.  Sigma and gamma are used as the Gaussian and Lorentzian
    scaler parameters.

    Functional Form:

        f(x; x0, sigma, gamma) = Re[w(z)] / Re[(w(z0)]

    Where:
        
        z = ((x - x0) + 1j * gamma) / (sigma * sqrt(2))
        z0 = (1j * gamma) / (sigma * sqrt(2))
        w(z) is the complex error function of z

    Parameters

    * x     Array of values at which to evalutate distribution.
    * x0    Center of the distribution
    * sigma Gaussian scale component of Voigt profile.  Variance of the Gaussian
            distribution.
    * gamma Lorentzian scale component of Voigt profile.  Half-width at
            half-maximum of the Lorentzian component. 

    """
    z = (x - x0 + 1j * gamma) / ( sigma * np.sqrt(2))
    z0 = (1j * gamma) / (sigma * np.sqrt(2))
    return scipy.special.wofz(z).real / scipy.special.wofz(z0).real


# Pseudo Voigt linehspae simulator functions
def sim_pvoigt_fwhm(x, x0, fwhm, eta): 
    """
    Simulate a Pseudo Voigt lineshape with unit height at the center.

    Simulates discrete points of the continuous Pseudo Voigt profile with unit 
    heigh at the center.  Full-width at half-maximum (FWHM) of the Gaussian and
    Lorentzian distribution are used as the scale parameter as well as eta, the
    mixing factor.  

    Functional Form:

        f(x; x0, fwhm, eta) = (1-eta) * G(x; x0, fwhm) + eta * L(x; x0, fwhm)

    Where:
        
        G(x; x0, fwhm) = exp( -(x-x0)^2 * 4 * ln(2)  / (fwhm^2))
        L(x; x0, fwhm) = (0.5 * fwhm)^2 / ((x-x0)^2 + (0.5 * fwhm)^2) 

    Parameters

    * x     Array of values at which to evalutate distribution.
    * x0    Center of the distribution.
    * fwhm  Full-width at half-maximum of the Pseudo Voigt profile.
    * eta   Lorentzian/Gaussian mixing parameter.

    """
    G = sim_gauss_fwhm(x, x0, fwhm)
    L = sim_lorentz_fwhm(x, x0, fwhm)
    return (1.0 - eta) * G + eta * L


########################
# 1D Lineshape classes #
########################

# A lineshape class defines methods used to fit and simulate one dimensional
# lineshapes, which can be used to build multidimensinal lineshapes.  These
# classes should have the following 6 methods:

# sim(self, M, p)   - Using parameters in p simulate a lineshape of length M.
# nparams(self, M)  - Determind the number of parameters needed for a length M 
#                     lineshape.
# guessp(self, sig) - Estimate parameters of signal sig, these should be 
#                     parameter which might be used for initial least-squares 
#                     fitting.
# pnames(self, M)   - Give names to the parameters of a lineshape of length M.
#
# add_edge(self, p, (min,max))  - take into account region limits at min,max
#                                 for parameters in p.
# remove_edge(self, p, (min,max))   - remove the effects of region limits 
#                                   min, max for parameters in p.


# location-scale lineshapes

class location_scale():
    """
    Class for building a 2 parameter location scale lineshape class
    """
    def nparam(self,M):
        return 2

    def add_edge(self,p,(min,max)):
        if p[0] is None:
            return p
        return p[0]-min,p[1]

    def remove_edge(self,p,(min,max)):
        if p[0] is None:
            return p
        return p[0]+min,p[1]

class gauss_sigma(location_scale):
    """
    Gaussian (normal) lineshape class with unit height at the mean and sigma
    scale parameter. See sim_gauss_sigma for functional form and parameters. 
    """
    name = "guassian"

    def sim(self,M,p):
        x = np.arange(M)
        x0, sigma = p
        return sim_gauss_sigma(x, x0, sigma)

    def guessp(self,sig):
        c,fwhm = center_fwhm(sig)
        return (c,fwhm/2.35482004503)

    def pnames(self,M):
        return ("x0","sigma")

class gauss_fwhm(location_scale):
    """
    Gaussian (normal) lineshape class with unit height at the mean and fwhm
    scale parameter. See sim_gauss_fwhm for functional form and parameters. 
    """
    name = "guassian"

    def sim(self,M,p):
        x = np.arange(M)
        x0, fwhm = p
        return sim_gauss_fwhm(x, x0, fwhm)

    def guessp(self,sig):
        c,fwhm = center_fwhm(sig)
        return (c,fwhm)

    def pnames(self,M):
        return ("x0","fwhm")

class lorentz_gamma(location_scale):
    """
    Lorentzian lineshape class with unit height at the center and gamma scale
    parameter.  See sim_lorentz_gamma for functional form and parameters.
    """
    name = "lorentz"

    def sim(self,M,p):
        x = np.arange(M)
        x0, gamma = p
        return sim_lorentz_gamma(x, x0, gamma)

    def guessp(self,sig):
        c,fwhm = center_fwhm(sig)
        return (c,fwhm/2.)

    def pnames(self,M):
        return("x0","gamma")

class lorentz_fwhm(location_scale):
    """
    Lorentzian lineshape class with unit height at the center and gamma scale
    parameter.  See sim_lorentz_fwhm for functional form and parameters.
    """
    name = "lorentz"

    def sim(self,M,p):
        x = np.arange(M)
        x0, fwhm = p
        return sim_lorentz_fwhm(x, x0, fwhm)

    def guessp(self,sig):
        c,fwhm = center_fwhm(sig)
        return (c,fwhm)

    def pnames(self,M):
        return("x0","fwhm")

# Voigt (location, 2 scale-like parameters) lineshapes. 

class location_2params():
    """
    Class for building a 3 parameter location, scale, other lineshape class
    """
    def nparam(self,M):
        return 3

    def add_edge(self,p,(min,max)):
        if p[0] is None:
            return p
        return p[0]-min,p[1]

    def remove_edge(self,p,(min,max)):
        if p[0] is None:
            return p
        return p[0]+min,p[1]

class voigt_fwhm(location_2params):
    """
    Voigt lineshape class with unit height at the center and full-width 
    half-maximum scale parameters.  See sim_voigt_fwhm for functional form and
    parameters.
    """
    name = "voigt"

    def sim(self, M, p):
        x = np.arange(M)
        x0, fwhm_g, fwhm_l = p
        return sim_voigt_fwhm(x, x0, fwhm_g, fwhm_l)
    
    def guessp(self, sig):
        c, fwhm = center_fwhm(sig)
        return (c, fwhm*0.5, fwhm*0.5)

    def pnames(self, M):
        return ("x0", "fwhm_gauss", "fwhm_lorentz")

class voigt_sigmagamma(location_2params):
    """
    Voigt lineshape class with unit height at the center and sigma, gamma scale
    parameters.  See sim_voigt_sigmagamma for functional form and parameters.
    """
    name = "voigt"

    def sim(self, M, p):
        x = np.arange(M)
        x0, sigma, gamma = p
        return sim_voigt_sigmagamma(x, x0, sigma, gamma)
    
    def guessp(self, sig):
        c, fwhm = center_fwhm(sig)
        return (c, fwhm/2.35482004503 * 0.5  , fwhm * 0.5 * 0.5)

    def pnames(self, M):
        return ("x0", "fwhm_gauss", "fwhm_lorentz")

class pvoigt_fwhm(location_2params):
    """
    Pseudo-Voigt lineshape class with unit height at the center and full-width 
    half-maximum scale parameter.  See sim_pvoigt_fwhm for functional form and
    parameters.
    """
    name = "pvoigt"

    def sim(self, M, p):
        x = np.arange(M)
        x0, fwhm, eta = p
        return sim_pvoigt_fwhm(x, x0, fwhm, eta)
    
    def guessp(self, sig):
        c, fwhm = center_fwhm(sig)
        return (c, fwhm, 0.5)

    def pnames(self, M):
        return ("x0", "fwhm", "eta")


# misc lineshape classes
class scale():
    """
    Scale lineshape class

    Simulates a lineshape with functional form:

    1.0, a0, a1, a2, ....

    Where a0, a1, ... are the parameters provided.

    """
    name = "scale"
    
    def sim(self,M,p):
        l = np.empty(M,dtype='float')
        l[0] = 1
        l[1:] = p
        return l

    def nparam(self,M):
        return int(M-1)

    def guessp(self,sig):
        return sig[1:]/sig[0]

    def pnames(self,M):
        return tuple(["a%i"%i for i in range(1,M)])
    
    def add_edge(self,p,(min,max)):
        return p

    def remove_edge(self,p,(min,max)):
        return p

# lineshape convience 
gauss = gauss_fwhm
lorentz = lorentz_fwhm
voigt = voigt_fwhm
pvoigt = pvoigt_fwhm

# lineshape class router

def ls_str2class(l):
    """ Convert lineshape string to lineshape class """
    if l == "gauss" or l == "g":
        return gauss()
    elif l == "lorentz" or l == "l":
        return lorentz()
    elif l == "scale" or l == "s":
        return scale()
    elif l == "voigt" or l == "v":
        return voigt()
    elif l == "pvoigt" or l == "pv":
        return pvoigt()
    else:
        raise ValueError("Unknown lineshape %s",(l))

# basic lineshape analysis

def center_fwhm(signal):
    """
    Estimate the center and full-width half max of a signal.
    """

    # negate the signal if it appears to be a negative peak
    if -signal.min() > signal.max():
        signal = -signal

    # the center is the highest point in the signal 
    center = signal.argmax()

    # find the points that bracket the first and last crossing of the
    # half max, then use linear extrapolation to find the location of the
    # half max on either side of the maximum.  The difference between these
    # two values is a good approximation of the full width at half max
    max = signal.max()
    hmax = max/2.

    top_args = np.nonzero(signal > hmax)[0]     # all points above half-max
    l_idx = top_args[0]     # index of left hand side above half-max
    r_idx = top_args[-1]    # index of right hand side above half-max

    # solve hmax = mx+b => x = y-b/m
    # for two points x_0 and x_1 this becomes y = (hmax-x_0)/(x_1-x_0)
    # to this value add the index of x_0 to get the location of the half-max.
    
    # left side
    if l_idx==0:
        left = l_idx    # this is a bad guess but the best we can do
    else:
        x_0,x_1 = signal[l_idx-1],signal[l_idx]
        left = l_idx-1+(hmax-x_0)/(x_1-x_0)

    # right side
    if r_idx == len(signal)-1:
        right = r_idx   # this is a bad guess but the best we can do
    else:
        x_0,x_1 = signal[r_idx],signal[r_idx+1]
        right = r_idx+(hmax-x_0)/(x_1-x_0)

    return center,right-left


def center_fwhm_bymoments(signal):
    """
    Estimate the center and full-width half max of a signal using moments
    """
    
    # calculate the zeroth, first and second moment
    x = np.arange(signal.size)
    m0 = signal.sum()
    m1 = (x*signal).sum()
    m2 = (x**2.*signal).sum()

    # mu (the center) is the first over the zeroth moment
    mu = m1/m0
    # sigma (the variance) is  sqrt( abs(m2)-mu**2)
    sigma = np.sqrt(np.abs(m2)-mu**2)

    return mu,sigma*2.3548
