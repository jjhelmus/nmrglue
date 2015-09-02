"""
Automated phase correction
These functions provide support for automatic phasing of NMR data. They
consist of the core `autops` function which performs the optimisation and
a set of private functions for calculating a spectral phase quality score
for a provided spectrum.
"""

import numpy as np
import scipy.optimize
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt

from .proc_base import ps


def autops(data, fn, p0=0.0, p1=0.0):
    """
    Automatic linear phase correction

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    fn : str or function
        Algorithm to use for phase scoring. Built in functions can be
        specified by one of the following strings: "acme", "peak_minima"
    p0 : float
        Initial zero order phase in degrees.
    p1 : float
        Initial first order phase in degrees.

    Returns
    -------
    ndata : ndarray
        Phased NMR data.

    """
    if not callable(fn):
        fn = {
            'peak_minima': _ps_peak_minima_score,
            'acme': _ps_acme_score,
        }[fn]

    opt = [p0, p1]
    opt = scipy.optimize.fmin(fn, x0=opt, args=(data, ))

    phasedspc = ps(data, p0=opt[0], p1=opt[1])

    return phasedspc


def _ps_acme_score(ph, data):
    """
    Phase correction using ACME algorithm by Chen Li et al.
    Journal of Magnetic Resonance 158 (2002) 164-168

    Parameters
    ----------
    pd : tuple
        Current p0 and p1 values
    data : ndarray
        Array of NMR data.

    Returns
    -------
    score : float
        Value of the objective function (phase score)

    """
    stepsize = 1

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    # Calculation of first derivatives
    ds1 = np.abs((data[1:]-data[:-1]) / (stepsize*2))
    p1 = ds1 / np.sum(ds1)

    # Calculation of entropy
    p1[p1 == 0] = 1

    h1 = -p1 * np.log(p1)
    h1s = np.sum(h1)

    # Calculation of penalty
    pfun = 0.0
    as_ = data - np.abs(data)
    sumas = np.sum(as_)

    if sumas < 0:
        pfun = pfun + np.sum((as_/2) ** 2)

    p = 1000 * pfun

    return h1s + p


def _ps_peak_minima_score(ph, data):
    """
    Phase correction using simple minima-minimisation around highest peak

    This is a naive approach but is quick and often achieves reasonable
    results.  The optimisation is performed by finding the highest peak in the
    spectra (e.g. TMSP) and then attempting to reduce minima surrounding it.

    Parameters
    ----------
    pd : tuple
        Current p0 and p1 values
    data : ndarray
        Array of NMR data.

    Returns
    -------
    score : float
        Value of the objective function (phase score)

    """

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    i = np.argmax(data)
    mina = np.min(data[i-100:i])
    minb = np.min(data[i:i+100])

    return np.abs(mina - minb)
    

def manual_ps(data):
    '''
    Interactive Phase correction for 1-D and 2-D datasets
    
    Usage:
    manual_ps(data)
    
    Needs:
    Matplotlib with an interactive backend (eg. QT)
    
    Parameters
    ----------
    data : ndarray for a fourier transformed dataset
           if array dimansions > 1
           data for the first dimension is used
           
           
    Other internal parameters used in interactive phase correction:
    -----------------------------------------------------
    pc0 = 1st order phase correction
    pc1 = 1st order phase correction
    piv = pivot point 

    While using other phase corrections modules, eg. proc_base.ps(),
    these parameters translate to p0 and p1 as follows:

    p0 = pc0 - pc1 * piv
    p1 = pc1
    

   Returns
   -------
   A global variable phcorr (tuple), 
   with phcorr[0] = p0 and phcorr[1] = p1


    How to use:
    -----------
    1. put matplotlib to an interactive backend
       %matplotlib qt (if you are using ipython)
    2. manual_ps(data)
    3. Use phcorr[0] and phcorr[1] as 0th and 1st order phase
       corrections in phase correction programs
       eg. phase_corrected_data = ng.proc_base.ps(data, p0=phcorr[0], 
                                                        p1=phcorr[1])      
    
    '''
     

    plt.subplots_adjust(left=0.25, bottom=0.30)
     
    
    if len(data.shape) > 1:
        data = data[0]
    
    global phcorr
    phcorr = (0, 0) # (p0, p1)
    
    interactive, = plt.plot(data, lw=1, color='black')
 
    axcolor = 'white'
    axpc0 = plt.axes([0.25, 0.10, 0.65, 0.03], axisbg=axcolor)
    axpc1 = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axpiv = plt.axes([0.25, 0.20, 0.65, 0.03], axisbg=axcolor)

    spc0 = Slider(axpc0, 'p0', -360, 360, valinit=0)
    spc1 = Slider(axpc1, 'p1', -360, 360, valinit=0)
    spiv = Slider(axpiv, 'pivot', 0, data.size, valinit=0)

    def update(val):
        pc0 = spc0.val * np.pi / 180
        pc1 = spc1.val * np.pi / 180
        pivot = spiv.val
        interactive.set_ydata(data * np.exp(1.0j *  
             (pc0 + (pc1 * np.arange(-pivot, -pivot + data.size) / data.size)))
             .astype(data.dtype))
        plt.draw()     
        
        global phcorr
        phcorr = ( spc0.val-spc1.val*spiv.val, spc1.val  ) 

    spc0.on_changed(update)
    spc1.on_changed(update)
    spiv.on_changed(update)        

    plt.show()
    
    
