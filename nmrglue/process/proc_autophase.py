"""
Automated phase correction
These functions provide support for automatic phasing of NMR data. They
consist of the core `autops` function which performs the optimisation and
a set of private functions for calculating a spectral phase quality score
for a provided spectrum.
"""

from __future__ import print_function

import numpy as np
import scipy.optimize

from .proc_base import ps


def psacme(data, p0=0.0, p1=0.0, gamma=5.0e-3, optreturn=False):
    """
    Phase correction using ACME algorithm by Chen Li et al.
    Journal of Magnetic Resonance 158 (2002) 164-168

    Parameters
    ----------
    ph : tuple
        Current p0 and p1 values
    data : ndarray
        Array of NMR data.
    gamma : float
        scale factor set to balance the penalty and entropy
    optreturn : boolean
        option to return the optimum p0 and p1

    Returns
    -------
    score : float
        Value of the objective function (phase score)
    opt : tuple
        phases returned by algoritm in a tuple

    """

    opt = [p0, p1]
    opt = scipy.optimize.fmin(_ps_acme_score, x0=opt, args=(data, gamma))

    phasedspc = ps(data, p0=opt[0], p1=opt[1])
    
    if optreturn:
        return tuple(opt), phasedspc
    else:
        return phasedspc
  

def pspmin( data, p0=0.0, p1=0.0, optreturn=False):
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
    opt : tuple
        phases returned by algoritm in a tuple

    """

    opt = [p0, p1]
    opt = scipy.optimize.fmin(_ps_acme_score, x0=opt, args=(data,))

    phasedspc = ps(data, p0=opt[0], p1=opt[1])
    
    if optreturn:
        return tuple(opt), phasedspc
    else:
        return phasedspc
    

def _ps_acme_score(ph, data, gamma):
    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    # Calculation of first derivatives
    ds1 = np.abs(np.diff(data, 1))
    p1 = ds1 / np.sum(ds1)

    # Calculation of entropy
    h1 = -p1 * np.log(p1)
    h1s = np.sum(h1)

    # Calculation of penalty
    fr = p1
    fr[fr >= 0] = 0
    fr[fr < 0] = 1
    pr = gamma * np.sum(fr * p1**2)

    return h1s + pr


def _ps_peak_minima_score(ph, data):
    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    i = np.argmax(data)
    mina = np.min(data[i-100:i])
    minb = np.min(data[i:i+100])

    return np.abs(mina - minb)


def manual_ps(data):
    """
    Manual Phase correction using matplotlib

    A matplotlib widget is used to manually correct the phase of a Fourier
    transformed dataset. If the dataset has more than 1 dimensions, the first
    trace will be picked up for phase correction.  Clicking the 'Set Phase'
    button will print the current linear phase parameters to the console.

    .. note:: Needs matplotlib with and interactive backend.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    p0, p1 : float
        Linear phase correction parameters. Zero and first order phase
        corrections in degrees calculated from pc0, pc1 and pivot displayed
        in the interactive window.

    Examples
    --------
    >>> import nmrglue as ng
    >>> p0, p1 = ng.process.proc_autophase.manual_ps(data)
    >>> # do manual phase correction and close window
    >>> phased_data = ng.proc_base.ps(data, p0=p0, p1=p1)

    """

    from matplotlib.widgets import Slider, Button
    import matplotlib.pyplot as plt

    plt.subplots_adjust(left=0.25, bottom=0.35)

    if len(data.shape) == 2:
        data = data[0, ...]
    elif len(data.shape) == 3:
        data = data[0, 0, ...]
    elif len(data.shape) == 4:
        data = data[0, 0, 0, ...]

    interactive, = plt.plot(data.real, lw=1, color='black')

    axcolor = 'white'
    axpc0 = plt.axes([0.25, 0.10, 0.65, 0.03], axisbg=axcolor)
    axpc1 = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axpiv = plt.axes([0.25, 0.20, 0.65, 0.03], axisbg=axcolor)
    axpst = plt.axes([0.25, 0.25, 0.15, 0.04], axisbg=axcolor)

    spc0 = Slider(axpc0, 'p0', -360, 360, valinit=0)
    spc1 = Slider(axpc1, 'p1', -360, 360, valinit=0)
    spiv = Slider(axpiv, 'pivot', 0, data.size, valinit=0)
    axps = Button(axpst, 'Set Phase', color=axcolor)

    def update(val):
        pc0 = spc0.val * np.pi / 180
        pc1 = spc1.val * np.pi / 180
        pivot = spiv.val
        interactive.set_ydata((data * np.exp(
            1.0j * (pc0 + (pc1 * np.arange(-pivot, -pivot + data.size) /
                    data.size))).astype(data.dtype)).real)
        plt.draw()

    def setphase(val):
        p0 = spc0.val-spc1.val*spiv.val/data.size
        p1 = spc1.val
        print(p0, p1)

    spc0.on_changed(update)
    spc1.on_changed(update)
    spiv.on_changed(update)
    axps.on_clicked(setphase)

    plt.show(block=True)

    p0 = spc0.val-spc1.val*spiv.val/data.size
    p1 = spc1.val
    return p0, p1
