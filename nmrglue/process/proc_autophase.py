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


def autops(data, fn, p0=0.0, p1=0.0, return_phases=False, **kwargs):
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
    opt = scipy.optimize.fmin(fn, x0=opt, args=(data,), **kwargs)

    phasedspc = ps(data, p0=opt[0], p1=opt[1])

    if return_phases:
        return phasedspc, opt
    else:
        return phasedspc


def _ps_acme_score(ph, data):
    """
    Phase correction using ACME algorithm by Chen Li et al.
    Journal of Magnetic Resonance 158 (2002) 164-168

    Parameters
    ----------
    ph : tuple
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

    return (h1s + p) / data.shape[-1] / np.max(data)


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


def manual_ps(data, notebook=False):
    """
    Manual Phase correction using matplotlib

    A matplotlib widget is used to manually correct the phase of a Fourier
    transformed dataset. If the dataset has more than 1 dimensions, the first
    trace will be picked up for phase correction.  Clicking the 'Set Phase'
    button will print the current linear phase parameters to the console.
    A ipywidget is provided for use with Jupyter Notebook to avoid changing
    backends. This can be accessed with notebook=True option in this function

    .. note:: Needs matplotlib with an interactive backend.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    notebook : Bool
        True for plotting interactively in Jupyter Notebook
        Uses ipywidgets instead of matplotlib widgets

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

    In  [1] # if you are using the Jupyter Notebook
    In  [2] ng.process.proc_autophase.manual_ps(data)
    Out [2] # do manual phase correction. p0 and p1 values will be updated
            # continuously as you do so and are printed below the plot
    In  [3] phased_data = ng.proc_base.ps(data, p0=p0, p1=p1)

    """

    if len(data.shape) == 2:
        data = data[0, ...]
    elif len(data.shape) == 3:
        data = data[0, 0, ...]
    elif len(data.shape) == 4:
        data = data[0, 0, 0, ...]

    if notebook:
        from ipywidgets import interact, fixed
        import matplotlib.pyplot as plt

        def phasecorr(dataset, phcorr0, phcorr1, pivot):
            fig, ax = plt.subplots(figsize=(10, 7))
            phaseddata = dataset * np.exp(
                1j * (phcorr0 + phcorr1 * (
                    np.arange(-pivot, -pivot+dataset.size)/dataset.size)))

            ax.plot(np.real(phaseddata))
            ax.set(ylim=(np.min(np.real(data))*2, np.max(np.real(data))*2))
            ax.axvline(pivot, color='r', alpha=0.5)
            plt.show()

            p0 = np.round(
                (phcorr0 - phcorr1 * pivot/dataset.size) * 360 / 2 / np.pi, 3)
            p1 = np.round(phcorr1*360/2/np.pi, 3)

            print('p0 =', p0, 'p1 =', p1)

        interact(
            phasecorr,
            dataset=fixed(data),
            phcorr0=(-np.pi, np.pi, 0.01),
            phcorr1=(-10*np.pi, 10*np.pi, 0.01),
            pivot=(0, data.size, 1))

    else:

        from matplotlib.widgets import Slider, Button
        import matplotlib.pyplot as plt

        plt.subplots_adjust(left=0.25, bottom=0.35)

        interactive, = plt.plot(data.real, lw=1, color='black')

        axcolor = 'white'
        axpc0 = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor=axcolor)
        axpc1 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axpiv = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
        axpst = plt.axes([0.25, 0.25, 0.15, 0.04], facecolor=axcolor)

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
