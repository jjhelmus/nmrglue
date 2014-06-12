"""
NMRPipe like processing functions for use with the
:py:mod:`nmrglue.fileio.pipe` module.

These functions attempt to mimic NMRPipe's processing functions but small
differences exist between to two implementations.  In particular when using
this module:

    * hdr=True overrides all values in the calling function.
    * A di flag is not used, rather the :py:func:`di` function should be used
      to delete the imaginary portion of a spectra.
    * x1, xn and other limits must be expressed in points. A unit conversion
      object function should be used before calling the processing function to
      calculate these values.
    * No functions implement the dmx or nodmx flags.

Additional differences from NMRPipe's functions are documented in the
individual processing functions.

The following functions have not been implemented and will raise a
NotImplemented exception:

    * ann      Fourier Analysis by Neural Net
    * ebs      EBS Reconstruction
    * mac      Macro Language Interpreter
    * mem      Maximum Entropy
    * ml       Maximum likelyhood frequency
    * poly     Polynomail baseline correction
    * xyz2zyx  3D matrix transpose
    * ztp      3D matrix transpose

"""

import numpy as np

# nmrglue modules
from ..fileio import pipe, fileiobase
from . import proc_base as p
from . import proc_bl
from . import proc_lp

pi = np.pi


###################
# Unit conversion #
###################


class unit_conversion(fileiobase.unit_conversion):
    """
    Unit converter class that returns NMRPipe like index values.  Useful
    when calling pipe_proc functions
    """
    # NMRPipe indexes from 1 to MAX instead on 0 to MAX-1
    # we need to modify two method to account for this off by one problem
    def __unit2pnt(self, val, units):
        return fileiobase.unit_conversion.__unit2pnt(self, val, units) + 1

    def __pnt2unit(self, val, units):
        return fileiobase.unit_conversion.__pnt2unit(self, val - 1, units)


def make_uc(dic, data, dim=-1):
    """
    Create a unit conversion object which accepts/returns NMRPipe indices.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters
    data : ndarray
        Array of NMR data.
    dim : int, optional
        Dimension number to create unit conversion object for.  Default is for
        last (direct) dimension.

    Returns
    -------
    uc : unit conversion object
        Unit conversion object for the given dimension which accepts and
        returns NMRPipe indices (starting from 1).

    """
    if dim == -1:
        dim = data.ndim - 1  # last dimention

    fn = "FDF" + str(int(dic["FDDIMORDER"][data.ndim - 1 - dim]))
    size = float(data.shape[dim])

    # check for quadrature in indirect dimentions
    if (dic[fn + "QUADFLAG"] != 1) and (dim != data.ndim - 1):
        size = size / 2.
        cplx = True
    else:
        cplx = False
    sw = dic[fn + "SW"]
    if sw == 0.0:
        sw = 1.0
    obs = dic[fn + "OBS"]
    if obs == 0.0:
        obs = 1.0
    car = dic[fn + "CAR"] * obs
    return unit_conversion(size, cplx, sw, obs, car)

########################
# Dictionary functions #
########################


def recalc_orig(dic, data, fn, axis=-1):
    """
    Recalculate the origin for given axis
    """
    # ORIG calculation
    s = float(data.shape[axis])

    # This really should check that the axis is not the last...
    if dic[fn + "QUADFLAG"] == 0 and axis != -1:
        s = int(s / 2.)

    # correct TPPI size in indirect dim when in time domain
    if dic["FD2DPHASE"] == 1 and fn != "FDF2" and dic[fn + "FTFLAG"] != 1:
        s = int(s / 2.)

    sw = dic[fn + "SW"]
    car = dic[fn + "CAR"]
    obs = dic[fn + "OBS"]
    s2 = float(dic[fn + "CENTER"])

    # DEBUG
    #print "Recalc of origin"
    #print "s:",s
    #print "axis:",axis
    #print "sw:",sw
    #print "car:",car
    #print "obs:",obs
    #print "s2:",s2

    dic[fn + "ORIG"] = car * obs - sw * ((s - s2) / s)
    return dic


def update_minmax(dic, data):
    """
    Update the MAX/MIN dictionary keys.
    """
    # maximum and minimum values
    dic["FDMAX"] = float(data.max().real)
    dic["FDDISPMAX"] = dic["FDMAX"]
    dic["FDMIN"] = float(data.min().real)
    dic["FDDISPMIN"] = dic["FDMIN"]
    dic["FDSCALEFLAG"] = 1.0    # FDMIN/MAX are valid
    return dic


def clean_minmax(dic):
    """
    Clean (set to zero) the MAX/MIN dictionary keys.
    """
    # maximum and minimum values
    dic["FDMAX"] = 0.0
    dic["FDDISPMAX"] = 0.0
    dic["FDMIN"] = 0.0
    dic["FDDISPMIN"] = 0.0
    dic["FDSCALEFLAG"] = 0.0    # FDMIN/MAX not valid
    return dic


#########################
# Apodization functions #
#########################


def apod(dic, data, qName=None, q1=1.0, q2=1.0, q3=1.0, c=1.0, start=1,
         size='default', inv=False, one=False, hdr=False):
    """
    Generic apodization.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    qName : {'SP', 'EM', 'GM', 'GMB', 'TM', 'TRI', 'JMOD'}
        Abbreviation of apodization function the apply. See the specific
        apodization function for a description.
    q1 : float
        First apodization function parameter. See specific apodization function
        for details.
    q2 : float
        Second apodization function parameter. See specific apodization
        function for details.
    q3 : float
        Third apodization function parameter. See specific apodization function
        for details.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with apodization applied.

    See Also
    --------
    em : Exponential apodization.
    gm : Lorentz-to-Gauss apodization.
    gmb : Modified Gaussian apodization.
    jmod : Exponentially damped J-modulation apodization.
    sp : Sine bell apodization.
    tm : Trapezoid apodization.
    tri : Triangular apodization.

    """
    if hdr:
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        qnum = dic[fn + "APODCODE"]
        qName = ["", "SP", "EM", "GM", "TM", "", "TRI", "GMB", "JMOD"][qnum]

    # Set apod codes here so that all three parameter are set
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    dic[fn + "APODQ1"] = q1
    dic[fn + "APODQ2"] = q2
    dic[fn + "APODQ3"] = q3

    if qName == "EM":
        return em(dic, data, q1, c, start, size, inv, one, hdr)
    elif qName == "GM":
        return gm(dic, data, q1, q2, q3, c, start, size, inv, one, hdr)
    elif qName == "GMB":
        return gmb(dic, data, q1, q2, c, start, size, inv, one, hdr)
    elif qName == "JMOD":
        return jmod(dic, data, q1, q2, q3, False, False, c, start, size, inv,
                    one, hdr)
    elif qName == "SP":
        return sp(dic, data, q1, q2, q3, c, start, size, inv, one, hdr)
    elif qName == "TM":
        return tm(dic, data, q1, q2, c, start, size, inv, one, hdr)
    elif qName == "TRI":
        return tri(dic, data, q1, q2, q3, c, start, size, inv, one, hdr)
    else:
        raise ValueError("qName must be SP, EM, GM, GMB, TM, TRI or JMOD")


def em(dic, data, lb=0.0, c=1.0, start=1, size='default', inv=False, one=False,
       hdr=False):
    """
    Exponential apodization.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    lb : float
        Exponential line broadening in Hz.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with exponential apodization applied.

    """
    start = start - 1   # arrays should start at 0

    # update dictionary
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data headers
        c = dic[fn + "C1"] + 1
        lb = dic[fn + "APODQ1"]
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 2.0
    dic[fn + "APODQ1"] = lb
    dic[fn + "APODQ2"] = 0.0
    dic[fn + "APODQ3"] = 0.0

    sw = dic[fn + "SW"]
    flb = lb / sw

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.em(data, lb=flb, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        data[..., start:stop] = p.em(data[..., start:stop], lb=flb, inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    dic = update_minmax(dic, data)
    return dic, data


def gm(dic, data, g1=0.0, g2=0.0, g3=0.0, c=1.0, start=1, size='default',
       inv=False, one=False, hdr=False):
    """
    Lorentz-to-Gauss apodization

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    g1 : float
        Inversion exponential width in Hz.
    g2 : float
        Gaussian broadening width in Hz
    g3 : float
        Location of Gaussian maximum, should be between 0.0 and 1.0.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with lorentz-to-gauss apodization applied.

    """
    start = start - 1  # arrays should start at 0
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data header
        g1 = dic[fn + "APODQ1"]
        g2 = dic[fn + "APODQ2"]
        g3 = dic[fn + "APODQ3"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 3.0
    dic[fn + "APODQ1"] = g1
    dic[fn + "APODQ2"] = g2
    dic[fn + "APODQ3"] = g3

    # calculate native parameters
    sw = dic[fn + "SW"]
    g1p = g1 / sw
    g2p = g2 / sw
    g3p = g3

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.gm(data, g1p, g2p, g3p, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size

        # pipe sets the maximum to the actual data maximum not
        # the maximum of the windowed region, so adj. g3p as necessary
        g3p = g3p * data.shape[-1] / (stop - start)
        #print start,stop,g1p,g2p,g3p
        data[..., start:stop] = p.gm(data[..., start:stop], g1p, g2p, g3p,
                                     inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    dic = update_minmax(dic, data)
    return dic, data


def gmb(dic, data, lb=0.0, gb=0.0, c=1.0, start=1, size='default', inv=False,
        one=False, hdr=False):
    """
    Modified Gaussian Apodization

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    lb : float
        Exponential apodization term in Hz.
    gb : float
        Gaussian apodization term in Hz.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a modified gaussian apodization applied.

    """
    start = start - 1  # arrays should start at 0
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data header
        lb = dic[fn + "APODQ1"]
        gb = dic[fn + "APODQ2"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 7.0
    dic[fn + "APODQ1"] = lb
    dic[fn + "APODQ2"] = gb
    dic[fn + "APODQ3"] = 0.0

    # calculate native parameters
    sw = dic[fn + "SW"]
    a = pi * lb / sw
    b = -a / (2.0 * gb * data.shape[-1])

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.gmb(data, a, b, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        data[..., start:stop] = p.gmb(data[..., start:stop], a, b, inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    dic = update_minmax(dic, data)
    return dic, data


def jmod(dic, data, off=0.0, j=0.0, lb=0.0, sin=False, cos=False, c=1.0,
         start=1, size='default', inv=False, one=False, hdr=False):
    """
    Exponentially Damped J-Modulation Apodization

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    off : float
        Starting location of J-modulation in a fractions of pi radians. This
        parameter is ignored if sin or cos parameters are True.
    j : float
        J-modulation in Hz.
    lb :
        Expoentntial line broadening in Hz.
    sin : bool
        True for sine modulation, off parameter is ignored.
    cos : bool
        True for cosine modulation, off parameter is ignored.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a exponentially damped J-modulation apodization
        applied.

    """
    start = start - 1  # arrays should start at 0
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if sin:
        off = 0.0
    if cos:
        off = 0.5

    if hdr:  # read apod values from data header
        off = dic[fn + "APODQ1"]
        j = dic[fn + "APODQ2"]
        lb = dic[fn + "APODQ3"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 8.0
    dic[fn + "APODQ1"] = off
    dic[fn + "APODQ2"] = j
    dic[fn + "APODQ3"] = lb

    # calculate native parameters
    sw = dic[fn + "SW"]
    e = pi * lb / sw
    end = off + j * (data.shape[-1] - 1) / sw

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.jmod(data, e, off, end, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        #print start,stop,e,off,end
        end = off + j * (stop - start - 1) / sw
        data[..., start:stop] = p.jmod(data[..., start:stop], e, off, end,
                                       inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c
    dic = update_minmax(dic, data)
    return dic, data


def sp(dic, data, off=0.0, end=1.0, pow=1.0, c=1.0, start=1, size='default',
       inv=False, one=False, hdr=False):
    """
    Sine bell apodization.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    off : float
        Starting location of sine-bell as a fraction of pi radians.
    end : float
        Ending location of sine-bell as a fraction of pi radians.
    pow : int
        Sine-bell power.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a sine-bell apodization applied.

    """
    start = start - 1  # arrays should start at 0
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data header
        off = dic[fn + "APODQ1"]
        end = dic[fn + "APODQ2"]
        pow = dic[fn + "APODQ3"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 1.0
    dic[fn + "APODQ1"] = off
    dic[fn + "APODQ2"] = end
    dic[fn + "APODQ3"] = pow

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.sp(data, off, end, pow, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        data[..., start:stop] = p.sp(data[..., start:stop], off, end, pow,
                                     inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    dic = update_minmax(dic, data)
    return dic, data


sine = sp   # wrapper for sine functions


def tm(dic, data, t1=0.0, t2=0.0, c=1.0, start=1, size='default', inv=False,
       one=False, hdr=False):
    """
    Trapezoid apodization.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    t1 : float
        Length in points of left side of the trapezoid.
    t2 : float
        Length in points of right side of the trapezoid.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a trapezoid apodization applied.

    """
    start = start - 1  # arrays should start at 0
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data header
        t1 = dic[fn + "APODQ1"]
        t2 = dic[fn + "APODQ2"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 4.0
    dic[fn + "APODQ1"] = t1
    dic[fn + "APODQ2"] = t2
    dic[fn + "APODQ3"] = 0.0

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.tm(data, t1=t1, t2=t2, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        data[..., start:stop] = p.tm(data[..., start:stop], t1=t1, t2=t2,
                                     inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    # check for NaN in array (when div by 0)
    if np.isnan(data).any():
        data = np.array(np.nan_to_num(data), dtype=data.dtype)

    dic = update_minmax(dic, data)
    return dic, data


def tri(dic, data, loc="auto", lHi=0.0, rHi=0.0, c=1.0, start=1,
        size='default', inv=False, one=False, hdr=False):
    """
    Triangular apodization

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    loc : int or "auto"
        Location in points of triangle apex. The default ("auto") is to place
        the apex in the middle.
    lHi : float
        Starting height of the left side of the triangle.
    rHi : float
        Starting height of the right side of the triangle.
    c : float
        First point scale value.
    start : int, optional
        Starting location of apodization window. Default is the first point, 1.
    size : int, optional
        Size of the apodization window. Default ('default') is the full size of
        the active dimension.
    inv : bool, optional
        True for inverse apodization, False for normal apodization.
    one : bool, optional
        True to set points outside of window to 1. False leaves points outside
        the apodization window as is.
    hdr : bool, optional
        True to read apodization parameters from the the parameters in dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a triangular apodization applied.

    Notes
    -----
    The right side of the apodization is differs slightly from NMRPipe's tri
    function.

    """
    start = start - 1  # arrays should start at 0

    if loc == "auto":
        loc = data.shape[-1] / 2

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if hdr:  # read apod values from data header
        loc = dic[fn + "APODQ1"]
        lHi = dic[fn + "APODQ2"]
        rHi = dic[fn + "APODQ3"]
        c = dic[fn + "C1"] + 1

    # update the dictionary
    dic[fn + "C1"] = c - 1.0

    # set the apod flags
    dic[fn + "APODCODE"] = 6.0
    dic[fn + "APODQ1"] = loc
    dic[fn + "APODQ2"] = lHi
    dic[fn + "APODQ3"] = rHi

    # apply apodization to data
    if start == 0 and size == 'default':
        data = p.tri(data, loc=loc, lHi=lHi, rHi=rHi, inv=inv)
    else:   # only part of the data window is apodized
        if size == 'default':
            stop = data.shape[-1]
        else:
            stop = start + size
        data[..., start:stop] = p.tri(data[..., start:stop], loc, lHi, rHi,
                                      inv=inv)
        if one is False:
            data[..., :start] = 0.0
            data[..., stop:] = 0.0

    # first point scaling
    if inv:
        data[..., 0] = data[..., 0] / c
    else:
        data[..., 0] = data[..., 0] * c

    dic = update_minmax(dic, data)
    return dic, data


###################
# Shift functions #
###################


def rs(dic, data, rs=0.0, sw=False):
    """
    Right shift and zero pad.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    rs : float
        Number of points to right shift. Negative values will left shift.
    sw : bool
        True to update chemical shift calibration parameters.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been right shifted.

    """
    if rs < 0:  # negative right shifts are left shifts
        return ls(dic, data, ls=-rs, sw=sw)

    data = p.rs(data, pts=rs)
    dic = update_minmax(dic, data)

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if sw and dic[fn + "FTFLAG"] == 1:
        # we are in freq domain and must update NDORIG and NDCENTER
        dic[fn + "CENTER"] = dic[fn + "CENTER"] + rs
        dic = recalc_orig(dic, data, fn)
    return dic, data


def ls(dic, data, ls=0.0, sw=False):
    """
    Left Shift and Zero Pad

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    ls : float
        Number of points to left shift. Negative values will right shift.
    sw : bool
        True to update chemical shift calibration parameters.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been left shifted.

    """
    if ls < 0:
        return rs(dic, data, rs=-ls, sw=sw)

    data = p.ls(data, ls)
    dic = update_minmax(dic, data)

    if sw:
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        if dic[fn + "FTFLAG"] == 1:   # freq domain
            # update NDORIG and NDCENTER
            dic[fn + "CENTER"] = dic[fn + "CENTER"] - ls
            dic = recalc_orig(dic, data, fn)
        else:   # time domain
            dic[fn + "APOD"] = data.shape[-1] - ls
            dic[fn + "TDSIZE"] = data.shape[-1] - ls

    return dic, data


def cs(dic, data, dir, pts=0.0, neg=False, sw=False):
    """
    Circular shift

    The syntax of this function is different from NMRPipe's CS function.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    dir : {'rs' or 'ls'}
        Direction to shift spectra, 'rs' for right shifting, 'ls' for left
        shifting.
    pts : float
        Number of points to shift.
    neg : bool
        True to negative points which are shifted.
    sw : bool
        True to update chemical shift calibration parameters.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been circular shifted.

    """
    if dir == "ls":
        pts = -pts
    elif dir != "rs":
        raise ValueError("dir must be ls or rs")

    data = p.cs(data, pts, neg=neg)
    dic = update_minmax(dic, data)

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if sw and dic[fn + "FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn + "CENTER"] = dic[fn + "CENTER"] + pts
        dic = recalc_orig(dic, data, fn)
    return dic, data


def fsh(dic, data, dir, pts, sw=True):
    """
    Frequency shift.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    dir : {'rs' or 'ls'}
        Direction to shift spectra, 'rs' for right shifting, 'ls' for left
        shifting.
    pts : float
        Number of points to shift.
    sw : bool
        True to update chemical shift calibration parameters.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been frequency shifted.

    Notes
    -----
    This function does not perfrom a Hilbert transfrom when data is complex,
    NMRPipe's FSH function appear to.  As such the results of the
    imaginary channel differs from NMRPipe. In addition MAX/MIN value are
    slightly different than those in NMRPipe.

    """
    if dir not in ["ls", "rs"]:
        raise ValueError("dir must be ls or rs")

    if np.iscomplexobj(data) is False:  # real data
        null, data = _ht(dict(dic), data, zf=True)
        del_imag = True
    else:   # imaginary data
        del_imag = False
        # NMRPipe always performs a hilbert transform
        # uncommenting the next two lines will match NMRPipe's fsh real
        # channel results, the imaginary channel is a mystery.
        #null,data = _ht(dict(dic),data,zf=True)
        #data = np.array(data,dtype="complex64")

    if dir == "ls":
        pts = -pts

    data = p.fsh(data, pts)

    dic = update_minmax(dic, data)

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    if dic[fn + "FTFLAG"] == 1 and sw:  # freq domain
        dic[fn + "CENTER"] = dic[fn + "CENTER"] + pts
        dic = recalc_orig(dic, data, fn)
    if del_imag is False:
        return dic, data
    else:
        return dic, data.real


##############
# Transforms #
##############


def ft(dic, data, auto=False, real=False, inv=False, alt=False, neg=False,
       null=False, bruk=False):
    """
    Complex Fourier transform.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    auto : bool
        True will choose mode automatically, not recomended.
    real : bool
        True to transform real-only data.
    inv : bool
        True to perform an inverse transform.
    alt : bool
        True to alternative the sign of points before transforming.
    neg : bool
        True will negate the imaginary channel before transforming.
    null : bool
        True will not apply transform but will update the parameter dictionary.
    bruk : bool
        True to process Redfield sequential data, this is the same as setting
        alt and real to True.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been Fourier transformed.

    Notes
    -----
    Choosing multiply conflicting modes can produces results different from
    NMRPipe's FT function.

    """
    size = data.shape[-1]

    # super-flags
    if auto:
        # turn off all flags
        real = False
        inv = False
        alt = False
        neg = False
        null = False
        bruk = False

        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        if dic[fn + "FTFLAG"] == 1.0:   # freq domain
            inv = True
        else:  # freq domain
            # Real, TPPI and Sequential data is real transform
            if dic["FDDIMCOUNT"] >= 2.:
                if ((dic["FD2DPHASE"] == 0 or dic["FD2DPHASE"] == 1) and
                        fn != "FDF2"):
                            real = True

            # sign and negation in AQSIGN
            if dic[fn + "AQSIGN"] == 1 or dic[fn + "AQSIGN"] == 2:
                alt = True

            if (dic[fn + "AQSIGN"] == 16 or dic[fn + "AQSIGN"] == 17 or
                    dic[fn + "AQSIGN"] == 18):
                        alt = True
                        neg = True

        # DEBUG
        #print "real:", real
        #print "inv:", inv
        #print "alt:", alt
        #print "neg:", neg
        #print "null:", null
        #print "bruk:", bruk

    if bruk:
        real = True
        alt = True

    if real:    # keep real data
        if np.iscomplexobj(data):
            data.imag = 0.0

    if alt:  # sign alternate
        if inv is False:    # inv with alt, alternates the inverse
            data[..., 1::2] = data[..., 1::2] * -1.

    if neg:  # negate the imaginary
        if np.iscomplexobj(data):
            data.imag = data.imag * -1.

    # update the dictionary
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    dic[fn + "AQSIGN"] = 0.0                       # we don't need sign alternation or negation anymore
    dic[fn + "FTFLAG"] = (dic[fn + "FTFLAG"] + 1) % 2   # toggle FT flag
    if dic[fn + "FTFLAG"] == 1:
        dic[fn + "FTSIZE"] = data.shape[-1]

    if null:
        # don't perform FT, just find min/max
        dic = update_minmax(dic, data)
        return dic, data

    # recast data if needed
    if data.dtype != "complex64":
        data = data.astype("complex64")

    if inv:  # inverse transform
        data = p.ifft_positive(data)
        if alt:
            data[..., 1::2] = data[..., 1::2] * -1
    else:
        data = p.fft_positive(data)

    if real:
        # return a complex array with double size
        data = np.array(data[..., size / 2:], dtype="complex64")

        # adjust quadrature
        dic[fn + "QUADFLAG"] = 0.0
        dic["FDQUADFLAG"] = 0.0

        # adjust size
        dic[fn + "APOD"] = dic[fn + "APOD"] / 2.0
        dic[fn + "TDSIZE"] = dic[fn + "TDSIZE"] / 2.0
        dic["FDSIZE"] = dic["FDSIZE"] / 2.0

    dic = update_minmax(dic, data)
    return dic, data


def rft(dic, data, inv=False):
    """
    Real Fourier transform.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    inv : bool
        True to perform an inverse transform.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been real Fourier transformed.

    Notes
    -----
    This function gives results which slightly differ from NMRPipe's RFT
    function in some cases.

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc

    # if dim is 2, direct dim is real and indirect is complex
    if (data.ndim == 2 and dic[fn + "QUADFLAG"] == 1 and
            dic[fn2 + "QUADFLAG"] == 0):
        data = data[::2]

    if inv:
        data = p.irft(data.real)
    else:
        data = p.rft(data.real)

    # update the dictionary
    dic[fn + "FTFLAG"] = (dic[fn + "FTFLAG"] + 1) % 2   # troggle FT flag
    dic[fn + "QUADFLAG"] = 1.0  # real data
    dic["FDQUADFLAG"] = 1.0    # real data
    dic = update_minmax(dic, data)
    return dic, data


def ha(dic, data, inv=False):
    """
    Hadamard transform.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    inv : bool
        True to perform an inverse transform.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been Hadamard transformed.

    Notes
    -----
    This function is slow.  Implemented a FWHT in proc_base would
    significantly improve the speed of this functions.

    """
    data = p.ha(data)

    if inv:
        data = data / data.shape[-1]

    dic = update_minmax(dic, data)
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    dic[fn + "FTFLAG"] = (dic[fn + "FTFLAG"] + 1) % 2   # troggle FT flag

    # calculation for dictionary updates
    s = data.shape[-1]
    s2 = s / 2.0 + 1
    dic[fn + "CENTER"] = s2
    dic = recalc_orig(dic, data, fn)
    dic["FDSIZE"] = s

    return dic, data


def ht(dic, data, mode="ps0-0", zf=False, td=False, auto=False):
    """
    Hilbert transform.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    mode : {'ps0-0', 'ps90-180'}
        Mirror image mode.
    zf : bool
        True to zero fill before transform for speed.
    td : bool
        True to set the time-domain parameter to half size.
    auto : bool
        True to select mode and zf parameters automatically from dic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been Hilbert transformed.

    Notes
    -----
    "ps90-180" mirror image mode gives different results than NMRPipe's HT
    function.

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    if auto:
        # when no cutting and 180 P1 correction
        if (dic[fn + "P1"] == 180.0 and dic[fn + "X1"] == 0.0 and
                dic[fn + "XN"] == 0.0):
            zf = False
            mode = "ps90-180"
        else:
            zf = True
            mode = "ps0-0"
    if mode not in ["ps0-0", "ps90-180"]:
        raise ValueError("mode must be ps0-0 or ps90-180")
    if mode == "ps90-180":
        # XXX determine how this works....
        pass
    if zf:
        N = 2 ** (np.ceil(np.log2(data.shape[-1])))  # not same as NMRPipe
    else:
        N = data.shape[-1]

    z = np.array(p.ht(data, N), dtype="complex64")
    dic = update_minmax(dic, data)

    # set the QUADFLAG as complex
    dic[fn + "QUADFLAG"] = 0.0
    if fn == "FDF2":
        dic["FDQUADFLAG"] = 0.0
    if td:
        dic[fn + "APOD"] = data.shape[-1] / 2.
    return dic, z


_ht = ht    # private function so ps can call the ht function


##########################
# Standard NMR Functions #
##########################


def di(dic, data):
    """
    Delete imaginaries

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data from which imaginaries have been removed.

    """
    data = p.di(data)
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if dic[fn + "QUADFLAG"] == 0.0:
        dic[fn + "QUADFLAG"] = 1

        if fn == "FDF2":
            # half number of point in indirect dim for States/
            if dic["FD2DPHASE"] != 0 and dic["FD2DPHASE"] != 1:
                dic["FDSPECNUM"] = dic["FDSPECNUM"] / 2.0
                dic["FDSLICECOUNT"] = dic["FDSPECNUM"]

        if dic["FDF1QUADFLAG"] == 1 and dic["FDF2QUADFLAG"]:
            dic["FDQUADFLAG"] = 1.0

    return dic, data


def ps(dic, data, p0=0.0, p1=0.0, inv=False, hdr=False, noup=False, ht=False,
       zf=False, exp=False, tc=0.0):
    """
    Phase shift

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    p0 : float
        Zero order phase in degrees.
    p1 : float
        First order phase in degrees.
    inv : bool
        True to perform inverse phase correction.
    hdr : bool
        True to use phasing parameters from dic.
    noup : bool
        True to not update phasing paramters in returned ndic.
    ht : bool
        True to perform a Hilbert transform to reconstruction imaginaries
        before phasing.
    zf : bool
        True to zero fill before applied Hilbert transform.
    exp : bool
        True to perform exponential phase correction.  False performs linear
        phase correction.
    tc : float, optional
        Exponential decay constant. User when exp is True.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been phased.

    Notes
    -----
    When inv is True this function will correctly invert an expoenential
    phase correction, NMRPipe's PS function does not. In addition, FDFNP0 and
    FDFNP1 are updated unless noup=True.  There are not rs and ls parameter, if
    the data need to be shifted before phasing use the :py:func:`rs` or
    :py:func`ls` function before using this function.

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if ht:  # Hilbert transform
        dic, data = _ht(dic, data, zf=zf)

    if hdr:  # read from header
        p0 = dic[fn + "P0"]
        p1 = dic[fn + "P1"]

    if exp:
        data = p.ps_exp(data, p0=p0, tc=tc, inv=inv)
    else:
        data = p.ps(data, p0=p0, p1=p1, inv=inv)

    if noup is False:
        dic[fn + "P0"] = p0
        dic[fn + "P1"] = p1

    dic = update_minmax(dic, data)

    return dic, data


def tp(dic, data, hyper=False, nohyper=False, auto=False, nohdr=False):
    """
    Transpose data (2D).

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    hyper : bool
        True to perfrom hypercomplex transpose.
    nohyper : bool
        True to supress hypercomplex transpose.
    auto : bool
        True to choose transpose mode automatically.
    nohdr : bool
        True to not update the transpose parameters in ndic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been transposed.

    """
    # XXX test if works with TPPI
    if nohyper:
        hyper = False

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc

    if auto:
        if (dic[fn + "QUADFLAG"] != 1) and (dic[fn2 + "QUADFLAG"] != 1):
            hyper = True
        else:
            hyper = False

    if hyper:   # Hypercomplex transpose need type recast
        data = np.array(p.tp_hyper(data), dtype="complex64")
    else:
        data = p.tp(data)
        if dic[fn2 + "QUADFLAG"] != 1 and nohyper is not True:
            # unpack complex as needed
            data = np.array(p.c2ri(data), dtype="complex64")

    # update the dimentionality and order
    dic["FDSLICECOUNT"] = data.shape[0]
    if data.dtype == 'float32':
        dic["FDSIZE"] = data.shape[1] / 2
    else:
        dic["FDSIZE"] = data.shape[1]


    dic["FDSPECNUM"] = dic["FDSLICECOUNT"]

    dic["FDDIMORDER1"], dic["FDDIMORDER2"] = (dic["FDDIMORDER2"],
                                              dic["FDDIMORDER1"])

    dic['FDDIMORDER'] = [dic["FDDIMORDER1"], dic["FDDIMORDER2"],
                         dic["FDDIMORDER3"], dic["FDDIMORDER4"]]

    if dic["FD2DPHASE"] == 0:
        dic['FDF1QUADFLAG'], dic['FDF2QUADFLAG'] = (dic['FDF2QUADFLAG'],
                                                    dic['FDF1QUADFLAG'])

    if nohdr is not True:
        dic["FDTRANSPOSED"] = (dic["FDTRANSPOSED"] + 1) % 2

    dic = clean_minmax(dic)
    return dic, data


ytp = tp    # alias for tp


xy2yx = tp  # alias for tp


def zf(dic, data, zf=1, pad="auto", size="auto", mid=False, inter=False,
       auto=False, inv=False):
    """
    Zero fill

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    zf : int, optional.
        Number of times to double the current dimensions size.
    pad : int
        Number of zeros to pad the data with.
    size : int
        Desired final size of the current dimension.
    mid : bool
        True to zero fill in the middle of the current dimension
    inter : bool
        True to zero fill between points.
    auto : bool
        True to round final size to nearest power of two.
    inv : bool
        True to extract the time domain data (remove zero filling).

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has zero filled.

    Notes
    -----
    Only one of the `zf`, `pad` and `size` parameter should be used, the other
    should be left as the default value.  If any of the `mid`, `inter`, `auto`
    and `inv` parameters are True other parameter may be ignored.

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if inv:  # recover original time domain points
        # calculation for dictionary updates
        s = dic[fn + "TDSIZE"]
        s2 = s / 2.0 + 1

        # update the dictionary
        dic[fn + "ZF"] = -1. * s
        dic[fn + "CENTER"] = s2
        dic = recalc_orig(dic, data, fn)
        dic["FDSIZE"] = s
        return dic, data[..., :s]

    if inter:   # zero filling between points done first
        data = p.zf_inter(data, zf)
        dic[fn + "SW"] = dic[fn + "SW"] * (zf + 1)
        zf = 0
        pad = 0  # NMRPipe ignores pad after a inter zf

    # set zpad, the number of zeros to be padded
    zpad = data.shape[-1] * 2 ** zf - data.shape[-1]

    if pad != "auto":
        zpad = pad
    if size != "auto":
        zpad = size - data.shape[-1]

    # auto is applied on top of other parameters:
    if auto:
        fsize = data.shape[-1] + zpad
        fsize = 2 ** (np.ceil(np.log(fsize) / np.log(2)))
        zpad = fsize - data.shape[-1]

    if zpad < 0:
        zpad = 0

    data = p.zf_pad(data, pad=zpad, mid=mid)

    # calculation for dictionary updates
    s = data.shape[-1]
    s2 = s / 2.0 + 1

    # update the dictionary
    dic[fn + "ZF"] = -1. * s
    dic[fn + "CENTER"] = s2
    if dic["FD2DPHASE"] == 1 and fn != "FDF2":   # TPPI data
        dic[fn + "CENTER"] = np.round(s2 / 2. + 0.001)
    dic = recalc_orig(dic, data, fn)
    dic["FDSIZE"] = s
    dic = update_minmax(dic, data)
    return dic, data


######################
# Baseline Functions #
######################


def base(dic, data, nl=None, nw=0, first=False, last=False):
    """
    Linear baseline correction.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    nl : list
        List of baseline node points.
    nw : int
        Node width in points.
    first : bool
        True to include the first point in the node list.
    last : bool
        True to include the last point in the node list.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a linear baseline correction applied.

    """
    if first:
        nl = [1] + nl
    if last:
        nl.append(data.shape[-1])

    # change values in node list to start at 0
    nl = [i - 1 for i in nl]

    data = proc_bl.base(data, nl, nw)
    dic = update_minmax(dic, data)
    return dic, data


def cbf(dic, data, last=10, reg=False, slice=slice(None)):
    """
    Constant baseline correction.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    last : float
        Percentage of trace to use in calculating baseline correction.
    reg : slice object, optional
        Python slice object describing region(s) from which to calculate the
        baseline correction. If False (default) the last parameter will be used
        to calculate the correction.
    slice : slice object
        Python slice describing regions to apply the baseline correction
        to.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a constant baseline correction applied.

    Notes
    -----
    The parameters of this function differ significantly from NMRPipe's cbf
    function.  The parameters `ref` and `slice` are Python slice objects
    if explicit correction regions are desired.  The `noseq` and `nodmx`
    parameters are not implemented.

    """
    if reg is not False:
        data = proc_bl.cbf_explicit(data, calc=reg, apply=slice)
    else:
        data = proc_bl.cbf(data, last, slice)
    dic = update_minmax(dic, data)
    return dic, data


def med(dic, data, nw=24, sf=16, sigma=5.0):
    """
    Median baseline correction

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    nw : int
        Median window size in points.
    sf : int
        Smoothing filter size in points.
    sigma : float
        Gaussian convolution width.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with an median baseline correction applied.

    Notes
    -----
    This function applies Friendrich's model-free baseline flatting algorithm
    (Friendrichs JBNMR 1995 5 147-153).  NMRPipe applies a different and
    unknown algorithm.

    """
    data = proc_bl.med(data, mw=nw, sf=sf, sigma=sigma)
    dic = update_minmax(dic, data)
    return dic, data


def sol(dic, data, mode="low", fl=16, fs=1, head=0):
    """
    Solvent filter

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    mode : {'low'}
        Filter mode.  Currenlty only 'low' is implemented.
    fl : int
        Length of filter in points.
    fs : {1, 2, 3}
        Shape of lowpass filter 1 : boxcar, 2: sine 3 : sine squared.
    head :
        Number of points to skip when applying filter.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a solvent filter applied.

    Notes
    -----
    This different from NMRPipe's SOL function in the only the low pass filter
    has been implemented. In addition the `mir`, `noseq` and `nodmx` parameters
    are not implemented.

    """
    if fs == 1:
        data[..., head:] = proc_bl.sol_boxcar(data[..., head:], w=fl * 2 + 1)
    elif fs == 2:
        data[..., head:] = proc_bl.sol_sine(data[..., head:], w=fl * 2 + 1)
    elif fs == 3:
        data[..., head:] = proc_bl.sol_sine2(data[..., head:], w=fl * 2 + 1)
    else:
        raise ValueError("fs must be 1, 2 or 3")
    dic = update_minmax(dic, data)
    return dic, data


###################
# Basic Utilities #
###################


def add(dic, data, r=0.0, i=0.0, c=0.0, ri=False, x1=1.0, xn='default'):
    """
    Add a constant

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    r : float
        Constant to add to real data.
    i : float
        Constant to add to imaginary data.
    c : float
        Constant to add to both real and imaginary data.
    ri : bool
        True to add real and imaginary data into real channel.
    x1 : int
        First point of region to add constant to.
    xn : int or 'default'
        Last point of region to add constant to.  'default' specifies the last
        point.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a constant added.

    Notes
    -----
    Parameter `c` is added to the real and imaginary data even when `r` and `i`
    are defined.  NMRPipe's ADD function ignores c when r or i are defined.

    """
    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if ri:
        data[..., mn:mx].real = p.add_ri(data[..., mn:mx])
    else:
        data[..., mn:mx] = p.add(data[..., mn:mx], r, i, c)

    dic = update_minmax(dic, data)
    return dic, data


def dx(dic, data):
    """
    Derivative by central difference.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Derivative of of NMR data.

    """
    data = p.dx(data)
    dic = update_minmax(dic, data)
    return dic, data


def ext(dic, data, x1="default", xn="default", y1="default", yn="default",
        round=1, left=False, right=False, mid=False, pow2=False, sw=True):
    """
    Extract a region.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    x1 : int or 'default'
        Starting point of the X-axis extraction. 'default' will start the
        extraction at the first point.
    xn : int or 'default'
        Ending point of the X-axis extraction. 'default' will stop the
        extraction at the last point.
    y1 : int or 'default'
        Starting point of the Y-axis extraction. 'default' will start the
        extraction at the first point.
    yn : int or 'default'
        Ending point of the Y-axis extraction. 'default' will stop the
        extraction at the last point.
    round : int
        Multiple to round extraction size to.
    left : bool
        True to extract the left half of the data.
    right : bool
        True to extract the right half of the data.
    mid : bool
        True to extract the central half of the data.
    pow2 : bool
        True will round the extracted size to the nearest power of 2.
    sw : bool
        True to update the sweep width and ppm calibration parameters,
        recommended.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Extracted region of NMR data.

    Notes
    -----
    The `time` parameter is not implemented.  Using multiple conflicting
    parameters may result in different results than NMRPipe.

    """
    # this function does not wrap proc_base.ext rather the slicing is
    # performed by this function

    # store old sizes
    old_x = float(data.shape[-1])
    if data.ndim == 2:
        old_y = float(data.shape[0])

    # slice find limits
    if x1 == "default":
        x_min = 0
    else:
        x_min = np.round(x1) - 1

    if xn == "default":
        x_max = data.shape[-1]
    else:
        x_max = np.round(xn)

    if left:
        x_min = 0
        x_max = int(data.shape[-1] / 2)

    if right:
        x_min = int(data.shape[-1] / 2)
        x_max = data.shape[-1]

    if mid:
        x_min = int(data.shape[-1] / 4)
        x_max = int(3 * data.shape[-1] / 4)

    r_x = round

    if pow2 and (x1 != "default" or xn != "default"):
        r_x = 2 ** np.ceil(np.log2(x_max - x_min))

    # round size to be multiple of r_x when axis is cut
    if x1 != "default" or xn != "default":
        remain_x = (x_min - x_max) % r_x     # -len_x%r_x
        x_min = x_min - np.floor(remain_x / 2)
        x_max = x_max + remain_x - np.floor(remain_x / 2)

    if x_min < 0:
        x_max = x_max - x_min
        x_min = 0.0

    if x_max > data.shape[-1]:
        x_min = x_min - (x_max - data.shape[-1])
        x_max = data.shape[-1]

    if data.ndim == 2:  # 2D array so we also have to modify y
        if y1 == "default":
            y_min = 0
        else:
            y_min = np.round(y1) - 1

        if yn == "default":
            y_max = data.shape[0]
        else:
            y_max = np.round(yn)

        r_y = round

        if pow2:
            r_y = 2 ** np.ceil(np.log2(y_max - y_min))

        # round only when axis is cut
        if y1 != "default" or yn != "default":
            remain_y = (y_min - y_max) % r_y
            y_min = y_min - np.floor(remain_y / 2)
            y_max = y_max + remain_y - np.floor(remain_y / 2)

        if y_min < 0:
            y_max = y_max - y_min
            y_min = 0.0

        if y_max > data.shape[0]:
            y_min = y_min - (y_max - data.shape[0])
            y_min = data.shape[0]

        #print "ymin:",y_min,"ymax:",y_max
        #print "xmin:",x_min,"xmax:",x_max

        data = data[y_min:y_max, x_min:x_max]
        if y_min != 1 and y_max != data.shape[0]:  # only update when sliced
            dic["FDSLICECOUNT"] = y_max - y_min
            dic["FDSPECNUM"] = y_max - y_min
        dic["FDSIZE"] = x_max - x_min

    else:       # 1D Array
        data = data[x_min:x_max]
        dic["FDSIZE"] = x_max - x_min

    # adjust sweep width and ppm calibration
    if sw:
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        s = data.shape[-1]

        if dic[fn + "FTFLAG"] == 0:   # time domain
            dic[fn + "CENTER"] = float(int(s / 2. + 1))
            dic[fn + "APOD"] = s
            dic[fn + "TDSIZE"] = s
            dic = recalc_orig(dic, data, fn)
        else:   # freq domain
            dic[fn + "X1"] = x_min + 1
            dic[fn + "XN"] = x_max
            dic[fn + "APOD"] = np.floor(dic[fn + "APOD"] * s / old_x)
            dic[fn + "CENTER"] = dic[fn + "CENTER"] - x_min
            dic[fn + "SW"] = dic[fn + "SW"] * s / old_x
            dic = recalc_orig(dic, data, fn)

        if data.ndim == 2:
            fn = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc
            s = data.shape[0]
            if dic[fn + "QUADFLAG"] == 0:
                s = s / 2

            if dic[fn + "FTFLAG"] == 0:  # time domain
                dic[fn + "CENTER"] = s / 2 + 1
                dic[fn + "APOD"] = s
                dic[fn + "TDSIZE"] = s
                dic = recalc_orig(dic, data, fn, -2)
            else:   # freq domain
                if y_min != 0:
                    dic[fn + "X1"] = y_min + 1
                if y_max != data.shape[0]:
                    dic[fn + "XN"] = y_max
                if y_min != 0 or y_max != data.shape[0]:
                    dic[fn + "APOD"] = np.floor(dic[fn + "APOD"] * s / old_y)
                dic[fn + "CENTER"] = dic[fn + "CENTER"] - y_min
                dic[fn + "SW"] = dic[fn + "SW"] * s / old_y
                dic = recalc_orig(dic, data, fn, -2)
    dic = update_minmax(dic, data)
    return dic, data


def integ(dic, data):
    """
    Integral by simple sum

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Integrated NMR data.

    """
    data = p.integ(data)
    dic = update_minmax(dic, data)
    return dic, data


def mc(dic, data, mode="mod"):
    """
    Modules or magnitude calculation.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    mode : {'mod' or 'pow'}
        'mod' to perform modules calculation, 'pow' to calculated square
        modules.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Modules of NMR data.

    """
    if mode == "mod":
        data = p.mc(data)
        dic["FDMCFLAG"] = 1.0
    elif mode == "pow":
        data = p.mc_pow(data)
        dic["FDMCFLAG"] = 2.0
    else:
        raise ValueError("mode must mod or pow")
    dic = update_minmax(dic, data)

    # change to mag. flags
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    dic[fn + "QUADFLAG"] = 1.0
    dic["FDQUADFLAG"] = 1.0
    return dic, data


def mir(dic, data, mode="left", invl=False, invr=False, sw=True):
    """
    Append mirror image.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    mode : {'left', 'right', 'center', 'ps90-180', pw0-0'}
        Type of mirror image to append.
    invl : bool
        True to negate left half.
    invr : bool
        True to negate right half.
    sw : bool
        True to update ppm parameters, recommended.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with mirror image appended.

    Notes
    -----
    Negations of selected region are applied regardless of the mode selected.

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if mode not in ['left', 'right', 'center', 'ps90-180', 'ps0-0']:
        raise ValueError("invalid mode")

    if dic[fn + "FTFLAG"] == 0:  # time domain
        if mode == "left":
            data = p.mir_left(data)
        if mode == "right":
            data = p.mir_right(data)
        if mode == "center":
            data = p.mir_center(data)
        if mode == "ps90-180":
            data = p.neg_edges(p.mir_center(data))
        if invr:
            data = p.neg_right(data)
        if invl:
            data = p.neg_left(data)
        dic["FDSIZE"] = dic["FDSIZE"] * 2
        if mode == "ps0-0":
            data = p.mir_center_onepoint(data)
            dic["FDSIZE"] = dic["FDSIZE"] - 1
    else:  # freq domain

        old_size = int(dic["FDSIZE"])

        if mode == "left":
            data = p.mir_left(data)
            dic[fn + "CENTER"] = old_size + dic[fn + "CENTER"]
        if mode == "right":
            data = p.mir_right(data)
            dic[fn + "CENTER"] = dic[fn + "CENTER"]
        if mode == "center":
            data = p.mir_center(data)
            dic[fn + "CENTER"] = dic[fn + "CENTER"] + old_size / 2.
        if mode == "ps90-180":
            data = p.neg_edges(p.mir_center(data))
            dic[fn + "CENTER"] = dic[fn + "CENTER"] + old_size / 2.
        if mode == "ps0-0":
            data = p.mir_center_onepoint(data)
            dic[fn + "CENTER"] = dic[fn + "CENTER"] + old_size
        if invr:
            data = p.neg_right(data)
        if invl:
            data = p.neg_left(data)

        # dictionary updates
        dic["FDSIZE"] = data.shape[-1]
        dic[fn + "APOD"] = dic["FDSIZE"]
        dic[fn + "FTSIZE"] = dic["FDSIZE"]
        dic[fn + "TDSIZE"] = dic["FDSIZE"]
        dic[fn + "ZF"] = -dic["FDSIZE"]
        s = dic["FDSIZE"]
        dic[fn + "SW"] = dic[fn + "SW"] * float(s) / float(old_size)
        dic = recalc_orig(dic, data, fn)  # recalculate origin

    dic = update_minmax(dic, data)
    return dic, data


def mult(dic, data, r=1.0, i=1.0, c=1.0, inv=False, hdr=False, x1=1.0,
         xn='default'):
    """
    Multiple by a constant.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    r : float
        Constant to multiply real data by.
    i : float
        Constant to multiply imaginary data by.
    c : float
        Constant to multiply both real and imaginary data by.
    inv : bool
        True to multiply by the inverse of the constant.
    hdr : bool
        True to use constant defined in dic.
    x1 : int
        First point of region to multiply constant by.
    xn : int or 'default'
        Last point of region to multiple constant by.  'default' specifies the
        last point.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been multiplied by a constant.

    Notes
    -----
    Parameter `c` is used even when `r` and `i` are defined.  NMRPipe's MULT
    function ignores c when r or i are defined.

    """
    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if hdr:  # read in C from header
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        c = dic[fn + "C1"]
        r = 1.0
        i = 1.0

    rf = (r * c)  # real factor
    cf = (i * c)  # complex factor
    if inv:
        rf = 1 / rf
        cf = 1 / cf

    data[..., mn:mx] = p.mult(data[..., mn:mx], r=rf, i=cf, c=1.0)
    dic = update_minmax(dic, data)
    return dic, data


def rev(dic, data, sw=True):
    """
    Reverse data.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    sw : bool
        True to update carrier parameters in ndic, recommended.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been revesed.

    """
    data = p.rev(data)
    dic = update_minmax(dic, data)
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    if sw and dic[fn + "FTFLAG"] == 1:
        # freq domain update NDORIG and NDCENTER
        dic[fn + "CENTER"] = dic[fn + "CENTER"] - 1
        dic = recalc_orig(dic, data, fn)
    return dic, data


def set(dic, data, r="a", i="a", c="a", x1=1.0, xn='default'):
    """
    Set data to a constant.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    r : float, or 'a'
        Constant to set real data to. "a" sets real data to 0 unless `c` is
        defined.
    i : float or 'a'
        Constant to set imaginary data to. "a" sets imaginary data to 0 unless
        `c` is defined.
    c : float
        Constant to set both real and imaginary data by. 'a' sets both channels
        to 0 unless r or i in defined.
    x1 : int
        First point of region to set to the constant.
    xn : int or 'default'
        Last point of region to set to the constant.  'default' specifies the
        last point.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been set to a constant.

    """
    mn = x1 - 1
    if xn == 'default':
        mx = data.shape[-1]
    else:
        mx = xn

    if r == "a" and i == "a" and c == "a":
        rc = 0
        ic = 0

    if c != "a":
        rc = c
        ic = c

    if r != "a":
        rc = r

    if i != "a":
        ic = i

    # this is so simple we do not use the proc_base functions
    data[..., mn:mx].real = rc
    if np.iscomplex(data).any():
        data[..., mn:mx].imag = ic

    dic = update_minmax(dic, data)
    return dic, data


def shuf(dic, data, mode=None):
    """
    Shuffle Utilities

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    modes : str
        Shuffle mode. See table below for valid modes.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data after shuffle.


    Notes
    -----
    'rr2ri' mode ignores any imaginary vector refusing to create a mis-sized
    vector.  'bswap' mode may results in NaN in the data.  'r2i' and 'i2r' not
    implemented.  All modes correctly update minimum and maximum values.  This
    behavor may differ from NMRPipe's SHUF function.

    Valid modes are:

    ======= ===================================
    string  Description
    ======= ===================================
    'ri2c'  Interleave real and imaginary data.
    'c2ri'  Seperate real and imaginary data.
    'ri2rr' Append real and imaginary data.
    'rr2ri' Unappend real and imaginary data.
    'exlr'  Exchange left and right halfs.
    'rolr'  Rotate left and right halfs.
    'swap'  Swap real and imaginary data.
    'bswap' Byte-swap data.
    'inv'   Do nothing.
    ======= ===================================

    """
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

    if mode == "ri2c":
        data = p.ri2c(data)  # interleave real and imaginary data
        # update the dictionary
        dic["FDQUADFLAG"] = 1.0
        dic[fn + "QUADFLAG"] = 1.0
        dic[fn + "APOD"] = data.shape[-1]
        dic[fn + "TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]
    elif mode == "c2ri":
        # seperate real and imaginary
        data = np.array(p.c2ri(data), dtype="complex64")
        # update the dictionary
        dic["FDQUADFLAG"] = 0.0
        dic[fn + "QUADFLAG"] = 0.0
        dic[fn + "APOD"] = data.shape[-1]
        dic[fn + "TDSIZE"] = data.shape[-1]
        dic["FDSIZE"] = data.shape[-1]
        dic["FDREALSIZE"] = data.shape[-1]
    elif mode == "ri2rr":
        data = p.ri2rr(data)    # appended imaginary data
        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0] / 2.0
            dic["FDSPECNUM"] = data.shape[0] / 2.0
        dic["FDQUADFLAG"] = 0.0
        dic[fn + "QUADFLAG"] = 1.0
        dic["FDSIZE"] = data.shape[-1]
    elif mode == "rr2ri":
        # unappend imaginary data (ignores imag data)
        data = np.array(p.rr2ri(data), dtype="complex64")
        # update the dictionary
        if data.ndim == 2:
            dic["FDSLICECOUNT"] = data.shape[0]
            dic["FDSPECNUM"] = data.shape[0]
        dic["FDQUADFLAG"] = 0.0
        dic[fn + "QUADFLAG"] = 0.0
        dic["FDSIZE"] = data.shape[-1]
    elif mode == "exlr":
        data = p.exlr(data)  # exchange left and right
    elif mode == "rolr":
        data = p.rolr(data)  # rotate left right halves
    elif mode == "swap":
        data = p.swap(data)
    elif mode == "bswap":
        data = p.bswap(data)
    elif mode == "r2i":
        raise NotImplementedError("Integer mode not implemented")
    elif mode == "i2r":
        raise NotImplementedError("Integer mode not implemented")
    elif mode == "inv":
        # This does not seem to do anything....
        #XXX check data with odd number of points
        pass
    else:
        raise ValueError("Invalid mode")
    # update the dictionary
    dic = update_minmax(dic, data)
    return dic, data


def sign(dic, data, ri=False, r=False, i=False, left=False, right=False,
         alt=False, abs=False, sign=False):
    """
    Sign manipulation utilities

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    ri : bool
        True to negate all data points.
    r : bool
        True to negate real data points.
    i : bool
        True to negate imaginary data points.
    left : bool
        True to negate the left half of the data.
    right : bool
        True to negate the right half of the data.
    alt : bool
        True to negate alternating data points.
    abs : bool
        True to replace both real and imaginary data with it's absolute value.
    sign : bool
        True to replace data with the sign (-1 or 1) of the data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data after sign manipulation.

    Notes
    -----
    All sign manupulation modes set True are applied in the order they appear
    in the function parameter list.

    """
    if ri:
        data = p.neg_all(data)
    if r:
        data = p.neg_real(data)
    if i:
        data = p.neg_imag(data)
    if left:
        data = p.neg_left(data)
    if right:
        data = p.neg_right(data)
    if alt:
        data = p.neg_alt(data)
    if abs:
        data = p.abs(data)
    if sign:
        data = p.sign(data)
    dic = update_minmax(dic, data)
    return dic, data


##################
# Misc Functions #
##################


def coadd(dic, data, cList=[1, 1], axis='x', time=False):
    """
    Co-addition of data

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    cList : list
        List of co-addition coefficients.
    axis : {'x', 'y'}
        Axis to co-add to and from.
    time : bool
        True will adjust time-domain parameters in dic to account for the size
        reduction.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data after co-addition.

    """
    if axis == 'x':
        data = p.coadd(data, cList, axis=-1)
        dic["FDSIZE"] = data.shape[-1]
        idx = 0
    elif axis == 'y':
        data = p.coadd(data, cList, axis=0)
        dic["FDSLICECOUNT"] = dic["FDSPECNUM"] = data.shape[0]
        idx = 1
    else:
        raise ValueError("axis must be x or y")

    dic = update_minmax(dic, data)
    if time:
        fn = "FDF" + str(int(dic["FDDIMORDER"][idx]))  # F1, F2, etc
        dic[fn + "APOD"] = np.floor(dic[fn + "APOD"] / len(cList))
        dic[fn + "TDSIZE"] = np.floor(dic[fn + "TDSIZE"] / len(cList))
    return dic, data


coad = coadd    # macro for coadd


def dev(dic, data):
    """
    Development function (does nothing)

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndic : dict
        Original dictionary of NMRPipe parameters.
    ndata : ndarray
        Original array of NMR data.

    """
    return dic, data


def img(dic, data, filter, dx=1.0, dy=1.0, kern=[1], conv=False, thres=None):
    """
    Image processing utilities

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    filter : {'median', 'min', 'max', 'amin', 'amax', 'range', 'avg', 'dev'}
        Image processing filter to apply. See table below for descriptions.
    dx : float
        Filter width along X-axis in points.
    dy : float
        Filter width along Y-axis in points.
    kern : list
        List of convolution filter kernel values, only used when conv is True.
    conv : bool
        True to apply convolution filter with kernel of kern.
    thres : float or None or True
        Threshold value.  Only points above this value will have the filter
        applied. None turns thresholding off and the filter will be applied
        to all points. True will set a threshold value of 0.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Imaged processed array of NMR data.

    Notes
    -----
    This function wraps when regions extend past the edges (NMRPipe doesn't).
    The filter is applied to both the real and imaginary channels

    Supported filters are:

    ======  ==================
    Name    Description
    ======  ==================
    median  Median
    min     Minimum
    max     Maximim
    amin    Absolute Minimum
    amax    Absolute Maximum
    range   Range
    avg     Average
    dev     Standard Deviation
    ======  ==================

    """
    # deal with thres by making a masked array
    if thres is not None:
        if thres is True:
            thres = 0.0  # default value of 0.0
        data = p.thres(data, thres)

    if conv:    # convolution with kernal
        data = p.conv(data, kern, m="wrap")
        dic = update_minmax(dic, data)
        return dic, data

    s = (2 * dy + 1, 2 * dx + 1)  # size tuple
    # the various filters
    if filter == "median":
        data = p.filter_median(data, s=s, m="wrap")
    elif filter == "min":
        data = p.filter_min(data, s=s, m="wrap")
    elif filter == "max":
        data = p.filter_max(data, s=s, m="wrap")
    elif filter == "amin":
        data = p.filter_amin(data, s=s, m="wrap")
    elif filter == "amax":
        data = p.filter_amax(data, s=s, m="wrap")
    elif filter == "range":
        data = p.filter_range(data, s=s, m="wrap")
    elif filter == "avg":
        data = p.filter_avg(data, s=s, m="wrap")
    elif filter == "dev":
        data = p.filter_dev(data, s=s, m="wrap")
    else:
        raise ValueError("Invalid filter")

    dic = update_minmax(dic, data)
    return dic, data


def null(dic, data):
    """
    Null function

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Unmodified array of NMR data.

    """
    dic = update_minmax(dic, data)  # Null actually does update this...
    return dic, data


def qart(dic, data, a=0.0, f=0.0, auto=False):
    """
    Scale Quad Artifacts

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    a : float
        Amplitude adjustment value.
    f : float
        Phase adjustment value.
    auto : bool
        True will perform a Gram-Schmidth orthorginalization to fund `a` and
        `f` automatically.  Provided `a` and `f` parameters are ignored.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with scaled quadrature artifacts.

    Notes
    -----
    Auto mode performs Gram-Schmidt orthogonalization, a different approach
    using a grid search is used in NMRPipe's QART function.

    """
    if auto:
        data = p.qart_auto(data)
    else:
        data = p.qart(data, a, f)

    dic = update_minmax(dic, data)
    return dic, data


def qmix(dic, data, ic=1, oc=1, cList=[0], time=False):
    """
    Complex mixing of input to outputs

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    ic : int
        Number of input channels
    oc : int
        Number of output channels
    cList : array_like
        Array or mixing coefficients.  This parameter must be able to be
        converted to an array and reshaped to (ic, oc).

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data after complex mixing.

    Notes
    -----
    `ic` and `oc` must evenly divide the data. This function refuses to make
    invalid length files, NMRPipe's qmix function will create such files.

    """
    ic = int(ic)
    oc = int(oc)

    if data.ndim != 2:
        raise ValueError("data must be 2D")

    if data.shape[0] % ic != 0 or data.shape[0] % oc != 0:
        raise ValueError("ic and oc must be divide the number of vectors")

    carr = np.array(cList, dtype='float').reshape(ic, oc)
    data = p.qmix(data, carr)

    #data = n
    dic = update_minmax(dic, data)
    dic["FDSPECNUM"] = data.shape[0]
    dic["FDSLICECOUNT"] = data.shape[0]

    if time:
        fn = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc
        dic[fn + "APOD"] = data.shape[0]
        dic[fn + "TDSIZE"] = data.shape[0]
        if dic[fn + "QUADFLAG"] == 0:
            dic[fn + "APOD"] = dic[fn + "APOD"] / 2.
            dic[fn + "TDSIZE"] = dic[fn + "TDSIZE"] / 2.
    return dic, data


def save(dic, data, name, overwrite=True):
    """
    Save the current vector.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    name : str
        Filename to save vector to.
    overwrite : bool
        True will overwrite existing files, False will raise a Warning if the
        file already exists.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Unmodified array of NMR data.

    Notes
    -----
    The resulting FDPIPECOUNT header parameter does not match the one created
    using NMRPipe's SAVE function.

    """
    dic["FDPIPECOUNT"] = 1.0

    pipe.write_single(name, dic, data, overwrite)

    dic["FDPIPECOUNT"] = 0.0
    dic = update_minmax(dic, data)
    return dic, data


def smo(dic, data, n=1, center=False):
    """
    Smooth data.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    n : int
        Size of smoothing window in points.
    center : bool
        True will perform perform a centering on the data (subtract the
        smoothed data).  False returns the smoothed data.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been smoothed or centered.

    """
    a = p.smo(data, n=n)
    # NMRPipe doesn't truely smooth the left edge of the vector
    for i in range(n):
        a[..., i] = data[..., 0:(n + i)].sum(axis=-1) / (n + 1 + i)
    if center:
        # to avoid the same error center without use proc_base functions
        a = data - a
    dic = update_minmax(dic, a)
    return dic, a


def zd(dic, data, wide=1.0, x0=1.0, slope=0, func=0, g=1):
    """
    Zero diagonal band.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    wide : int
        Width of the diagonal band in points.
    x0 : int
        Starting location of the diagonal band in points.
    slope : float
        Slope of the diagonal band (X/Y ratio). A value of 0 will determine the
        slope automatically.
    func : {0, 1, 2, 3}
        Function to perform zero-ing with. 0 for a boxcar window, 1 for a
        triangle window, 2 for a sine bell, 3 for a Gaussian.
    g : float
        Width of Gaussian window in points. Only used if func is 3.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a diagonal band zero-ed.

    """
    if x0 == 0:      # pipe takes x0=0 to be x0=1
        x0 = 1.0

    if slope == 0:    # Auto Mode
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc
        sw1 = dic[fn + "SW"]
        sw2 = dic[fn2 + "SW"]
        slope = data.shape[-1] * sw1 / (data.shape[0] * sw2)

    if func == 0:
        data = p.zd_boxcar(data, wide, x0 - 1, slope)
    elif func == 1:
        data = p.zd_triangle(data, wide, x0 - 1, slope)
    elif func == 2:
        data = p.zd_sinebell(data, wide, x0 - 1, slope)
    elif func == 3:
        data = p.zd_gaussian(data, wide, x0 - 1, slope, g)
    else:
        raise ValueError("func parameter must be 0, 1, 2 or 3")
    dic = update_minmax(dic, data)
    return dic, data


###############################
# Linear Prediction Functions #
###############################


def lp(dic, data, pred="default", x1="default", xn="default", ord=8, mode='f',
       append='after', bad_roots='auto', mirror=None, fix_mode='on',
       method='tls'):
    """
    Linear Prediction

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    pred : int
        Number of points to predict, "default" chooses the vector size for
        forward prediction, 1 for backward prediction
    x1 : int or 'default'
        First point in 1D vector to use to extract LP filter. 'default' will
        use the first or last point depending on the mode.
    xn : int or 'default'
        Last point in 1D vector to use to extract LP filter. 'default' will use
        the first or last point depending on the mode.
    ord : int
        Prediction order, number of LP coefficients used in prediction.
    mode : {'f', 'b', 'fb'}
        Mode to generate LP filter, 'f' for forward,'b' for backward, 'fb' for
        forward-backward.
    append : {'before' or 'after'}
        Location to append predicted data, 'before' or 'after' the existing
        data.
    bad_roots {'incr', 'decr', None, 'auto'} :
        Type of roots which are will be marked as bad and stabilized. Choices
        are 'incr' for increasing roots, 'decr' for decreasing roots, or None
        for not root stabilization. The default 'auto' will set this parameter
        based upon the LP `mode` parameter: 'f' and 'fb' will results in an
        'incr' parameter. 'b' in 'decr'.
    mirror : {'90-180', '0-0', None}
        Mirror mode, option are '90-180' for a one point shifted mirror image,
        '0-0' for an exact mirror image, and None for no mirror imaging of the
        data.
    fix_mode : {'on', 'reflect'}
        Method used to stabilize bad roots, 'on' moves bad roots onto the unit
        circle, 'reflect' reflect bad roots across the unit circle.
    method : {'svd', 'qr', 'choleskey', 'tls'}
        Method to use to calculate the LP filter.

    Notes
    -----
    The results from this function do not match NMRPipe's LP function.  Also
    some additional parameter and different parameter in this function.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with linear prediction applied.

    """
    # check parameter
    if mirror not in [None, '90-180', '0-0']:
        raise ValueError("mirror must be None, '90-180' or '0-0'")

    # pred default values
    if pred == "default":
        if mode == "after":
            pred = data.shape[-1]   # double the number of points
        else:
            pred = 1    # predict 1 point before the data

    # remove first pred points if appending before data
    if append == "before":
        data = data[..., pred:]

    # create slice object
    if x1 == "default":
        x_min = 0
    elif mode == "before":
        x_min = x1 - pred - 1
    else:
        x_min = x1 - 1

    if xn == "default":
        x_max = data.shape[-1]
    else:
        x_max = xn - 1
    sl = slice(x_min, x_max)

    # mirror mode (remap to proc_lp names
    mirror = {None: None, '90-180': '180', '0-0': '0'}[mirror]

    # mode, append, bad_roots, fix_mode, and method are passed unchanged
    # use LP-TLS for best results
    data = proc_lp.lp(data, pred, sl, ord, mode, append, bad_roots, fix_mode,
                      mirror, method)

    # calculation for dictionary updates
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    s = data.shape[-1]
    s2 = s / 2.0 + 1

    # update the dictionary
    dic[fn + "CENTER"] = s2
    if dic["FD2DPHASE"] == 1 and fn != "FDF2":   # TPPI data
        dic[fn + "CENTER"] = np.round(s2 / 2. + 0.001)
    dic = recalc_orig(dic, data, fn)
    dic["FDSIZE"] = s
    dic[fn + "APOD"] = s
    dic[fn + "TDSIZE"] = s

    dic = update_minmax(dic, data)
    return dic, data


lpc = lp        # macro to lp


def lp2d(dic, data, xOrd=8, yOrd=8, xSize="default", ySize="default",
         xMirror='0', yMirror='0', fix_pts=True, method='svd'):
    """
    2D Linear Prediction using LP2D procedure

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    xOrd : int
        X dimension linear prediction order.
    yOrd : int
        Y dimension linear prediction order.
    xSize : int
        New size of Y-axis, 'default' doubles the current size.
    ySize : int
        New size of Y-axis, 'default' double the current size.
    xMirror : {'0', '180'}
    '   Mode in which the mirror image of the X-axis should be formed.  '0'
        indicated no delay, '180' for a half-point delay.
    yMirror : {'0', '180'}
        Mode in which the mirror image of the Y-axis should be formed.
    fix_pts : bool
        True to reduce predicted points with magnitude larger than the largest
        data point. False leaved predicted points unaltered.
    method : {'svd', 'qr', 'cholesky', 'tls'}
        Method used to calculate the LP prediction filter.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with 2D linear prediction applied.

    Notes
    -----
    This function applies the LP2D procedure as described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 98, 192-199.
    to the data matrix. The parameters and algorith used in NMRPipe's LP2D
    function are not well documented and are not replicated here.

    """
    # determind how many points to predict in each dimension
    if xSize == "default":
        xpred = data.shape[1]
    else:
        xpred = xSize - data.shape[1]
    if ySize == "default":
        ypred = data.shape[0]
    else:
        ypred = ySize - data.shape[0]

    # predict the X (last) axis.
    data = proc_lp.lp2d(data, xpred, yOrd, xOrd, yMirror, fix_pts, method)

    # transpose the data matrix, predict Y axis, tranpose back
    data = data.T
    data = proc_lp.lp2d(data, ypred, xOrd, yOrd, xMirror, fix_pts, method)
    data = data.T

    # update dictionary
    # x-axis updates
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    s = data.shape[1]
    s2 = s / 2.0 + 1

    # update the dictionary
    dic[fn + "CENTER"] = s2
    if dic["FD2DPHASE"] == 1 and fn != "FDF2":   # TPPI data
        dic[fn + "CENTER"] = np.round(s2 / 2. + 0.001)
    dic = recalc_orig(dic, data, fn)
    dic["FDSIZE"] = s
    dic[fn + "APOD"] = s
    dic[fn + "TDSIZE"] = s

    # y-axis updates
    fn = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc
    s = data.shape[0]
    s2 = s / 2.0 + 1

    # update the dictionary
    dic[fn + "CENTER"] = s2
    if dic["FD2DPHASE"] == 1 and fn != "FDF2":   # TPPI data
        dic[fn + "CENTER"] = np.round(s2 / 2. + 0.001)
    dic = recalc_orig(dic, data, fn)
    dic[fn + "APOD"] = s
    dic[fn + "TDSIZE"] = s

    dic = update_minmax(dic, data)
    return dic, data

#############################
# Not Implemented Functions #
#############################


def ann(dic, data):
    """
    Fourier Analysis by Neural Net
    """
    raise NotImplementedError


def ebs(dic, data):
    """
    EBS Reconstruction
    """
    raise NotImplementedError


def mac(dic, data):
    """
    Macro Language Interpreter
    """
    raise NotImplementedError


def mem(dic, data):
    """
    Maximum Entropy Reconstruction
    """
    raise NotImplementedError


def ml(dic, data):
    """
    Maximum Likelihood Frequency Map
    """
    raise NotImplementedError


def poly(dic, data):
    """
    Polynomial Baseline Correction
    """
    raise NotImplementedError


def xyz2zyx(dic, data):
    """
    3D Matrix transpose
    """
    raise NotImplementedError


def ztp(dic, data):
    """
    3D Matrix Transpose
    """
    raise NotImplementedError
