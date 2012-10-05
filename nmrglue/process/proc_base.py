"""
A collection of NMR spectral processing functions which operate on the last
dimension (1) of 2D arrays.  These functions are wrapped by other processing
modules but can also be used directly.  All parameter are assumed to be in
units of points unless otherwise noted.
"""

# TODO determine which of these work on N-dimension and which assume 2D

import numpy as np
import scipy.signal
import scipy.linalg

pi = np.pi


#########################
# Apodization functions #
#########################


def em(data, lb=0.0, inv=False):
    """
    Exponential apodization

    Functional form of apodization window:

    .. math::
        em(x_i) = \\exp(-pi * i * lb)

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    lp : float
        Exponential line broadening, in units of points.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    apod = np.exp(-pi * np.arange(data.shape[-1]) * lb).astype(data.dtype)
    if inv:
        apod = 1 / apod   # invert apodization
    return apod * data


def gm(data, g1=0.0, g2=0.0, g3=0.0, inv=False):
    """
    Lorentz-to-Gauss apodization

    Functional form of apodization window:

    .. math::
        gm(x_i) = \\exp(e - g^2)

    Where:

    .. math::
        e = pi * i * g1 \\\\
        g = 0.6 * pi * g2 * (g3 * (size - 1) - i)

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    g1 : float
        Inverse exponential width.
    g2 : float
        Gaussian broadening width.
    g3 : float
        Location of Gaussian maximum.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    e = pi * np.arange(size) * g1
    g = 0.6 * pi * g2 * (g3 * (size - 1) - np.arange(size))
    apod = np.exp(e - g * g).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


def gmb(data, a=0.0, b=0.0, inv=False):
    """
    Modified gaussian apodization

    Functional form of apodization window:

    .. math::
        gmb(x_i) = \\exp(-a * i - b * i^2)

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    a : float
        Exponential term in apodization.
    b : float
        Gaussian term in apodization.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    apod = np.exp(-a * np.arange(size) - b *
                  np.arange(size) ** 2).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


def jmod(data, e=0.0, off=0.0, end=0.0, inv=False):
    """
    Exponentially damped J-modulation apodization

    Functional form of apodization window:

    .. math::
        jmod(x_i) = \\exp(-e) * \\sin(\\frac{pi * off + pi * (end - off) * i}
                                            {size - 1})

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    e : float
        Exponential apodization term.
    off : float
        Start of J-modulation in fractions of pi radians (180 degrees).
    end : float
        End of J-modulation in fractions of pi radians (180 degrees).
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    apod = (np.exp(-e * np.arange(size)).astype(data.dtype) *
            np.sin(pi * off + pi * (end - off) * np.arange(size) /
                (size - 1)).astype(data.dtype))
    if inv:
        apod = 1 / apod
    return apod * data


def sp(data, off=0, end=1.0, pow=1.0, inv=False):
    """
    Shifted sine-bell apodization

    Functional form of apodization window:

    .. math::
        sp(x_i) = \\sin(\\frac{pi * off + pi * (end - off) * i}
                              {size - 1})^{pow}

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    off : float
        Offset start of sine-bell as a percentage of the trace (0 -> 1).
    end : float
        End of of sine-bell as a percentage of the trace (0 -> 1).
    pow : float
        Power to raise sine-bell to.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    apod = np.power(np.sin(pi * off + pi * (end - off) * np.arange(size) /
                    (size - 1)).astype(data.dtype), pow).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


sine = sp


def tm(data, t1=0.0, t2=0.0, inv=False):
    """
    Trapezoid Apodization

    Functional form of apodization:

    ============ ================================
    Range        Description
    ============ ================================
    0:t1         Linear increases from 0.0 to 1.0
    t1:size - t2 Flat with value of 1.0
    -t2:         Linear descrease from 1.0 to 0.0
    ============ ================================

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    t1 : int
        Length of left ramp in points.
    t2 : int
        Length of right ramp in points.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    apod = np.concatenate((np.linspace(0, 1, t1), np.ones(size - t1 - t2),
                           np.linspace(1, 0, t2))).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


def tri(data, loc="auto", lHi=0.0, rHi=0.0, inv=False):
    """
    Triangle apodization.

    Functional form of apodization window.

    ============ ========================================
    Range        Description
    ============ ========================================
    0:loc        Linear increase/decrease from lHi to 1.0
    loc:         Linear increase/decrease from 1.0 to rHi
    ============ ========================================

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    loc : int or "auto"
        Location of triangle apex. 'auto' set the apex to the middle of the
        trace (default).
    lHi : float
        Starting height of the left side of the triangle.
    rHi : float
        Starting height of the right side of the triangle.
    inv : bool, optional
        True for inverse apodization.  False (default) for standard.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with apodization applied.

    """
    size = data.shape[-1]
    if loc == "auto":
        loc = int(size / 2.)
    apod = np.concatenate((np.linspace(lHi, 1., loc), np.linspace(1., rHi,
                            size - loc + 1)[1:])).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


###################
# Shift functions #
###################


def rs(data, pts=0.0):
    """
    Right shift and zero fill.

    Parameters
    ----------
    data :
        Array of NMR data.
    pts : int
        Number of points to right shift.

    Returns
    -------
    ndata : ndarray
        Array of NMR data right shifted and zero filled.

    See Also
    --------
    roll : shift without zero filling.

    """
    data = np.roll(data, int(pts), axis=-1)
    data[..., :int(pts)] = 0
    return data


def ls(data, pts=0.0):
    """
    Left shift and fill with zero

    Parameters
    ----------
    data :
        Array of NMR data.
    pts : int
        Number of points to left shift.

    Returns
    -------
    ndata : ndarray
        Array of NMR data left shifted and zero filled.

    See Also
    --------
    roll : shift without zero filling.

    """
    data = np.roll(data, -int(pts), axis=-1)
    data[..., -int(pts):] = 0
    return data


def cs(data, pts=0.0, neg=False):
    """
    Circular shift

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pts : int
        Number of points to shift.  Positive value will right shift the data,
        negative values will left shift the data.
    neg : bool
        True to negate the shifted points.

    Returns
    -------
    ndata : ndarray
        Array of shifted NMR data.

    """
    return roll(data, pts, neg)


def roll(data, pts=0.0, neg=False):
    """
    Roll axis

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pts : int
        Number of points to shift.  Positive value will right shift the data,
        negative values will left shift the data.
    neg : bool
        True to negate the shifted points.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with last axis rolled.

    """
    data = np.roll(data, int(pts), axis=-1)
    if neg:
        if pts > 0:
            data[..., :pts] = -data[..., :pts]
        else:
            data[..., pts:] = -data[..., pts:]
    return data


def fsh(data, pts):
    """
    Frequency shift by Fourier transform. Negative signed phase correction.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pts : float
        Number of points to frequency shift the data.  Positive value will
        shift the spectrum to the right, negative values to the left.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with last axis rolled.

    """
    s = float(data.shape[-1])

    # inverse fft -> first order phase correction -> fft
    # idata = icomplexft(data)
    # pdata =np.exp(-2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*icomplexft(data)
    # data = complexft(pdata)

    # inplace processing
    return fft(np.exp(-2.j * pi * pts * np.arange(s) /
                        s).astype(data.dtype) * ifft(data))


def fsh2(data, pts):
    """
    Frequency Shift by Fourier transform. Postive signed phase correction.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pts : float
        Number of points to frequency shift the data.  Positive value will
        shift the spectrum to the right, negative values to the left.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with last axis rolled.

    """
    s = float(data.shape[-1])
    return fft_positive(np.exp(2.j * pi * pts * np.arange(s) /
                        s).astype(data.dtype) * ifft_positive(data))


##############
# Transforms #
##############


def nmr_reorder(data):
    """
    Reorder spectrum after FT transform to NMR order (swap halves and reverse).
    """
    s = data.shape[-1]
    return np.append(data[..., int(s / 2)::-1], data[..., s:int(s / 2):-1],
                     axis=-1)


def swap_halves(data):
    """
    Swap the halves of a spectrum,
    """
    s = data.shape[-1]
    return np.append(data[..., int(s / 2):], data[..., :int(s / 2)], axis=-1)


# Fourier based Transforms
def rft(x):
    """
    Real Fourier transform.

    Parameters
    ----------
    x : ndarray
        Array of NMR data.

    Returns
    -------
    y : ndarray
        Array of NMR data with real Fourier transform applied.

    """
    # XXX figure out what exactly this is doing...
    s = x.shape[-1]
    xp = np.zeros(x.shape, dtype="complex64")
    xp[..., 1:s / 2] = x[..., 1:-1:2] + x[..., 2::2] * 1.j
    xp[..., 0] = x[..., 0] / 2.
    xp[..., s / 2] = x[..., -1] / 2.
    return np.array(nmr_reorder(np.fft.fft(2 * xp, axis=-1).real),
                    dtype="float32")


def irft(xp):
    """
    Inverse real fourier transform

    Parameters
    ----------
    x : ndarray
        Array of NMR data.

    Returns
    -------
    y : ndarray
        Array of NMR data with an inverse real Fourier transform applied.

    """
    # XXX figure out what exactly this is doing
    s = xp.shape[-1]
    xp = np.fft.ifft(nmr_reorder(xp))   # re-order, inverse FT

    # output results
    x = np.zeros(xp.shape, dtype="float32")

    # unpack ifft data
    x[..., 1:-1:2] = xp[..., 1:s / 2].real
    x[..., 2::2] = xp[..., 1:s / 2].imag
    x[..., 0] = xp[..., 0].real
    x[..., -1] = xp[..., s / 2].real
    return x


# Fourier transforms
def fft(data):
    """
    Fourier transform, NMR ordering of results.

    There are a number of definitions of the discrete Fourier transform
    the version used in this function is as follows.

    .. math::

        A_k = \\sum_{m=0}^{n-1} a_m \\exp\\left\\{-2\\pi i{mk \\over n}
        \\right\\}\\qquad k = 0,\\ldots,n-1.

    With the inverse DFT in the :py:func:`ifft` function defined as follows.

    .. math::

        a_m = \\frac{1}{n} \\sum_{k=0}^{n-1} A_k \\exp \\left\\{2\\pi
        i{mk\\over n} \\right\\}\\qquad n = 0,\\ldots,n-1.

    Two alternative definitions are also supported by nmrglue. One in
    which both the sum in the fft and ifft are multiplied by
    :math:`\\frac{1}{\\sqrt{n}}` which results in a pair of transforms in
    which the total power contained in the the signals before and after the
    transforms are equal.  This is the type transforms used in the
    Rowland NMR Toolkit. This type of transform is performed by the
    :py:func:`fft_norm` and :py:func:`ifft_norm` functions.

    The second definition changes the sign of the exponent to be positive while
    keeping the normalization factors the same.  This type of transform is
    performed by the NMRPipe processing package and the functions
    :py:func:`fft_positive` and :py:func:`ifft_positive`.

    All of the Fourier transforms perfromed by nmrglue return results in 'NMR
    order', in which the two half of the spectrum have been swapped and
    reversed.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    ifft : Inversion Fourier transform.
    fft_norm : Norm (power) conserving Fourier transform.
    fft_positive : Forier transform with a positive exponential.


    """
    return np.fft.fftshift(np.fft.fft(data, axis=-1).astype(data.dtype), -1)


def fft_norm(data):
    """
    Fourier transform, total power preserved, NMR ordering of results

    This function is  similar to the transform performed by The Rowland NMR
    Toolkit's FFT function.

    See :py:func:`fft` for documentation of the transformation applied by this
    function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    ifft_norm : Inversion Fourier transform.
    fft : None-norm (power) conserving Fourier transform.
    fft_positive : Forier transform with a positive exponential.

    """
    return fft(data) / np.sqrt(float(data.shape[-1]))


def fft_positive(data):
    """
    Fourier transform with positive exponential, NMR ordering of results

    This is similar to the transform performed by NMRPipe's FFT function.

    See :py:func:`fft` for documentation of the transformation applied by this
    function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    ifft_positive : Inversion Fourier transform.
    fft_norm : Norm (power) conserving Fourier transform.
    fft_positive : Forier transform with a positive exponential.

    """
    # a positive exponential is the same as a IFFT, but we need to undo
    # the 1/N scaling
    s = float(data.shape[-1])
    return (np.fft.fftshift(np.fft.ifft(data, axis=-1).astype(data.dtype), -1)
                    * s)


def ifft(data):
    """
    Inverse fourier transform, NMR ordering of results.

    See :py:func:`fft` for documentation of the transformation applied by this
    function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Inverse fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    fft : Fourier transform.
    ifft_norm : Norm (power) conserving inverse Fourier transform.
    ifft_positive : inverse Forier transform with a positive exponential.

    """
    return np.fft.ifft(np.fft.ifftshift(data, -1), axis=-1).astype(data.dtype)


def ifft_norm(data):
    """
    Inverse fourier transform, total power preserved, NMR ordering of results

    This is similar to the transform performed by the Rowland NMR Toolkit's
    IFFT function.

    See :py:func:`fft` for documentation of the transformation applied by this
    function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Inverse fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    fft_norm : Norm (power) conserving Fourier transform.
    ifft : Non-norm (power) conserving inverse Fourier transform.
    ifft_positive : inverse Forier transform with a positive exponential.

    """
    return ifft(data) * np.sqrt(float(data.shape[-1]))


def ifft_positive(data):
    """
    Inverse fourier transform with positive exponential, NMR ordered results.

    This is similar to the transform performed by NMRPipe's FFT function with
    the -inv flag

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Inverse fourier transform of NMR data in 'NMR order'.

    See Also
    --------
    fft_postive : Fourier transform with a positive exponential.
    ifft_norm : Norm (power) conserving inverse Fourier transform.
    ifft : Inverse Forier transform with a negative exponential.

    """
    # a inverse fft with positive exponential in the FFT definition is the
    # same as a FFT with negative exponentials, but with a 1/N scaling factor
    s = 1.0 / float(data.shape[-1])
    return (np.fft.fft(np.fft.ifftshift(data, -1), axis=-1).astype(data.dtype)
             * s)


# Hadamard Transform functions
def int2bin(n, digits=8):
    """
    Integer to binary string
    """
    return "".join([str((n >> y) & 1) for y in range(digits - 1, -1, -1)])


def bin2int(s):
    """
    binary string to integer
    """
    o = 0
    k = len(s) - 1
    for i, v in enumerate(s):
        o = o + int(v) * 2 ** (k - i)
    return o


def gray(n):
    """
    Calculate n-bit gray code
    """
    g = [0, 1]
    for i in range(1, int(n)):
        mg = g + g[::-1]   # mirror the current code
        # first bit 0/2**u for mirror
        first = [0] * 2 ** (i) + [2 ** (i)] * 2 ** (i)
        g = [mg[j] + first[j] for j in range(2 ** (i + 1))]
    return g


def ha(data):
    """
    Hadamard Transform

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Hadamard transform of NMR data.

    Notes
    -----
    This function is very slow.  Implement a Fast Walsh-Hadamard Transform
    with sequency/Walsh ordering (FWHT_w) will result in much faster tranforms.

    http://en.wikipedia.org/wiki/Walsh_matrix
    http://en.wikipedia.org/wiki/Fast_Hadamard_transform

    """
    # implementation is a proof of concept and EXTEMEMLY SLOW

    # determind the order and final size of input vectors
    ord = int(np.ceil(np.log2(data.shape[-1])))  # Walsh/Hadamard order
    max = 2 ** ord

    # zero fill to power of 2
    pad = max - data.shape[-1]
    zdata = zf(data, pad)

    # Multiple each vector by the hadamard matrix
    nat = np.zeros(zdata.shape, dtype=zdata.dtype)
    H = scipy.linalg.hadamard(max)
    nat = np.dot(zdata, H)
    nat = np.array(nat, dtype=data.dtype)

    # Bit-Reversal Permutation
    s = [int2bin(x, digits=ord)[::-1] for x in range(max)]
    brp = [bin2int(x) for x in s]
    brp_data = np.take(nat, brp, axis=-1)

    # Gray code permutation (bit-inverse)
    gp = gray(ord)
    gp_data = np.take(brp_data, gp, axis=-1)

    return gp_data


def ht(data, N=None):
    """
    Hilbert transform.

    Reconstruct imaginary data via hilbert transform.

    Parameters
    ----------
    data : ndarrat
        Array of NMR data.
    N : int or None
        Number of Fourier components.

    Returns
    -------
    ndata : ndarray
        NMR data which has been Hilvert transformed.

    """
    # XXX come back and fix this when a sane version of scipy.signal.hilbert
    # is included with scipy 0.8

    # create an empty output array
    fac = N / data.shape[-1]
    z = np.empty(data.shape, dtype=(data.flat[0] + data.flat[1] * 1.j).dtype)
    if data.ndim == 1:
        z[:] = scipy.signal.hilbert(data.real, N)[:data.shape[-1]] * fac
    else:
        for i, vec in enumerate(data):
            z[i] = scipy.signal.hilbert(vec.real, N)[:data.shape[-1]] * fac

    # correct the real data as sometimes it changes
    z.real = data.real
    return z


##########################
# Standard NMR Functions #
##########################


def di(data):
    """
    Delete imaginary from data

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Array of NMR data without imaginaries.

    """
    return data.real


def ps(data, p0=0.0, p1=0.0, inv=False):
    """
    Linear phase correction

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    p0 : float
        Zero order phase in degrees.
    p1 : float
        First order phase in degrees.
    inv : bool, optional
        True for inverse phase correction

    Returns
    -------
    ndata : ndarray
        Phased NMR data.

    """
    p0 = p0 * pi / 180.  # convert to radians
    p1 = p1 * pi / 180.
    size = data.shape[-1]
    apod = np.exp(1.0j * (p0 + (p1 * np.arange(size) / size))
                  ).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


def ps_exp(data, p0=0.0, tc=0.0, inv=False):
    """
    Exponential Phase Correction

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    p0 : float
        Zero order phase in degrees.
    tc : float
        Exponential decay constant.
    inv : bool, optional
        True for inverse phase correction

    Returns
    -------
    ndata : ndarray
        Phased NMR data.

    """
    p0 = p0 * pi / 180.  # convert to radians
    size = data.shape[-1]
    apod = np.exp(1.0j * (p0 * np.exp(-tc * np.arange(size) / size))
                  ).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data


def tp(data, hyper=False):
    """
    Transpose data.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    hyper : bool
        True if hypercomplex data.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with axes transposed.

    """
    if hyper:
        return tp_hyper(data)
    else:
        return data.transpose()


ytp = tp


xy2yx = tp


def tp_hyper(data):
    """
    Hypercomplex tranpose.

    Use when both dimension are complex.

    Parameters
    ----------
    data : ndarray
        Array of hypercomplex NMR data.

    Returns
    -------
    ndata : ndarray
        Array of hypercomplex NMR data with axes transposed.
    """
    return c2ri(ri2c(data).transpose())


def zf_inter(data, pts=1):
    """
    Zero fill between points.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pts : int
        Number zeros to add between points.

    Returns
    -------
    ndata : ndarray
        Array of NMR data to which `pts` zero have been added between all
        points.

    """
    size = list(data.shape)
    size[-1] = (pts + 1) * size[-1]
    z = np.zeros(size, dtype=data.dtype)
    z[..., ::pts + 1] = data[..., :]
    return z


def zf_pad(data, pad=0, mid=False):
    """
    Zero fill by padding with zeros.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    pad : int
        Number of zeros to pad data with.
    mid : bool
        True to zero fill in middle of data.

    Returns
    -------
    ndata : ndarray
        Array of NMR data to which `pad` zeros have been appended to the end or
        middle of the data.

    """
    size = list(data.shape)
    size[-1] = pad
    z = np.zeros(size, dtype=data.dtype)

    if mid:
        h = int(data.shape[-1] / 2.0)
        return np.concatenate((data[..., :h], z, data[..., h:]), axis=-1)
    else:
        return np.concatenate((data, z), axis=-1)


zf = zf_pad


def zf_double(data, n, mid=False):
    """
    Zero fill by doubling original data size once or multiple times.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    n : int
        Number of times to double the size of the data.
    mid : bool
        True to zero fill in the middle of data.

    Returns
    -------
    ndata : ndarray
        Zero filled array of NMR data.

    """
    return zf_pad(data, int((data.shape[-1] * 2 ** n) - data.shape[-1]), mid)


def zf_size(data, size, mid=False):
    """
    Zero fill to given size.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    size : int
        Size of data after zero filling.
    mid : bool
        True to zero fill in the middle of data.

    Returns
    -------
    ndata : ndarray
        Zero filled array of NMR data.

    """
    return zf_pad(data, int(size - data.shape[-1]), mid)


def largest_power_of_2(value):
    """
    Find the nearest power of two equal to or larger than a value.

    Parameters
    ----------
    value : int
        Value to find nearest power of two equal to or larger than.

    Returns
    -------
    pw : int
        Power of 2.

    """
    return int(pow(2, np.ceil(np.log(value) / np.log(2))))


def zf_auto(data, mid=False):
    """
    Zero fill to next largest power of two.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    mid : bool
        True to zero fill in the middle of data.

    Returns
    -------
    ndata : ndarray
        Zero filled array of NMR data.

    """
    return zf_size(data, largest_power_of_2(data.shape[-1]), mid)


####################
# Basic Untilities #
####################


# Add Constant
def add(data, r=0.0, i=0.0, c=0.0):
    """
    Add constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    r : float
        Constant to add to real data.
    i : float
        Constant to add to imaginary data, ignored if no data does not have
        imaginary data.
    c : float
        Constant to add to both real and imaginary data.  This constant is only
        added to the real data if no imaginary data is present.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with constant added.

    """
    data.real = data.real + r + c
    if np.iscomplexobj(data):
        data.imag = data.imag + i + c
    return data


def add_ri(data):
    """
    Add real and imaginary components.

    Parameters
    ----------
    data : ndarray
        Array of complex NMR data.

    Returns
    -------
    ndata : ndarray
        Sum of real and imaginary component of NMR data.

    """
    return data.real + data.imag


# Derivative
def dx(data):
    """
    Derivative by central difference

    Edges are takes as difference between nearest points

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Derivate of NMR data.

    """
    z = np.empty_like(data)
    z[..., 0] = data[..., 1] - data[..., 0]    # first point
    z[..., -1] = data[..., -1] - data[..., -2]  # last point
    z[..., 1:-1] = data[..., 2:] - data[..., :-2]  # interior
    return z


# Extract Region
def ext(data, x0=None, xn=None, y0=None, yn=None):
    """
    Extract region.

    Data should be two dimensional. Axes are labeled Y, X.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    x0 : int
        X-axis extract region start.
    xn : int
        X-axis extract region end.
    y0 : int
        Y-axis extract region start.
    yn : int
        Y-axis extract region end.

    Returns
    -------
    ndata : ndarray
        Extracted region.

    """
    return data[y0:yn, x0:xn]


def ext_left(data):
    """
    Extract the left half of spectrum.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Left half of spectrum.

    """
    return data[..., 0:int(data.shape[-1] / 2.)]


def ext_right(data):
    """
    Extract the right half of the spectrum.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Right half of spectrum.

    """
    return data[..., int(data.shape[-1] / 2.):]


def ext_mid(data):
    """
    Extract the middle of the spectrum.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Middle half of spectrum.

    """
    return data[..., int(data.shape[-1] * 1. / 4.):
                                int(data.shape[-1] * 3. / 4.)]


# Integrate
def integ(data):
    """
    Integrate by cumulative sum.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.

    Returns
    -------
    ndata : ndarray
        Integrated NMR data.

    """
    return np.cumsum(data, axis=-1)


# Modulus Calculation
def mc(data):
    """
    Modulus calculation.

    Calculates sqrt(real^2 + imag^2)
    """
    return np.sqrt(data.real ** 2 + data.imag ** 2)


def mc_pow(data):
    """
    Modulus calculation. Squared version.

    Calculated real^2+imag^2
    """
    return data.real ** 2 + data.imag ** 2


# Mirror
def mir_left(data):
    """
    Append a mirror image of the data on the left.
    """
    return np.append(data, data[..., ::-1], axis=-1)


def mir_right(data):
    """
    Append a mirror image of the data on the right.
    """
    return np.append(data[..., ::-1], data, axis=-1)


def mir_center(data):
    """
    Append a mirror image of the data in the center.
    """
    s = data.shape[-1]
    return np.concatenate((data[..., s / 2:], data, data[..., :s / 2]),
                          axis=-1)


def mir_center_onepoint(data):
    """
    Append a mirror image of the data in the center with a one point shift
    amd negate appended imaginary data.
    """
    s = int(data.shape[-1])
    data = np.concatenate((data[..., s - 1:0:-1], data), axis=-1)
    if np.iscomplexobj(data):
        data.imag[..., :s - 1] = -data.imag[..., :s - 1]
    return data


# Multiply by a constant
def mult(data, r=1.0, i=1.0, c=1.0):
    """
    Multiply by a constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    r : float
        Constant to multiply real channel by.
    i : float
        Constant to multiply imaginary channel by. If data has no imaginary
        component, this parameter is ignored.
    c : float
        Constant to multiply both channels by.

    Returns
    -------
    ndata : ndarray
        Array of NMR data which has been multiplied by a constant.

    """
    data.real = data.real * r * c
    if np.iscomplexobj(data):
        data.imag = data.imag * i * c
    return data


# Reverse
def rev(data):
    """
    Reverse data.
    """
    return data[..., ::-1]


# Set to a Constant
def set(data, c):
    """
    Set to a constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    c : float or complex
        Constant to set data to.

    Returns
    -------
    ndata : ndarray
        Array of constants.

    Notes
    -----
    data is set in place, if original data is required make a copy before
    calling this function.

    """
    data[..., :] = c
    return data


def set_complex(data, v):
    """
    Set real and imaginary portions of data to a constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    v : float
        Constant to set both real and imaginary component of `data` to.

    Returns
    -------
    ndata : ndarray
        Array with real and imaginary components set to v.

    Notes
    -----
    data is set in place, if original data is required make a copy before
    calling this function.

    """
    data.real = v
    if np.iscomplexobj(data):
        data.imag = v
    return data


def set_real(data, v):
    """
    Set real component of data to a constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    v : float
        Constant to set both real component of `data` to.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with real components set to v.

    Notes
    -----
    data is set in place, if original data is required make a copy before
    calling this function.

    """
    data.real = v
    return data


def set_imag(data, v):
    """
    Set imaginary portion of data to a constant.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    v : float
        Constant to set both imaginary component of `data` to.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with imaginary components set to v.

    Notes
    -----
    data is set in place, if original data is required make a copy before
    calling this function.

    """
    if np.iscomplexobj(data):
        data.imag = v
    return data


# Shuffle Utilities
def ri2c(data):
    """
    Interleave real and imaginary data into a real array.
    """
    s = list(data.shape)
    s[-1] = s[-1] * 2
    n = np.empty(s, data.real.dtype)
    n[..., ::2] = data.real
    n[..., 1::2] = data.imag
    return n


def interleave_complex(data):
    """
    Unpack complex data into an interleaved real, imaginary array.
    """
    return ri2c(data)


def unpack_complex(data):
    """
    Unpacks complex array into real array (interleaves values).
    """
    return ri2c(data)


def c2ri(data):
    """
    Seperate interleaved real, imaginary data into complex array.

    Assumes data is real only, ignores imaginary portion of data.

    """
    # make a 1,1 array to determind dtype
    temp = np.array(data.flat[0] + data.flat[1] * 1j)
    s = list(data.shape)
    s[-1] = int(s[-1] / 2)
    n = np.empty(s, temp.dtype)
    del(temp)
    n.real = data.real[..., ::2]
    n.imag = data.real[..., 1::2]
    return n


def seperate_interleaved(data):
    """
    Seperate interleaved real, imaginary data into complex array.
    """
    return c2ri(data)


def pack_complex(data):
    """
    Packs interleaved real array into complex array.
    """
    return c2ri(data)


def decode_States(data):
    """
    Decode data collected using States (seperates interleaved data).
    """
    return c2ri(data)


def ri2rr(data):
    """
    Append imaginary data to end of real data, returning a real array.
    """
    s = list(data.shape)
    half = int(s[-1])
    s[-1] = half * 2
    n = np.empty(s, data.real.dtype)
    n[..., :half] = data.real
    n[..., half:] = data.imag
    return n


append_imag = ri2rr


def rr2ri(data):
    """
    Unappend real and imaginary data returning a complex array.
    """
    # make a 1,1 array to determind dtype
    temp = np.array(data.flat[0] + data.flat[1] * 1.j)
    s = list(data.shape)
    half = int(s[-1] / 2.0)
    s[-1] = half
    n = np.empty(s, temp.dtype)
    del(temp)
    n.real = data[..., :half]
    n.imag = data[..., half:]
    return n


unappend_imag = rr2ri


def exlr(data):
    """
    Exchange left and right halves of array.
    """
    half = int(data.shape[-1] / 2)
    n = np.empty_like(data)
    n[..., :half] = data[..., half:]
    n[..., half:] = data[..., :half]
    return n


exchange_lr = exlr


def rolr(data):
    """
    Rotate left and right halves of array.
    """
    half = int(data.shape[-1] / 2)
    n = np.empty_like(data)
    n[..., :half] = data[..., (half - 1)::-1]
    n[..., half:] = data[..., :(half - 1):-1]
    return n


rotate_lr = rolr


def swap(data):
    """
    Swap real and imaginary data.
    """
    n = np.empty_like(data)
    n.real = data.imag
    n.imag = data.real
    return n


swap_ri = swap


def bswap(data):
    """
    Byteswap data
    """
    return data.byteswap()


byte_swap = bswap


# Sign Manipulation Utilities
def neg_left(data):
    """
    Negate left half.
    """
    data[..., 0:int(data.shape[-1] / 2.)] = \
            -data[..., 0:int(data.shape[-1] / 2.)]
    return data


def neg_right(data):
    """
    Negate right half.
    """
    data[..., int(data.shape[-1] / 2.):] = \
            -data[..., int(data.shape[-1] / 2.):]
    return data


def neg_middle(data):
    """
    Negate middle half.
    """
    data[..., int(data.shape[-1] * 1. / 4.):int(data.shape[-1] * 3. / 4.)] = \
        -data[..., int(data.shape[-1] * 1. / 4.):int(data.shape[-1] * 3. / 4.)]
    return data


def neg_edges(data):
    """
    Negate edge half (non-middle) of spectra.
    """
    data[..., :int(data.shape[-1] * 1. / 4)] = \
            -data[..., :int(data.shape[-1] * 1. / 4)]
    data[..., int(data.shape[-1] * 3. / 4):] = \
            -data[..., int(data.shape[-1] * 3. / 4):]
    return data


def neg_all(data):
    """
    Negate data
    """
    return -data


def neg_real(data):
    """
    Negate real data
    """
    data.real = -data.real
    return data


def neg_imag(data):
    """
    Negate imaginary data
    """
    data.imag = -data.imag
    return data


def neg_even(data):
    """
    Negate even points
    """
    data[..., ::2] = -data[..., ::2]
    return data


def neg_odd(data):
    """
    Negate odd points
    """
    data[..., 1::2] = -data[..., 1::2]
    return data


def neg_alt(data):
    """
    Negate alternate (odd) points.
    """
    return neg_odd(data)


def abs(data):
    """
    Replace data with absolute value of data (abs of real, imag seperately)
    """
    data.real = np.abs(data.real)
    data.imag = np.abs(data.imag)
    return data


def sign(data):
    """
    Replace data with sign (-1 or 1) of data (seperately on each channel)
    """
    data.real = np.sign(data.real)
    data.imag = np.sign(data.imag)
    return data


##################
# Misc Functions #
##################


# Coadd data
def coadd(data, clist, axis=-1):
    """
    Coadd data.

    Reduce data along axis by blocks data and multiplying by
    coefficients in clist.  Incomplete blocks are discarded.

    Parameters
    ----------
    data : 2D ndarray
        Array of NMR data data.
    clist : list
        List of Coefficients
    axis : {0, 1, -1}
        Axis to reduce. 0 for Y, 1 or -1 for X.

    """
    # there is probably a efficient way to do this with tile and inner
    # or scipy.ndimage.generic_filter

    # algorith creates a empty array, then fills it element wise
    # with each factor from clist times the blocks selected

    s = list(data.shape)    # data shape
    k = len(clist)          # length of coefficient list

    if axis == 1 or axis == -1:   # 'x' axis
        s[-1] = np.floor(float(s[-1]) / k)
        n = np.zeros(s, dtype=data.dtype)
        m = s[-1] * k   # last element read
        for i in range(k):
            n = n + clist[i] * data[..., i:m:k]
    else:   # 'y' axis
        s[0] = np.floor(float(s[0]) / k)
        n = np.zeros(s, dtype=data.dtype)
        m = s[0] * k
        for i in range(k):
            n = n + clist[i] * data[i:m:k]
    return n


coad = coadd


# Image Processing
def thres(data, thres=0.0):
    """
    Mark values less than thres as invalid (for use with filters)

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    thres : float
        Threshold value.

    Returns
    -------
    ndata : masked ndarray
        Masked array of NMR data with values less than thres masked.

    """
    return np.ma.masked_less(data, thres)


def conv(data, kern=[1.], m="wrap", c=0.0):
    """
    Convolute data with kernel.

    Real and imaginary components of data are convolved seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    kern : array_like
        Convolution kernel
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of convolution are dealt with.
    c : float
        Constant value for use when m is  'constant'.

    Returns
    -------
    ndata : ndarray
        NMR data which has been convolved with the given kernel.

    """
    kern = np.array(kern)
    data.real = scipy.ndimage.convolve(data.real, weights=kern, mode=m, cval=c)
    data.imag = scipy.ndimage.convolve(data.imag, weights=kern, mode=m, cval=c)
    return data


convolute = conv


def corr(data, kern=[1.], m="wrap", c=0.0):
    """
    Correlate data with a kernel (weights).

    Real and imaginary components of data are correlated seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    kern : array_like
        Correlation kernel (weights).
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of correlation are dealt with.
    c : float
        Constant value for use when m is  'constant'.

    Returns
    -------
    ndata : ndarray
        NMR data which has been correlated with the given kernel.

    """
    kern = np.array(kern)
    data.real = scipy.ndimage.correlate(data.real, weights=kern, mode=m,
                                        cval=c)
    data.imag = scipy.ndimage.correlate(data.imag, weights=kern, mode=m,
                                        cval=c)
    return data


correlate = corr


def filter_median(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a median filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    data.real = scipy.ndimage.median_filter(data.real, size=s, mode=m, cval=c)
    data.imag = scipy.ndimage.median_filter(data.imag, size=s, mode=m, cval=c)
    return data


def filter_min(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a minimum filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    data.real = scipy.ndimage.minimum_filter(data.real, size=s, mode=m, cval=c)
    data.imag = scipy.ndimage.minimum_filter(data.imag, size=s, mode=m, cval=c)
    return data


def filter_max(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a maximum filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    data.real = scipy.ndimage.maximum_filter(data.real, size=s, mode=m, cval=c)
    data.imag = scipy.ndimage.maximum_filter(data.imag, size=s, mode=m, cval=c)
    return data


def filter_percentile(data, percentile, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a percentile filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    percentile : float
        Filter percentile parameter.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    data.real = scipy.ndimage.percentile_filter(data.real, percentile, size=s,
                                                mode=m, cval=c)
    data.imag = scipy.ndimage.percentile_filter(data.imag, percentile, size=s,
                                                mode=m, cval=c)
    return data


def filter_rank(data, rank, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a rank filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    rank : int
        Filter rank parameter.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    data.real = scipy.ndimage.rank_filter(data.real, rank, size=s, mode=m,
                                                cval=c)
    data.imag = scipy.ndimage.rank_filter(data.imag, rank, size=s, mode=m,
                                                cval=c)
    return data


# These filter are much slower as they use the generic filter functions...


def filter_amin(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply an absolute minimum filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = amin_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_amax(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply an absolute maximum filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = amax_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_range(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a range filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = range_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_avg(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply an average filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = avg_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_dev(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a standard deviation filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = std_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_sum(data, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a summation filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = sum_flt
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


def filter_generic(data, filter, s=(1, 1), m="wrap", c=0.0):
    """
    Apply a generic filter.

    Real and imaginary components are filtered seperately.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    filter : functions
        Python function which takes an array and returns a single value.
    s : tuple
        Shape or size taken for each step of the filter.
    m : {'reflect', 'constant', 'nearest', 'wrap'}
        How edges of are dealt with.
    c : float
        Constant value for use when m is 'constant'.

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.

    """
    flt = filter
    data.real = scipy.ndimage.generic_filter(data.real, flt, size=s, mode=m,
                                             cval=c)
    data.imag = scipy.ndimage.generic_filter(data.imag, flt, size=s, mode=m,
                                             cval=c)
    return data


# filter functions

def amin_flt(arr):
    return  arr[np.abs(arr).argmin()]


def amax_flt(arr):
    return arr[np.abs(arr).argmax()]


def range_flt(arr):
    return arr.max() - arr.min()


def avg_flt(arr):
    return arr.avg()


def std_flt(arr):
    return arr.std()


def sum_flt(arr):
    return arr.sum()


# Scale Quad Artifacts
def qart(data, a=0.0, f=0.0):
    """
    Scale Quad Artifacts.

    Replaces imaginary data with (1 + a) * data.imag + f * data.real

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    a : float
        Amplitude adjustment.
    f : float
        Phase adjustment.

    Returns
    -------
    ndata: ndarray
        Array of NMR data with quadrature artifacts scaled.

    """
    data.imag = (1 + a) * data.imag + f * data.real
    return data


def qart_auto(data):
    """
    Scale quad artifacts by values from Gram-Schmidt orthogonalization.
    """
    a, f = gram_schmidt(data)
    return qart(data, a, f)


def gram_schmidt(data):
    """
    Calculate Gram-Schmidt orthogonalization parameters.
    """
    # method similar to Hock and Stern, "NMR Data Processing" p.61
    # sum of correlation between data.real and data.imag
    C = (data.real * data.imag).sum()

    # total power in real channel
    R = (data.real * data.real).sum()

    # remove correlation from imag channel
    idata = data.imag - (C / R) * data.real

    # total power in uncorrelated imag channel
    S = (idata * idata).sum()

    # imag(data'') = R/S*imag(data')
    # imag(data')  = imag(data)-C/R * real(data)
    # therefore:
    # imag(data'') = R/S*imag(data) - R*C/(S*R) * real(data)
    # so A = R/S, B=-C/(S)
    return(R / S, -C / S)


# Complex Mixing
def qmix(data, carr):
    """
    Mix input and output channels provided coefficient array.

    Parameters
    ----------
    data : 2D ndarray
        Array of NMR data.
    carr : array_liek
        Array of coefficients for mixing.  The size of carr must evenly divide
        qmix.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with channels mixed.

    """
    carr = np.array(carr).transpose()
    ic = carr.shape[1]  # input channels
    oc = carr.shape[0]  # output channels

    if data.shape[0] % ic != 0 or data.shape[0] % oc != 0:
        raise ValueError("Coefficient array does not evenly divide data")

    # create an empty blank output array
    s = list(data.shape)
    s[0] = s[0] * float(oc) / float(ic)
    n = np.empty(s, data.dtype)

    # remix each block
    for i in xrange(int(data.shape[0] / float(ic))):
        block = data[i * ic:(i + 1) * ic]
        n[i * oc:(i + 1) * oc] = np.dot(carr, block)

    return n


# Smooth and Center
def smo(data, n):
    """
    Smooth data.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    n : int
        Size of smoothing window (+/- points)

    Returns
    -------
    ndata : ndarray
        Array of smoothed NMR data.

    """
    # XXX this can probably be accomplished by a median_filter
    n = int(n)

    # a is the accumulator
    a = np.copy(data)

    for i in range(1, n + 1):
        a = a + rs(data, i) + ls(data, i)

    # divide the interior by 2*n+1 to get mean
    a[..., n:-n] = a[..., n:-n] / (2 * n + 1)
    # divide the left edges by 2n+1-i where i is the distance from the interior
    for i in range(1, n + 1):
        a[..., n - i] = a[..., n - i] / (2. * n + 1 - i)
    # divide the right edge similarly
    for i in range(-n, 0):
        a[..., i] = a[..., i] / (n - i)
    return a


def center(data, n):
    """
    Center data.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    n : int
        Size of centering window (+/- points)

    Returns
    -------
    ndata : ndarray
        Array of centered NMR data.

    """
    return data - smo(data, n)


def zd(data, window, x0=0.0, slope=1.0):
    """
    Zero Diagonal band with generic window function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    window : ndarray
        Window to apply to diagonal band.
    wide : int
        Diagonal band half width in points.
    x0 : int
        Starting location of diagonal band in points.
    slope : float
        Slope of diagonal band.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with diagonal band set to zero.


    """
    width = len(window)         # full width
    wide = (width - 1.) / 2     # half width
    rows = data.shape[0]        # rows in data
    cols = data.shape[-1]       # columns in data
    c_start = x0 + slope        # start of center diagonal band

    # last row to apply window to is last row or where we run off the grid
    max_r = int(min(rows, np.floor((cols - c_start + wide) / slope) + 1))

    # apply window to band row by row
    for r in xrange(max_r):  # r from 0 to max_r-1
        w_min = 0           # window min
        w_max = width       # window max

        c_mid = int(r * slope + (c_start))  # middle of diagonal band
        c_min = c_mid - wide
        c_max = c_mid + wide + 1

        if c_min < 0:
            w_min = int(-c_min)
            c_min = 0
        if c_max > cols:
            w_max = int(w_max - (c_max - cols))
            c_max = cols

        data[r, c_min:c_max] = data[r, c_min:c_max] * window[w_min:w_max]

    return data


def zd_boxcar(data, wide=1.0, x0=0.0, slope=1.0):
    """
    Zero diagonal band with a boxcar function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    wide : int
        Diagonal band half width in points.
    x0 : int
        Starting location of diagonal band in points.
    slope : float
        Slope of diagonal band.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with diagonal band set to zero.

    """
    window = np.zeros(2 * wide + 1)
    return zd(data, window, x0=x0, slope=slope)


def zd_triangle(data, wide=1.0, x0=0.0, slope=1.0):
    """
    Zero diagonal band with triangle function

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    wide : int
        Diagonal band half width in points.
    x0 : int
        Starting location of diagonal band in points.
    slope : float
        Slope of diagonal band.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with diagonal band set to zero.

    """
    window = np.append(np.linspace(1, 0, wide + 1),
                    np.linspace(0, 1, wide + 1)[1:])
    return zd(data, window, x0=x0, slope=slope)


def zd_sinebell(data, wide=1.0, x0=0.0, slope=1.0):
    """
    Zero diagonal band with sinebell function.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    wide : int
        Diagonal band half width in points.
    x0 : int
        Starting location of diagonal band in points.
    slope : float
        Slope of diagonal band.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with diagonal band set to zero.

    """
    window = 1 - np.sin(np.linspace(0, pi, 2 * wide + 1))
    return zd(data, window, x0=x0, slope=slope)


def zd_gaussian(data, wide=1.0, x0=0.0, slope=1.0, g=1):
    """
    Zero Diagonal band with gaussian function

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    wide : int
        Diagonal band half width in points.
    x0 : int
        Starting location of diagonal band in points.
    slope : float
        Slope of diagonal band.
    g : float
        Width of Gaussian function.

    Returns
    -------
    ndata : ndarray
        Array of NMR data with diagonal band set to zero.

    """
    tln2 = np.sqrt(2 * np.log(2))
    window = 1 - scipy.signal.gaussian(2 * wide + 1, g / tln2)
    return zd(data, window, x0=x0, slope=slope)
