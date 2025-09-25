"""
Linear Prediction (LP) functions for extrapolating and modeling NMR signals.
"""


__developer_info__ = """
Notes
^^^^^
This module contains functions for performing linear prediction on NMR  data.
The algorithms used were selected for simplicity to show how linear prediction
works not for computational speed nor stability.  Locations where significant
improvements can be made to improve speed or stability are indicated with
SPEED and STABILITY within the source code with discussion following.

The notation for the Linear Prediction equation, coefficients, roots, etc.
closely match those in "NMR Data Processing" by Hoch snd Stern. This book was
references for many of the algorithms in this module.

Reduced order LP-SVD and LP-TLS methods are not implemented but should
be easy to add if desired.  See :py:func:`find_lproots_hsvd` for an example.

"""

import numpy as np
import scipy.linalg


#############################################
# linear prediction extrapolation functions #
#############################################


def lp(data, pred=1, slice=slice(None), order=8, mode="f", append="after",
       bad_roots="auto", fix_mode="on", mirror=None, method="svd"):
    """
    Linear prediction extrapolation of 1D or 2D data.

    Parameters
    ----------
    data : ndarray
        1D or 2D NMR data with the last (-1) axis in the time domain.
    pred : int
        Number of points to predict along the last axis.
    slice : slice object, optional
        Slice object which selects the region along the last axis to use in LP
        equation.  The default (slice(None)) will use all points.
    order : int
        Prediction order, number of LP coefficients calculated.
    mode : {'f', 'b', 'fb' or 'bf'}
        Mode to generate LP filter. 'f' for forward,'b' for backward, fb for
        'forward-backward and 'bf' for backward-forward.
    append : {'before', 'after'}
        Location to append the data, either 'before' the current data, or
        'after' the existing data. This is independent of the `mode` parameter.
    bad_roots : {'incr', 'decr', None, 'auto'}
        Type of roots which to consider bad and to stabilize.  Option are those
        with increasing signals 'incr' or decreasing signals 'decr'.  None will
        perform no root stabilizing.  The default ('auto') will set the
        parameter based on the `mode` parameter.  'f' or 'fb' `mode` will
        results in a 'incr' `bad_roots` parameter, 'b' or 'bf` in 'decr'
    fix_mode : {'on', 'reflect'}
        Method used to stabilize bad roots, 'on' to move the roots onto the
        unit circle, 'reflect' to reflect bad roots across the unit circle.
        This parameter is ignored when `bad_roots` is None.
    mirror : {None, '0', '180'}
        Mode to form mirror image of data before processing.  None will
        process the data trace as provided (no mirror image). '0' or '180'
        forms a mirror image of the sliced trace to calculate the LP filter.
        '0' should be used with data with no delay, '180' with data
        with an initial half-point delay.
    method : {'svd', 'qr', 'choleskey', 'tls'}
        Method to use to calculate the LP filter. Choices are a SVD ('svd'), QR
        ('qr'), or Choleskey ('choleskey') decomposition, or Total Least
        Squares ('tls').

    Returns
    -------
    ndata : ndarray
        NMR data with `pred` number of points linear predicted and appended to
        the original data.

    Notes
    -----
    When given 2D data a series of 1D linear predictions are made to
    each row in the array, extending each by pred points. To perform a 2D
    linear prediction using a 2D prediction matrix use :py:func:`lp2d`.

    In forward-backward or backward-forward mode root stabilizing is done
    on both sets of signal roots as calculated in the first mode direction.
    After averaging the coefficient the roots are again stabilized.

    When the append parameter does not match the LP mode, for example
    if a backward linear prediction (mode='b') is used to predict points
    after the trace (append='after'), any root fixing is done before reversing
    the filter.

    """
    if data.ndim == 1:
        return lp_1d(data, pred, slice, order, mode, append, bad_roots,
                     fix_mode, mirror, method)
    elif data.ndim == 2:
        # create empty array to hold output
        s = list(data.shape)
        s[-1] = s[-1] + pred
        new = np.empty(s, dtype=data.dtype)
        # vector-wise 1D LP
        for i, trace in enumerate(data):
            new[i] = lp_1d(trace, pred, slice, order, mode, append, bad_roots,
                           fix_mode, mirror, method)
        return new
    else:
        raise ValueError("data must be a one or two dimensional array")


# method specific extrapolation
def lp_svd(data, pred=1, slice=slice(None), order=8, mode="f", append="after",
           bad_roots="auto", fix_mode="on", mirror=None):
    """
    Linear Prediction extrapolation of 1D or 2D data using SVD decomposition.

    See :py:func:`lp` for documentation.

    """
    return lp(data, pred, slice, order, mode, append, bad_roots, fix_mode,
              mirror, method="svd")


def lp_qr(data, pred=1, slice=slice(None), order=8, mode="f", append="after",
          bad_roots="auto", fix_mode="on", mirror=None):
    """
    Linear Prediction extrapolation of 1D or 2D data using QR decomposition.

    See :py:func:`lp` for documentation

    """
    return lp(data, pred, slice, order, mode, append, bad_roots, fix_mode,
              mirror, method="qr")


def lp_cho(data, pred=1, slice=slice(None), order=8, mode="f", append="after",
           bad_roots="auto", fix_mode="on", mirror=None):
    """
    Linear Prediction extrapolation of 1D or 2D data using Cholesky
    decomposition.

    See :py:func:`lp` for documentation
    """
    return lp(data, pred, slice, order, mode, append, bad_roots, fix_mode,
              mirror, method="cholesky")


def lp_tls(data, pred=1, slice=slice(None), order=8, mode="f", append="after",
           bad_roots="auto", fix_mode="on", mirror=None):
    """
    Linear Prediction extrapolation of 1D or 2D data using Total Least Squares.

    See :py:func:`lp` for documentation.

    """
    return lp(data, pred, slice, order, mode, append, bad_roots, fix_mode,
              mirror, method="tls")


# underlying 1D extrapolation
def lp_1d(trace, pred=1, slice=slice(None), order=8, mode="f", append="after",
          bad_roots="auto", fix_mode="on", mirror=None, method="svd"):
    """
    Linear Prediction extrapolation of 1D data.

    Parameters
    ----------
    trace : ndarray
        1D NMR data in the time domain.
    pred : int
        Number of points to predict along the last axis.
    slice : slice object, optional
        Slice object which selects the region along the last axis to use in LP
        equation.  The default (slice(None)) will use all points.
    order : int
        Prediction order, number of LP coefficients calculated.
    mode : {'f', 'b', 'fb' or 'bf'}
        Mode to generate LP filter. 'f' for forward,'b' for backward, fb for
        'forward-backward and 'bf' for backward-forward.
    append : {'before', 'after'}
        Location to append the data, either 'before' the current data, or
        'after' the existing data. This is independent of the `mode` parameter.
    bad_roots : {'incr', 'decr', None, 'auto'}
        Type of roots which to consider bad and to stabilize.  Option are those
        with increasing signals 'incr' or decreasing signals 'decr'.  None will
        perform no root stabilizing.  The default ('auto') will set the
        parameter based on the `mode` parameter.  'f' or 'fb' `mode` will
        results in a 'incr' `bad_roots` parameter, 'b' or 'bf` in 'decr'
    fix_mode : {'on', 'reflect'}
        Method used to stabilize bad roots, 'on' to move the roots onto the
        unit circle, 'reflect' to reflect bad roots across the unit circle.
        This parameter is ignored when `bad_roots` is None.
    mirror : {None, '0', '180'}
        Mode to form mirror image of data before processing.  None will
        process the data trace as provided (no mirror image). '0' or '180'
        forms a mirror image of the sliced trace to calculate the LP filter.
        '0' should be used with data with no delay, '180' with data
        with an initial half-point delay.
    method : {'svd', 'qr', 'choleskey', 'tls'}
        Method to use to calculate the LP filter. Choices are a SVD ('svd'), QR
        ('qr'), or Choleskey ('choleskey') decomposition, or Total Least
        Squares ('tls').

    Returns
    -------
    ntrace : ndarray
        NMR data with `pred` number of points linear predicted and appended to
        the original data.

    Notes
    -----
    In forward-backward or backward-forward mode root stabilizing is done
    on both sets of signal roots as calculated in the first mode direction.
    After averaging the coefficient the roots are again stabilized.

    When the append parameter does not match the LP mode, for example
    if a backward linear prediction (mode='b') is used to predict points
    after the trace (append='after'), any root fixing is done before reversing
    the filter.


    See Also
    --------
    lp : 1D or 2D linear prediction extrapolation.

    """
    # check for bad arguments
    if mode not in ["f", "b", "fb", "bf"]:
        raise ValueError("mode must be 'f', 'b', 'fb', or 'bf'")
    if append not in ["before", "after"]:
        raise ValueError("append must be 'before' or 'after'")
    if bad_roots not in [None, "incr", "decr", "auto"]:
        raise ValueError("bad_roots must be None, 'auto', 'incr' or 'decr'")
    if fix_mode not in ["on", "reflect"]:
        raise ValueError("fix_mode must be 'on' or 'reflect'")
    if mirror not in [None, "0", "180"]:
        raise ValueError("mirror must be None, '0' or '180'")
    if method not in ['svd', 'qr', 'cholesky', 'tls']:
        raise ValueError("Invalid method")
    if trace.ndim != 1:
        raise ValueError("trace must be 1D")

    # bad_roots auto mode
    if bad_roots == "auto":
        if mode in ("f", "fb"):
            bad_roots = "incr"
        else:
            bad_roots = "decr"

    x = trace[slice]    # extract region to use for finding LP coefficients

    if mirror is not None:  # make mirror image if selected
        x = make_mirror(x, mirror)

    if mode == "fb":
        a = find_lpc_fb(x, order, bad_roots, fix_mode, method)
        mode = "f"
    elif mode == "bf":
        a = find_lpc_bf(x, order, bad_roots, fix_mode, method)
        mode = "b"
    else:
        # form the LP equation matrix and vector
        D, d = make_Dd(x, order, mode)
        a = find_lpc(D, d, method)  # determine the LP prediction filter

    # stabilize roots if needed
    if bad_roots is not None:           # stabilize roots if needed
        poles = find_roots(a, mode)  # find roots (poles)
        poles = fix_roots(poles, bad_roots, fix_mode)  # fix roots
        # reverse filter when calculated filter is in wrong direction
        if (mode == "b" and append == "after") or (mode == "f" and
                                                   append == "before"):
            poles = [1. / pole for pole in poles]
            mode = {'f': 'b', 'b': 'f'}[mode]
        a = find_coeff(poles, mode)  # find LP filter from roots
    else:
        # reverse filter when calculated filter is in wrong direction
        if (mode == "b" and append == "after") or (mode == "f" and
                                                   append == "before"):
            a = reverse_filter(a, mode)
    # extrapolate the trace using the prediction filter
    ntrace = extrapolate(trace, a, pred, append)
    return ntrace


##################
# LP2D functions #
##################


def lp2d(data, pred, P, M, mirror='0', fix_points=True, method='svd'):
    """
    Perform a forward 2D linear prediction extrapolation on data.

    Use the 2D linear prediction algorithm presented in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 98, 192-199.
    to extend the last (1) axis by `pred` points. A PxM prediction matrix, C,
    is formed by solving the modified linear prediction equation given by:

    data[n,m] = /sigma_{l=0}^{P-1} /sigma_{k=1}^M C_{l,k}*data[n-l,m-k]

    For all valid points in data. This prediction matrix together with the
    data matrix with a mirror image appended is used to extend the last (1)
    axis by pred points resulting in a new array of size [N_0, N_1+pred] where
    N_0 and N_1 are the sizes of the original data.  To linear predict
    both dimensions this function should be used twice with a transpose
    between the uses.

    Backward linear prediction using this method is not possible as the
    method depends on being able to mirror the data before the first collected
    point. In backwards mode this would correspond to being able to correctly
    determine points after the last point which cannot be determined using the
    mirror method. A backward prediction matrix can be calculated but would
    not prove useful.

    The forward-backward averaging of the linear prediction coefficients is
    not possible as there is no characteristic polynomial to root and reflect.
    Therefore the backward prediction matrix cannot be reversed.

    Parameters
    ----------
    data : ndarray
        2D NMR data (time domain for last axes).
    pred : int
        Number of points to predict along the last (1) axes.
    P : int
        Prediction matrix length along the non-predicted (0) axis.
    M : int
        Prediction matrix length along the predicted (1) axis.
    mirror : {'0' or '180'}
        Method to use for forming the mirror image of the non-predicted axis.
        '0' indicated no initial delay, '180' for a half-point delay.
    fix_points : bool
        True to reduce predicted points with magnitude larger than the largest
        data point. False leaved predicted points  unaltered.
    method : {'svd', 'qr', 'cholesky', 'tls'}
        Method used to calculate the LP prediction matrix.  See :py:func:`lp`
        for a description of these methods.

    Returns
    -------
    ndata : ndarray
        2D NMR data with `pred` points appended to the last (1) axes.

    Notes
    -----
    The axes in this function are reversed as compared to the JMR paper.

    """
    # check parameters
    if data.ndim != 2:
        raise ValueError("data must be a 2D array")
    if mirror not in ['0', '180']:
        raise ValueError("mirror must be '0' or '180'")
    if method not in ['svd', 'qr', 'cholesky', 'tls']:
        raise ValueError("method must be 'svd', 'qr', 'cholesky' or 'tls'")

    # form lp2d equation matrix and vector
    D, d = make_lp2d_Dd(data, P, M, 'f')

    # Solve lp2d equation to find the prediction matrix
    c = find_lpc(D, d, method)
    C = c.reshape(P, M)

    # extrapolate the 2D data using the prediction matrix
    return extrapolate_2d(data, C, pred, fix_points, mirror)


def extrapolate_2d(x, C, pred, fix_points, mirror):
    """
    Extrapolate points along the 1st axis using the lp2d algorithm.
    """
    # find the prediction matrix shape and flatten it
    P, M = C.shape
    c = C.flatten()

    # find data parameters
    x_max = x.max()
    N_0, N_1 = x.shape

    # create a empty matrix
    if mirror == "0":
        new = np.empty((2 * N_0 - 1, N_1 + pred), dtype=x.dtype)
        plane = N_0 - 1   # index of first non-mirrored point
    else:
        new = np.empty((2 * N_0, N_1 + pred), dtype=x.dtype)
        plane = N_0     # index of first non-mirrored plane

    last = new.shape[0]     # number of rows in new matrix

    # fill the matrix with the mirrored version of each column
    for i in range(N_1):
        new[:, i] = make_mirror(x[:, i], mirror)

    # fill each new column with predicted values
    # i,j give coordinates of top-left corner of PxM reading matrix
    # after filling a column, replace the whole column with the mirrored column

    for j in range(N_1 - M, N_1 - M + pred):  # column index loop
        for i in range(plane - P + 1, last - P + 1):  # row index loop
            new[i + P - 1, j + M] = np.sum(new[i:i + P, j:j + M].flat * c)

            if fix_points:  # reduce predicted point is needed
                if new[i + P - 1, j + M] > x_max:
                    new[i + P - 1, j + M] = ((x_max * x_max) /
                                             (new[i + P - 1, j + M]))

        # fill the column with the mirrored column so it can be read in the
        # next iteration of the loop
        new[:, j + M] = make_mirror(new[plane:, j + M], mirror)

    return new[plane:]


def make_lp2d_Dd(x, P, M, mode='f'):
    """
    Form the lp2d equation matrix and vector.
    """
    # Form the D and d' matrix and vector in the 2DLP equation:
    # D*c = d'
    # where c is a flattened prediction matrix ordered:
    # (P-1,M) (P-1,M-1) ... (P-1,1) (P-2,M-1) ... (P-2,1) (P-3,M-1) ... ...
    # (2,1) (1,M-1) ... (1,1) (0,M-1) ... (0,1)
    if mode == 'b':
        # backward mode not implemented
        # this would have the d' value as the
        # top left corner
        # this can be done with the same code
        # after reversing x,
        # x = x[::-1,::-1]
        raise NotImplemented

    # Build D and d' row by row by flattening a PxM region of the x matrix
    # starting at 0, 0 and moving down to the bottom of the matrix, then moving
    # back to the top and over one row, and continue till all data has been
    # read. d' is filled with the element next to the right corner of this
    # PxM region.

    N_0, N_1 = x.shape  # length of the matrix

    count_P = N_0 - P + 1   # number of valid starting position vertically
    # number of valid starting position horizontally
    # taking into account the element next to the
    # bottom right corner is the predicted value.
    count_M = N_1 - M

    # create an empty D matrix
    D = np.empty((count_P * count_M, P * M), dtype=x.dtype)
    d = np.empty((count_P * count_M, 1), dtype=x.dtype)

    # fill D and d' row by row (i,j) give coordinates of top-left corner of
    # PxM reading matrix
    for j in range(count_M):
        for i in range(count_P):
            D[j * count_P + i] = x[i:i + P, j:j + M].flat
            d[j * count_P + i] = x[i + P - 1, j + M]
    return D, d


##############################################
# Cadzow/Minimum variance signal enhancement #
##############################################


def cadzow(data, M, K, niter, min_var=False):
    """
    Perform a (row wise) Cadzow-like signal enhancement on 1D or 2D data.

    Performs a Cadzow-like signal enhancement with optional adjustment
    of singular values using the minimum variance method as described in:
    Chen, VanHuffel, Decanniere, VanHecke, JMR, 1994, 109A, 46-55.

    For 2D data performs independent enhancement on each row of data array.

    Parameters
    ----------
    data : ndarray
        1D or 2D NMR data to enhance.
    M : int
        Large prediction order.  For best results should be between
        K + 5 and 2 * K.
    K : int
        Reduced prediction order.
    niter : int
        Number if iteration of the Cadzow procedure to perform.
    min_var : bool
        True to adjust retained singular values using the  minimum variance
        method. False does not correct the singular values and is the
        Cadzow method.

    Returns
    -------
    ndata : ndarray
        Array of enhanced data

    """
    if data.ndim == 1:
        for i in range(niter):
            data = cadzow_single(data, M, K, min_var)
        return data
    elif data.ndim == 2:
        for trace in data:
            for i in range(niter):
                trace = cadzow_single(trace, M, K, min_var)
        return data
    else:
        raise ValueError("data must be a 1D or 2D array")


def cadzow_single(x, M, K, min_var=False):
    """
    Perform a single iteration of Cadzow signal enhancement on a 1D vector

    See :py:func:`cadzow` for documentation.

    """
    # variable names based upon Chen et al, JMR 1994 109A 46
    # form the Hankel data matrix X
    N = len(x)
    L = N - M + 1
    X = scipy.linalg.hankel(x[:L], x[L - 1:])

    # Compute the SVD of X
    U, s, Vh = scipy.linalg.svd(X)

    # correct the singular values and truncate the rank K
    Ul = np.asmatrix(U[:, :K])
    # first K columns of V are first K rows of Vh
    Vlh = np.asmatrix(Vh[:K, :])
    sl = s[:K]

    if min_var:  # adjust singular values using minimum variance method
        # estimate variance of noise singular values
        s2 = (1. / (M - K)) * np.power(s[K:], 2).sum()
        sl = np.array([l - s2 / l for l in sl])

    Sl = np.asmatrix(np.diag(sl))

    # compute enhanced data vector for rank-reduced data matrix
    Xp = Ul * Sl * Vlh

    xp = np.empty_like(x)
    for i, v in enumerate(range(M - 1, -L, -1)):
        # the anti-diagonal is the diagonal with rows reversed
        xp[i] = np.diag(Xp[:, ::-1], v).mean()
    return xp


###################################################
# Linear Prediction parametric modeling functions #
###################################################


def lp_model(trace, slice=slice(None), order=8, mode="f", mirror=None,
             method="svd", full=False):
    """
    Use Linear Prediction to model 1D NMR time domain data.

    Parameters
    ----------
    trace : 1D ndarray
        One dimensional time domain NMR data to model.
    slice : slice object, optional
        Slice object which selects the region along the last axis to use in LP
        equation.  The default, slice(None), will use all points.
    order : int
        Prediction order, number of LP coefficients calculated.
    mode : {'f', 'b'}
        Mode to generate LP filter. 'f' for forward,'b' for backward.
    mirror : {None, '0', '180'}
        Mode to form mirror image of data before processing.  None will
        process the data trace as provided (no mirror image). '0' or '180'
        forms a mirror image of the sliced trace to calculate the LP filter.
        '0' should be used with data with no delay, '180' with data
        with an initial half-point delay.
    method : {'svd', 'qr', 'choleskey', 'tls'}
        Method to use to calculate the LP filter. Choices are a SVD ('svd'), QR
        ('qr'), or Choleskey ('choleskey') decomposition, or Hankel SVD
        ('hsvd').
    full : bool
        True to return amplitudes and phases calculated by performing a least
        squares fitting to the data after LP modeling.  False will return only
        the damping (relaxation) factors and signal frequencies.

    Returns
    -------
    damp : list
        List of damping (relaxation) factors found from LP modeling.
    freq : list
        List of signal frequencies found from LP modeling.
    amp : list, optional
        List of signal amplitudes found by least squares fitting of data after
        LP modeling, only returned when `full` parameter is True.
    phase : list, optional.
        List of signal phases found by least squares fitting of data after LP
        modeling, only returned when `full` parameter is True.


    Notes
    -----
    When backward LP is used the signal roots are reflected before calculating
    model parameters.

    """
    # check for bad arguments
    if mode not in ["f", "b"]:
        raise ValueError("mode must be 'f' or 'b'")
    if method not in ['svd', 'qr', 'cholesky', 'tls', 'hsvd']:
        raise ValueError("Invalid method")
    if trace.ndim != 1:
        raise ValueError("trace must be a 1D array")

    x = trace[slice]    # extract region to use for finding LP coefficients

    if mirror is not None:  # make mirror image if requested
        x = make_mirror(x, mirror)

    # calculate LP coefficient and factor to find poles
    if method in ['svd', 'qr', 'cholseky', 'tls']:
        D, d = make_Dd(x, order, mode)   # form the LP equation elements
        a = find_lpc(D, d, method)      # find LP coefficients
        poles = find_roots(a, mode)     # find roots

    elif method == "hsvd":
        poles = find_lproots_hsvd(x, M=order, K=order, mode=mode, zmethod='sm')
    else:
        raise ValueError("Invalid method")

    # reverse poles if we have backward poles
    if mode == "b":
        poles = [1. / pole for pole in poles]

    # determine the damping factor and frequencies from the roots
    damp = [root2damp(pole) for pole in poles]
    freq = [root2freq(pole) for pole in poles]

    if full is False:
        return damp, freq

    # perform Least Squares fitting to determine amplitudes and phases.

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
    B = np.vstack([poles ** (i) for i in range(len(x))])
    z, resid, rank, s = np.linalg.lstsq(B, np.array(x))

    # Now the z_n = amp_n*exp(phase_n*i), use this to determine the amplitudes
    # and phases
    amp = [cof2amp(cof) for cof in z]
    phase = [cof2phase(cof) for cof in z]

    return damp, freq, amp, phase


###############################################################
# functions to determine signal parameters from LP parameters #
###############################################################


def root2damp(pole):
    """
    Calculate the damping factor from a LP root
    """
    # damping factor is the radius in the complex plane -1/pi*ln(|pole|)
    return -1. / (np.pi) * np.log(np.abs(pole))


def root2freq(pole):
    """
    Calculate the frequency from a LP root
    """
    # frequency is the angle from the x-axis to the point in the complex plane
    # arg(pole) / (2*pi) = atan2(imag,real) / (2*pi)
    return np.arctan2(pole.imag, pole.real) / (2. * np.pi)


def cof2amp(z):
    """
    Calculate a signal amplitude from a model coefficient
    """
    # z = amp*exp(phase*i) so amp is abs(z)
    return np.abs(z)


def cof2phase(z):
    """
    Calculate a signal phase from a model coefficient
    """
    # z = amp*exp(phase(i) so phase is the arg(z)
    return np.arctan2(z.imag, z.real)


##############################
# data preparation functions #
##############################


def make_D(x, order, mode):
    """
    Make the LP equation D matrix (Da = d')
    """
    L = len(x) - order
    if mode == "f":
        return scipy.linalg.hankel(x[:L], x[L - 1:-1])
    elif mode == "b":
        return scipy.linalg.hankel(x[1:L + 1], x[L:])
    else:
        raise ValueError("mode must be 'f' or 'b'")


def make_little_d(x, order, mode):
    """
    Make the LP equation d' vector (Da = d')
    """
    if mode == "f":
        return x[order:].reshape(len(x) - order, 1)
    elif mode == "b":
        L = len(x) - order
        return x[:L].reshape(L, 1)
    else:
        raise ValueError("mode must be 'f' or 'b'")


def make_Dd(x, order, mode):
    """
    make the LP equation D matrix and d' vector (Da=d')
    """
    return make_D(x, order, mode), make_little_d(x, order, mode)


def make_mirror(x, mode):
    """
    Make a mirror image trace.

    Reflects trace over zero as described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1990, 90, 405

    When mode is "0" (no initial delay) form the an array with length 2N-1:
        x_n-1 ... x_1 x_0 x_1 ... x_n-1

    When mode is "180" (half point delay) form an array with length 2N:
        x_n-1 .. x_1 x_0 x_0 x_1 ... x_n-1


    Parameters
    ----------
    x : ndarray
        1D array to form mirrored trace from.
    mode : {'180', '0'}
        Mirror mode, see above.

    """
    if mode == "0":
        return np.concatenate((x[:0:-1], x))
    elif mode == "180":
        return np.concatenate((x[::-1], x))


###########################################################
# LP prediction filter calculation functions (find_lpc_*) #
###########################################################

# the coefficients returned from these functions depend on the mode of the
# prediction.  Forward LP returns coefficients ordered m, m-1, ...1
# Backward LP returns 1, 2, ..., m where m is the order of the prediction.


def find_lpc(D, d, method):
    """
    Find linear prediction filter using a provided method.
    """
    if method == "svd":
        return find_lpc_svd(D, d)
    elif method == "qr":
        return find_lpc_qr(D, d)
    elif method == "cholesky":
        return find_lpc_cholesky(D, d)
    elif method == "tls":
        return find_lpc_tls(D, d)
    else:
        raise ValueError("invalid method")


def find_lpc_svd(D, d):
    """
    Find linear prediction filter using single value decomposition.
    """
    L = D.shape[0]
    m = D.shape[1]
    U, s, Vh = scipy.linalg.svd(D)  # SVD decomposition
    U, Vh = np.asmatrix(U), np.asmatrix(Vh)   # make U and Vh matrices
    Si = pinv_diagsvd(s, m, L)      # construct the pseudo-inverse sigma matrix
    return np.array(Vh.H * Si * U.H * d)


# the next 3 lines and the pinv_diagsvd function were adapted from the
# scipy.linalg.pinv2 function - jjh
eps = np.finfo('float').eps
feps = np.finfo('single').eps
_array_precision = {'f': 0, 'd': 1, 'F': 0, 'D': 1}


def pinv_diagsvd(s, m, L):
    """
    Construct the pseudo-inverse of the sigma matrix from singular values
    """
    t = s.dtype.char
    cond = {0: feps * 1e3, 1: eps * 1e6}[_array_precision[t]]
    cutoff = s[0] * cond
    Si = np.zeros((m, L), t)
    for i in range(len(s)):
        if s[i] > cutoff:
            Si[i, i] = 1.0 / np.conj(s[i])
    return Si


def find_lpc_qr(D, d):
    """
    Find linear prediction filter using QR decomposition.
    """
    q, r = scipy.linalg.qr(D)
    q, r = np.asmatrix(q), np.asmatrix(r)

    # SPEED
    # the next line is slow and the use of pinv2 should be avoided as
    # pseudo inversion of r involves a computationally expensive SVD
    # decomposition which is not needed.  Rather r*x = q.H*d should be
    # solved for x using LAPACK's ZTRTRS function (or similar function with
    # different prefix).  This is not currently available in scipy/numpy and
    # therefore is not used here.
    return scipy.linalg.pinv2(r) * q.H * d


def find_lpc_cholesky(D, d):
    """
    Find linear prediction filter using a Cholesky decomposition.
    """
    # form the normal equation (D.H*D)*a = D.H*d
    # SPEED
    # this can be improved by using the Hankel nature of D
    D = np.asmatrix(D)
    DhD = np.asmatrix(np.dot(D.H, D))
    Dhd = np.asmatrix(np.dot(D.H, d))

    c, lower = scipy.linalg.cho_factor(DhD)     # Compute Cholesky decomp.
    return scipy.linalg.cho_solve((c, lower), Dhd)  # solve normal equation


def find_lpc_tls(D, d):
    """
    Find linear prediction filter using the Total Least Squares method
    """
    m = D.shape[1]  # the order of the prediction
    E = np.append(D, d, axis=1)     # form the augmented data matrix
    U, s, Vh = scipy.linalg.svd(E)  # SVD decomposition of augmented matrix
    V = np.conj(Vh.T)               # Hermetian transpose
    return (-1. / V[m, m] * V[:m, m]).reshape((m, 1))


def find_lpc_fb(x, order, bad_roots, fix_mode, method):
    """
    Determine LP coefficient using forward-backward linear prediction.

    Averages LP coefficients generated from solving the forward and backward
    linear prediction equations after reversing the roots of characteristic
    polynomial of the backward solution.  Method is described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 100, 202-207.

    Description of parameters can be found in :py:func:`lp`.

    """
    # find forward LP coefficients
    D, d = make_Dd(x, order, 'f')
    a = find_lpc(D, d, method)

    # stabilize roots if needed
    if bad_roots is not None:
        poles = find_roots(a, 'f')
        poles = fix_roots(poles, bad_roots, fix_mode)
        a = find_coeff(poles, 'f')
    # store the forward coefficients
    forward_a = a.copy()

    # find the backwards LP coefficients
    D, d = make_Dd(x, order, 'b')
    a = find_lpc(D, d, method)

    # find poles, reverse poles
    poles = find_roots(a, 'b')
    poles = [1. / pole for pole in poles]
    # stabilize roots if needed
    if bad_roots is not None:
        poles = fix_roots(poles, bad_roots, fix_mode)
    # find the backward predicted, forward ordered coefficients
    backward_a = find_coeff(poles, 'f')

    # average the forward and backward coefficients
    return (forward_a + backward_a) / 2.


def find_lpc_bf(x, order, bad_roots, fix_mode, method):
    """
    Determine LP coefficient using backward-forward linear prediction.

    Averages LP coefficients generated from solving the forward and backward
    linear prediction equations after reversing the roots of characteristic
    polynomial of the forward solution.  Similar to method described in:
    G. Zhu and A. Bax, Journal of Magnetic Resonance, 1992, 100, 202-207.

    Description of parameters can be found in :py:func:`lp` function.

    """
    # find backward LP coefficients
    D, d = make_Dd(x, order, 'b')
    a = find_lpc(D, d, method)

    # stabilize roots if needed
    if bad_roots is not None:
        poles = find_roots(a, 'b')
        poles = fix_roots(poles, bad_roots, fix_mode)
        a = find_coeff(poles, 'b')
    # store the forward coefficients
    backward_a = a.copy()

    # find the forward LP coefficients
    D, d = make_Dd(x, order, 'f')
    a = find_lpc(D, d, method)

    # find poles, reverse poles
    poles = find_roots(a, 'f')
    poles = [1. / pole for pole in poles]
    # stabilize roots if needed
    if bad_roots is not None:
        poles = fix_roots(poles, bad_roots, fix_mode)
    # find the forward predicted, backward ordered coefficients
    forward_a = find_coeff(poles, 'b')

    # average the forward and backward coefficients
    return (forward_a + backward_a) / 2.


#####################################
# root finding and fixing functions #
#####################################


def find_lproots_hsvd(x, M, K, mode, zmethod='sm'):
    """
    Find LP roots (poles) using the HSVD method

    Perform a HSVD linear prediction to determine signal roots (poles) as
    described in:
    Barkhuijsen, DeBeer, and Van Ormondt, JMR, 1987, 73, 553

    Parameters x, M and K are the same as those described in the above article.
    zmethod refer to the method used to calculate Z', either a least-squares
    method (lstsq) can be used to solve U_b*Z'=U_t or the Sherman-Morrison
    formula (sm) can be used to avoid the full matrix inversion with equation
    [12] being used to find Z'. The Sherman-Morrison method should be faster
    with similar precision.

    Parameters
    ----------
    x : 1D ndarray
        1D trace of NMR data in the time domain, the FID.
    M : int
        Length (M+1) of data matrix to form.
    K : int
        Reduced prediction order (number of signal roots) Must be less than the
        smaller of M + 1 or len(x) - M.
    mode : {'f', 'b'}
        Mode to perform LP. 'f' for forward,'b' for backward.
    zmethod : {'lstsq', 'sm'}
        Method used to find Z' 'lstsq' for least squares, 'sm' for
        Sherman-Morrison.

    Returns
    -------
    y : ndarray
        Array of signal roots (poles)

    """
    # check parameters
    if mode not in ['f', 'b']:
        raise ValueError("mode must be 'f' or 'b'")
    if zmethod not in ['lstsq', 'sm']:
        raise ValueError("zmethod must be 'lstsq' or 'sm'")
    if K > min(M + 1, len(x) - M):
        raise ValueError("K must be less than min(M + 1, len(x) - M)")

    # form the data matrix X
    N = len(x)
    L = N - M - 1

    if mode == "f":
        X = scipy.linalg.hankel(x[:L + 1], x[L:])
    else:
        # for backward LP we need to make the hankel matrix:
        # x_N-1 x_N-2 ... x_N-M-1
        # x_N-2 x_N-3 ... x_N-M-2
        # ...
        # x_M   x_M-1 ... x_0
        X = scipy.linalg.hankel(x[:M - 1:-1], x[M::-1])

    # SVD of data matrix and truncation of U to form Uk
    U, s, Vh = scipy.linalg.svd(X)
    Uk = np.asmatrix(U[:, :K])   # truncated U matrix of rank K
    Ub = Uk[:-1]            # Uk with bottom row removed
    Ut = Uk[1:]             # Uk with top row removed

    # calculate the Z' matrix
    if zmethod == 'lstsq':  # solve Ub*Z' = Ut using least-squares
        Zp, resid, rank, s = scipy.linalg.lstsq(Ub, Ut)
    else:
        # solve using equation [12]:
        # Z' = (Ek + (u*uh / (1-uh*u)) ) * Ub.H*Ut
        uh = Uk[-1]  # bottom row of Uk
        u = uh.H
        Zp = (np.eye(K, dtype=u.dtype) + (u * uh / (1. - uh * u))) * Ub.H * Ut

    # diagonalization (find eigenvalues) of Z' to yield roots
    return scipy.linalg.eigvals(Zp)


def find_roots(a, mode="f"):
    """
    Find LP roots (poles) from a set of LP coefficients.

    Parameters
    ----------
    a : array
        LP coefficients.
    mode : {'f', 'b'}
        Mode of LP coefficients.  'f' for coefficients ordered m, m - 1,..., 1.
        'b' for coefficients ordered 1, 2, ...., m

    Returns
    -------
    roots : array
        LP roots (poles)

    """
    if mode not in ['f', 'b']:
        raise ValueError("mode must be 'f' or 'b'")

    # STABILITY
    # the algorithm here is that used by numpy roots, build the companion
    # matrix and find its eigenvalues.  These values should be polished for
    # better numerical accuracy.

    # np.roots expects a array, p, with coefficients
    # p[0] * x**n + p[1] * x**(n-1] + ... + p[n-1]*x + p[n]
    # in forward mode LP the coefficients are ordered m,m-1,...1
    # in backward mode LP the coefficient are ordered is 1,2,...,m
    # To calculate the roots, create a leading 1.0+0.0j and reverse if needed.
    p = np.empty(len(a) + 1, dtype=a.dtype)
    p[0] = (1.0 + 0.0j)

    if mode == "f":      # reverse for forward LP
        p[1:] = -a.flat[::-1]
    else:   # backward LP
        p[1:] = -a.flat[:]
    return np.roots(p)


def find_coeff(poles, mode="f"):
    """
    Find LP coefficients from a set of LP roots (poles).

    Parameters
    ----------
    poles : ndarray
        Array of LP roots (poles)
    mode : {'f', 'b'}
        Mode in which LP coefficients should be returned.  'f' for coefficients
        ordered m, m - 1,..., 1. 'b' for coefficients ordered 1, 2, ...., m.

    Returns
    -------
    c : ndarray
        LP coefficients ordered according to `mode`.

    """
    # STABILITY
    # the algorithm used here is numpy poly function which convolves
    # the roots to find the coefficients, the accuracy of this method
    # depends on the dtype of the poles parameter.
    if mode not in ['f', 'b']:
        raise ValueError("mode must be 'f'or 'b'")

    if mode == 'f':  # reverse resulting coefficients
        return np.squeeze(-np.poly(poles)[:0:-1])
    else:   # keep coefficients as is
        return np.squeeze(-np.poly(poles)[1:])


def reverse_filter(a, mode):
    """
    Reverse a filter (change forward LP to backwards LP).
    """
    nmode = {'f': 'b', 'b': 'f'}[mode]
    # find roots, replace each root with 1/root, then recalculate filter
    return find_coeff([1. / pole for pole in find_roots(a, mode)], nmode)


def fix_roots(poles, fix_roots="incr", fix_mode="reflect"):
    """
    Fix (stabilize) LP roots.

    Parameters
    ----------
    poles : ndarray
        Array of LP roots (poles).
    fix_roots : {'incr', 'decr'}
        Type of roots which to consider bad and to stabilize.  Either those
        with increasing signals 'incr' or decreasing signals 'decr'.
    fix_mode : {'on', 'reflect'}
        Method used to stabilize bad roots, 'on' to move the roots onto the
        unit circle, 'reflect' to reflect bad roots across the unit circle.

    Returns
    -------
    npoles : ndarray
        Array of stabilized LP roots (poles).

    """
    if fix_roots not in ["incr", "decr"]:
        raise ValueError("fix_roots must be 'incr' or 'decr'")
    if fix_mode not in ["on", "reflect"]:
        raise ValueError("fix_mode must be 'on' or 'reflect'")

    if fix_roots == "incr":     # remove increasing signals
        for i, pole in enumerate(poles):
            if np.abs(pole) > 1:
                # print("Fixing root:",i)
                if fix_mode == "on":
                    poles[i] = pole / np.abs(pole)
                else:
                    poles[i] = 1 / np.conj(pole)
    else:   # remove decreasing signals
        for i, pole in enumerate(poles):
            if np.abs(pole) < 1:
                # print("Fixing root:",i)
                if fix_mode == "on":
                    poles[i] = pole / np.abs(pole)
                else:
                    poles[i] = 1 / np.conj(pole)
    return poles


###########################
# Extrapolation functions #
###########################


def extrapolate(trace, a, pred, append):
    """
    Extrapolate points using LP prediction filter.

    Parameters
    ----------
    trace : 1D ndarray
        1D array to extrapolate from and append to.
    a : ndarray
        LP coefficients, must be ordered according to direction of
        extrapolation.
    pred : int
        Number of points to predict using LP.
    append : {'a', 'b'}
        Location to append new points, 'a' for after the current data, 'b' for
        before the current data.

    Returns
    -------
    ntrace : 1D ndarray
        1D array with extrapolated points appended

    """
    m = len(a)      # LP order
    M = len(trace)   # number of points in original trace
    ntrace = np.empty((M + pred), dtype=trace.dtype)

    if append not in ["after", "before"]:
        raise ValueError("append must be 'a' or 'b'")

    if append == "after":   # append after trace
        ntrace[:M] = trace
        for i in range(pred):
            ntrace[M + i] = np.sum(np.multiply(ntrace[M - m + i:M + i],
                                   a.flat))
        return ntrace

    if append == "before":   # append before trace
        ntrace[-M:] = trace
        for i in range(pred):
            ntrace[pred - i - 1] = np.sum(np.multiply(
                ntrace[pred - i:pred + m - i], a.flat))
        return ntrace
