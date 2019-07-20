"""
Functions for fitting and simulating arbitrary dimensional lineshapes commonly
found in NMR experiments
"""

from __future__ import print_function

import numpy as np

from .leastsqbound import leastsqbound
from .analysisbase import squish
from .lineshapes1d import ls_str2class
from ..fileio import table

pi = np.pi


# table packing/unpacking
def add_to_table(rec, columns, column_names):
    """
    Add (append) multiple columns to a records array.

    Parameters
    ----------
    rec : recarray
        Records array (table).
    columns : list of ndarrays
        List of columns data to append to table.
    column_names : list of str
        List of names of columns.

    Returns
    -------
    nrec : recarray
        Records array with columns added

    """
    for col, col_name in zip(columns, column_names):
        rec = table.append_column(rec, col, name=col_name)
    return rec


def pack_table(pbest, abest, iers, rec, param_columns, amp_column,
               ier_column=None):
    """
    Pack fitting parameters into table

    Parameters
    ----------
    pbest : list
        List of best-fit parameters.  See :py:func:`fit_NDregion` for format.
    abest : list
        List of best-fit amplitudes.
    iers : list
        List of fitting error return values.
    rec : recarray
        Records array (table) to save fitting parameters into. Updated with
        fitting parameter in place.
    param_columns : list
        List of parameter columns in rec. Format is the same as pbest.
    amp_columns : str
        Name of amplitude column in rec.
    ier_column : str or None, optional
        Name of column in rec to save iers to. None will not record this in the
        table.

    """
    # pack the amplitudes
    rec[amp_column] = abest

    # pack the parameters
    for dbest, dcolumns in zip(zip(*pbest), param_columns):
        for p, c in zip(zip(*dbest), dcolumns):
            rec[c] = p

    # pack the iers
    if ier_column is not None:
        rec[ier_column] = iers


def unpack_table(rec, param_columns, amp_column):
    """
    Unpack initial fitting parameters from a table.

    Parameters
    ----------
    rec : recarray
        Records array (table) holding parameters.
    param_columns : list
        List of column names which hold lineshape parameters.  See
        :py:func:`fit_NDregion` for format.
    amp_column : str
        Name of columns in rec holding initial amplitudes.

    Returns
    -------
    params : list
        List of initial parameter in the format required for
        :py:func:`fit_NDregion`.
    amps : list
        List of initial peak amplitudes.

    """
    params = zip(*[zip(*[rec[c] for c in dc]) for dc in param_columns])
    amps = rec[amp_column]
    return params, amps


def estimate_scales(spectrum, centers, box_width, scale_axis=0):
    """
    Estimate scale parameter for peaks in a spectrum.

    Parameters
    ----------
    spectrum : array_like
        NMR spectral data. ndarray or emulated type which can be sliced.
    centers : list
        List of N-tuples indicating peak centers.
    box_width : tuple
        N-tuple indicating box width to add and subtract from peak centers to
        form region around peak to fit.
    scale_axis : int
        Axis number to estimate scale parameters for.

    Returns
    -------
    scales : list
        List of estimated scale parameters.

    """
    shape = spectrum.shape
    bcenters = np.round(np.array(centers)).astype('int')
    scales = []
    # loop over the box centers
    for bc in bcenters:

        # calculate box limits
        bmin = [max(c - w, 0) for c, w in zip(bc, box_width)]
        bmax = [min(c + w + 1, s) for c, w, s in zip(bc, box_width, shape)]

        # cut the spectrum and squish
        s = tuple([slice(mn, mx) for mn, mx in zip(bmin, bmax)])
        scale = squish(spectrum[s], scale_axis)
        scale = scale / scale[0]
        scales.append(scale[1:])
    return scales


# User facing fit/simulation functions
def fit_spectrum(spectrum, lineshapes, params, amps, bounds, ampbounds,
                 centers, rIDs, box_width, error_flag, verb=True, **kw):
    """
    Fit a NMR spectrum by regions which contain one or more peaks.

    Parameters
    ----------

    spectrum : array_like
        NMR data. ndarray or emulated type, must be slicable.
    lineshape :list
        List of lineshapes by label (str) or a lineshape class. See
        :py:func:`fit_NDregion` for details.
    params : list
        P-length list (P is the number of peaks in region) of N-length lists
        of tuples where each each tuple is the optimiztion starting parameters
        for a given peak and dimension lineshape.
    amps : list
        P-length list of amplitudes.
    bounds : list
        List of bounds for parameter of same shape as params.  If none of the
        parameters in a given dimension have limits None can be used,
        otherwise each dimension should have a list or tuple of (min,max) or
        None for each parameter. min or max may be None when there is no
        bounds in a given direction.
    ampbounds : list
        P-length list of bounds for the amplitude with format similar to
        bounds.
    centers : list
        List of N-tuples indicating peak centers.
    rIDs : list
        P-length list of region numbers.  Peak with the same region number
        are fit together.
    box_width : tuple
        Tuple of length N indicating box width to add and subtract from peak
        centers to form regions around peak to fit.
    error_flag : bool
        True to estimate errors for each lineshape parameter and amplitude.
    verb : bool, optional
        True to print a summary of each region fit, False (the default)
        supresses all printing.
    **kw : optional
        Additional keywords passed to the scipy.optimize.leastsq function.

    Returns
    -------
    params_best : list
        Optimal values for lineshape parameters with same format as params
        input parameter.
    amp_best : list
        List of optimal peak amplitudes.
    param_err : list, only returned when error_flag is True
        Estimated lineshape parameter errors with same format as params.
    amp_err : list, only returned when error_flag is True
        Estimated peak amplitude errors.
    iers : list
        List of interger flag from scipy.optimize.leastsq indicating if the
        solution was found for a given peak.  1,2,3,4 indicates that a
        solution was found. Other indicate an error.

    """
    pbest = [[]] * len(params)
    pbest_err = [[]] * len(params)
    abest = [[]] * len(params)
    abest_err = [[]] * len(params)
    iers = [[]] * len(params)
    shape = spectrum.shape

    ls_classes = []
    for l in lineshapes:
        if isinstance(l, str):
            ls_classes.append(ls_str2class(l))
        else:
            ls_classes.append(l)

    cIDs = set(rIDs)    # region values to loop over

    for cID in cIDs:
        cpeaks = [i for i, v in enumerate(rIDs) if v == cID]

        # select the parameter
        cparams = [params[i] for i in cpeaks]
        camps = [amps[i] for i in cpeaks]
        cbounds = [bounds[i] for i in cpeaks]
        campbounds = [ampbounds[i] for i in cpeaks]
        ccenters = [centers[i] for i in cpeaks]

        # find the box edges
        bcenters = np.round(np.array(ccenters)).astype('int')
        bmin = bcenters - box_width
        bmax = bcenters + box_width + 1

        # correct for spectrum edges
        for i in range(len(shape)):
            bmin[:, i][np.where(bmin[:, i] < 0)] = 0
        for i, v in enumerate(shape):
            bmax[:, i][np.where(bmax[:, i] > v)] = v

        # find the region limits
        rmin = edge = np.array(bmin).min(0)
        rmax = np.array(bmax).max(0)

        # cut the spectrum
        s = tuple([slice(mn, mx) for mn, mx in zip(rmin, rmax)])
        region = spectrum[s]

        # add edge to the box limits
        ebmin = bmin - edge
        ebmax = bmax - edge

        # create the weight mask array
        wmask = np.zeros(region.shape, dtype='bool')
        for bmn, bmx in zip(ebmin, ebmax):
            s = tuple([slice(mn, mx) for mn, mx in zip(bmn, bmx)])
            wmask[s] = True

        # add edges to the initial parameters
        ecparams = [[ls.add_edge(p, (mn, mx)) for ls, mn, mx, p in
                    zip(ls_classes, rmin, rmax, g)] for g in cparams]

        # TODO make this better...
        ecbounds = [[list(zip(*[ls.add_edge(b, (mn, mx)) for b in zip(*db)]))
                    for ls, mn, mx, db in zip(ls_classes, rmin, rmax, pb)]
                    for pb in cbounds]

        # fit the region
        t = fit_NDregion(region, ls_classes, ecparams, camps, ecbounds,
                         campbounds, wmask, error_flag, **kw)
        if error_flag:
            ecpbest, acbest, ecpbest_err, acbest_err, ier = t
            cpbest_err = [[ls.remove_edge(p, (mn, mx)) for ls, mn, mx, p in
                          zip(ls_classes, rmin, rmax, g)] for g in ecpbest_err]
        else:
            ecpbest, acbest, ier = t

        # remove edges from best fit parameters
        cpbest = [[ls.remove_edge(p, (mn, mx)) for ls, mn, mx, p in
                   zip(ls_classes, rmin, rmax, g)] for g in ecpbest]

        if verb:
            print("-----------------------")
            print("cID:", cID, "ier:", ier, "Peaks fit", cpeaks)
            print("fit parameters:", cpbest)
            print("fit amplitudes", acbest)

        for i, pb, ab in zip(cpeaks, cpbest, acbest):
            pbest[i] = pb
            abest[i] = ab
            iers[i] = ier

        if error_flag:
            for i, pb, ab in zip(cpeaks, cpbest_err, acbest_err):
                pbest_err[i] = pb
                abest_err[i] = ab

    if error_flag is False:
        return pbest, abest, iers

    return pbest, abest, pbest_err, abest_err, iers


def fit_NDregion(region, lineshapes, params, amps, bounds=None,
                 ampbounds=None, wmask=None, error_flag=False, **kw):
    """
    Fit a N-dimensional region.

    Parameters
    ----------

    region : ndarray
        Region of a NMR data to fit.
    lineshape :list
        List of lineshapes by label (str) or a lineshape class. See
        Notes for details.
    params : list
        P-length list (P is the number of peaks in region) of N-length lists
        of tuples where each each tuple is the optimiztion starting parameters
        for a given peak and dimension lineshape.
    amps : list
        P-length list of amplitudes.
    bounds : list
        List of bounds for parameter of same shape as params.  If none of the
        parameters in a given dimension have limits None can be used,
        otherwise each dimension should have a list or tuple of (min,max) or
        None for each parameter. min or max may be None when there is no
        bounds in a given direction.
    ampbounds : list
        P-length list of bounds for the amplitude with format similar to
        bounds.
    wmask : ndarray, optional
        Array with same shape as region which is used to weight points in the
        error calculation, typically a boolean array is used to exclude
        certain points in the region.  Default of None will include all
        points in the region equally in the error calculation.
    centers : list
        List of N-tuples indicating peak centers.
    error_flag : bool
        True to estimate errors for each lineshape parameter and amplitude.
    **kw : optional
        Additional keywords passed to the scipy.optimize.leastsq function.

    Returns
    -------
    params_best : list
        Optimal values for lineshape parameters with same format as params
        input parameter.
    amp_best : list
        List of optimal peak amplitudes.
    param_err : list, only returned when error_flag is True
        Estimated lineshape parameter errors with same format as params.
    amp_err : list, only returned when error_flag is True
        Estimated peak amplitude errors.
    iers : list
        List of interger flag from scipy.optimize.leastsq indicating if the
        solution was found for a given peak.  1,2,3,4 indicates that a
        solution was found. Other indicate an error.


    Notes
    -----

    The lineshape parameter:

    Elements of the lineshape parameter list can be string indicating the
    lineshape of given dimension or an instance of a lineshape class
    which provide a sim method which takes two arguments, the first being the
    length of the lineshape the second being a list of lineshape parameters,
    and returns a simulated lineshape as well as a nparam method which when
    given the length of lineshape returns the number of parameters needed to
    describe the lineshape. Currently the following strings are allowed:

    * 'g' or 'gauss'    Gaussian (normal) lineshape.
    * 'l' or 'lorentz'  Lorentzian lineshape.
    * 'v' or 'voigt'    Voigt lineshape.
    * 'pv' or 'pvoight' Pseudo Voigt lineshape
    * 's' or 'scale'    Scaled lineshape.

    The first four lineshapes (Gaussian, Lorentzian, Voigt and Pseudo Voigt)
    all take a FWHM scale parameter.

    The following are all valid lineshapes parameters for a 2D Gaussian peak:

    * ['g','g']
    * ['gauss','gauss']
    * [ng.lineshapes1d.gauss(),ng.lineshapes1d.gauss()]

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
        if isinstance(l, str):
            ls_classes.append(ls_str2class(l))
        else:
            ls_classes.append(l)

    # determind the number of parameter in each dimension
    dim_nparam = [c.nparam(l) for l, c in zip(shape, ls_classes)]

    # parse params
    n_peaks = len(params)
    p0 = []
    for i, guess in enumerate(params):  # peak loop
        if len(guess) != ndim:
            err = "Incorrect number of params for peak %i"
            raise ValueError(err % (i))

        for j, dim_guess in enumerate(guess):    # dimension loop
            if len(dim_guess) != dim_nparam[j]:
                err = "Incorrect number of parameters in peak %i dimension %i"
                raise ValueError(err % (i, j))

            for g in dim_guess:  # parameter loop
                p0.append(g)

    # parse the bounds parameter
    if bounds is None:   # No bounds
        peak_bounds = [[(None, None)] * i for i in dim_nparam]
        bounds = [peak_bounds] * n_peaks

    if len(bounds) != n_peaks:
        raise ValueError("Incorrect number of parameter bounds provided")

    # build the parameter bound list to be passed to f_NDregion
    p_bounds = []
    for i, peak_bounds in enumerate(bounds):  # peak loop

        if peak_bounds is None:
            peak_bounds = [[(None, None)] * i for i in dim_nparam]

        if len(peak_bounds) != ndim:
            err = "Incorrect number of bounds for peak %i"
            raise ValueError(err % (i))

        for j, dim_bounds in enumerate(peak_bounds):    # dimension loop

            if dim_bounds is None:
                dim_bounds = [(None, None)] * dim_nparam[j]

            if len(dim_bounds) != dim_nparam[j]:
                err = "Incorrect number of bounds for peak %i dimension %i"
                raise ValueError(err % (i, j))

            for k, b in enumerate(dim_bounds):    # parameter loop
                if b is None:
                    b = (None, None)

                if len(b) != 2:
                    err = "No min/max for peak %i dim %i parameter %i"
                    raise ValueError(err % (i, j, k))

                p_bounds.append(b)

    # parse amps parameter
    if len(amps) != n_peaks:
        raise ValueError("Incorrect number of amplitude guesses provided")
    p0 = list(amps) + p0  # amplitudes appended to front of p0

    # parse ampbounds parameter
    if ampbounds is None:
        ampbounds = [(None, None)] * n_peaks

    if len(ampbounds) != n_peaks:
        raise ValueError("Incorrect number of amplitude bounds")

    to_add = []
    for k, b in enumerate(ampbounds):
        if b is None:
            b = (None, None)

        if len(b) != 2:
            err = "No min/max for amplitude bound %i"
            raise ValueError(err % (k))
        to_add.append(b)
    p_bounds = to_add + p_bounds    # amplitude bound at front of p_bounds

    # parse the wmask parameter
    if wmask is None:   # default is to include all points in region
        wmask = np.ones(shape, dtype='bool')
    if wmask.shape != shape:
        err = "wmask has incorrect shape:" + str(wmask.shape) +   \
              " should be " + str(shape)
        raise ValueError(err)

    # DEBUGGING
    # print("--------------------------------")
    # print(region)
    # print(ls_classes)
    # print(p0)
    # print(p_bounds)
    # print(n_peaks)
    # print(dim_nparam)
    # print("=================================")
    # for i,j in zip(p0,p_bounds):
    #     print(i, j)

    # include full_output=True when errors requested
    if error_flag:
        kw["full_output"] = True

    # perform fitting
    r = f_NDregion(region, ls_classes, p0, p_bounds, n_peaks, wmask, **kw)

    # DEBUGGING
    # print(r)

    # unpack results depending of if full output requested
    if "full_output" in kw and kw["full_output"]:
        p_best, cov_xi, infodic, mesg, ier = r
    else:
        p_best, ier = r

    # unpack and repack p_best
    # pull off the ampltides
    amp_best = p_best[:n_peaks]

    # split the remaining parameters into n_peaks equal sized lists
    p_list = split_list(list(p_best[n_peaks:]), n_peaks)

    # for each peak repack the flat parameter lists to reference by dimension
    param_best = [make_slist(l, dim_nparam) for l in p_list]

    # return as is if no errors requested
    if error_flag is False:
        return param_best, amp_best, ier

    # calculate errors
    p_err = calc_errors(region, ls_classes, p_best, cov_xi, n_peaks, wmask)

    # unpack and repack the error p_err
    # pull off the amplitude errors
    amp_err = p_err[:n_peaks]

    # split the remaining errors into n_peaks equal sized lists
    pe_list = split_list(list(p_err[n_peaks:]), n_peaks)

    # for each peak repack the flat errors list to reference by dimension
    param_err = [make_slist(l, dim_nparam) for l in pe_list]

    return param_best, amp_best, param_err, amp_err, ier


def sim_NDregion(shape, lineshapes, params, amps):
    """
    Simulate an N-dimensional region with one or more peaks.

    Parameters
    ----------
    shape : tuple of ints
        Shape of region.
    lineshapes : list
        List of lineshapes by label (str) or a lineshape class. See
        :py:func:`fit_NDregion` for additional documentation.
    params : list
        P-length list (P is the number of peaks in region) of N-length lists
        of tuples where each each tuple is lineshape parameters for a given
        peak and dimension.
    amps : list
        P-length of peak amplitudes.

    Returns
    -------
    sim : ndarray with shape, shape.
        Simulated region.

    """
    # parse the user-friendly input into a format digestable by s_NDregion
    # parse the shape
    ndim = len(shape)

    # parse the lineshape parameters
    if len(lineshapes) != ndim:
        raise ValueError("Incorrect number of lineshapes provided")

    ls_classes = []
    for l in lineshapes:
        if isinstance(l, str):
            ls_classes.append(ls_str2class(l))
        else:
            ls_classes.append(l)

    # determind the number of parameters in each dimension.
    dim_nparam = [c.nparam(l) for l, c in zip(shape, ls_classes)]

    # parse the params parameter
    n_peaks = len(params)
    p = []
    for i, param in enumerate(params):
        if len(param) != ndim:
            err = "Incorrect number of parameters for peak %i"
            raise ValueError(err % (i))
        for j, dim_param in enumerate(param):
            if len(dim_param) != dim_nparam[j]:
                err = "Incorrect number of parameters in peak %i dimension %i"
                raise ValueError(err % (i, j))

            for g in dim_param:
                p.append(g)

    # parse the amps parameter
    if len(amps) != n_peaks:
        raise ValueError("Incorrect number of amplitudes provided")
    p = list(amps) + p   # amplitudes appended to front of p

    # DEBUGGING
    # print("p",p)
    # print("shape",shape)
    # print("ls_classes",ls_classes)
    # print("n_peaks",n_peaks)

    return s_NDregion(p, shape, ls_classes, n_peaks)


def make_slist(l, t_sizes):
    """
    Create a list of tuples of given sizes from a list

    Parameters
    ----------
    l : list or ndarray
        List or array to pack into shaped list.
    t_sizes : list of ints
        List of tuple sizes.

    Returns
    -------
    slist : list of tuples
        List of tuples of lengths given by t_sizes.

    """
    out = []  # output
    start = 0
    for s in t_sizes:
        out.append(l[start:start + s])
        start = start + s
    return out


def split_list(l, N):
    """ Split list l into N sublists of equal size """
    step = int(len(l) / N)
    div_points = range(0, len(l) + 1, step)
    return [l[div_points[i]:div_points[i + 1]] for i in range(N)]


def calc_errors(region, ls_classes, p, cov, n_peaks, wmask):
    """
    Calcuate the parameter errors from the standard errors of the estimate.

    Parameters
    ----------
    region : ndarray
        Region which was fit.
    ls_classes : list
        List of lineshape classes.
    p : ndarray
        Fit parameters.
    cov : ndarray
        Covariance matrix from least squares fitting.
    n_peaks : int
        Number of peaks in the region.

    Returns
    -------
    errors : ndarray
        Array of standard errors of parameters in p.

    """
    # calculate the residuals
    resid = err_NDregion(p, region, region.shape, ls_classes, n_peaks, wmask)
    SS_err = np.power(resid, 2).sum()   # Sum of squared residuals
    n = region.size  # size of sample XXX not sure if this always makes sense
    k = p.size - 1   # free parameters
    st_err = np.sqrt(SS_err / (n - k - 1))    # standard error of estimate
    if cov is None:   # indicate that parameter errors cannot be calculated.
        return [None] * len(p)
    return st_err * np.sqrt(np.diag(cov))


# internal functions
def s_NDregion(p, shape, ls_classes, n_peaks):
    """
    Simulate an N-dimensional region with one or more peaks.

    Parameters
    ----------
    p : list
        List of parameters, must be a list, modified by function.
    shape : tuple of ints
        Shape of region.
    ls_classes : list
        List of lineshape classes.
    n_peaks : int
        Number of peaks in region.

    Returns
    -------
    r : ndarray
        Simulated region.

    """
    # split the parameter list into a list of amplitudes and peak param lists
    As = [p.pop(0) for i in range(n_peaks)]
    ps = split_list(p, n_peaks)

    # simulate the first region
    A, curr_p = As.pop(0), ps.pop(0)
    r = s_single_NDregion([A] + curr_p, shape, ls_classes)

    # simulate any additional regions
    for A, curr_p in zip(As, ps):
        r = r + s_single_NDregion([A] + curr_p, shape, ls_classes)

    return r


def s_single_NDregion(p, shape, ls_classes):
    """
    Simulate an N-dimensional region with a single peak.

    This function is called repeatly by s_NDregion to build up a full
    simulated region.

    Parameters
    ----------
    p : list
        List of parameters, must be a list.
    shape : tuple
        Shape of region.
    ls_classes : list
        List of lineshape classes.

    Returns
    -------
    r : ndarray
        Simulated region.

    """
    A = p.pop(0)    # amplitude is ALWAYS the first parameter
    r = np.array(A, dtype='float')

    for length, ls_class in zip(shape, ls_classes):
        # print("Making lineshape of", ls_class.name, "with length:", length)
        s_p = [p.pop(0) for i in range(ls_class.nparam(length))]
        ls = ls_class.sim(length, s_p)
        # print("Lineshape is:", ls)
        r = np.kron(r, ls)   # vector direct product flattened
    return r.reshape(shape)


def err_NDregion(p, region, shape, ls_classes, n_peaks, wmask):
    """
    Error function for an N-dimensional region, called by :py:func:`f_NDregion`
    """
    sim_region = s_NDregion(list(p), shape, ls_classes, n_peaks)
    return ((region - sim_region) * wmask).flatten()


def f_NDregion(region, ls_classes, p0, p_bounds, n_peaks, wmask, **kw):
    """
    Fit an N-dimensional regions containing one or more peaks.

    Region is fit using a contrained Levenberg-Marquard optmization algorithm.
    See :py:func:`fit_NDregion` for additional documentation.

    Parameters
    ----------
    region : ndarray
        Region to fit.
    ls_classes : list
        List of lineshape classes.
    p0 : ndarray
        Initial parameters.
    p_bounds : list of tuples
        List of (min, max) bounds for each element of p0.
    n_peaks : int
        Number of peaks in the simulated region.
    wmask : ndarray
        Array with same shape as region which is used to weight points in the
        error calculation, typically a boolean array is used to exclude
        certain points in the region.
    **kw : optional
        Additional keywords passed to the scipy.optimize.leastsq function.

    See Also
    --------
    fit_NDregion : Fit N-dimensional region with user friendly parameter.

    """
    args = (region, region.shape, ls_classes, n_peaks, wmask)
    p_best = leastsqbound(err_NDregion, p0, bounds=p_bounds, args=args, **kw)
    return p_best
