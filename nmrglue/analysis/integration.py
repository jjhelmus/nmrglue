import numpy as np


def integrate(data, unit_conv, limits, unit='ppm', noise_limits=None,
              norm_to_range=None, calibrate=1.0):
    r"""
    Integrate one 1D data array within limits given in units. Data points must
    be equally spaced.

    Functional form of integration is:

    .. math::
        value = \sum_a^b s(x_{i}) dx

    Where:
    s is the signal, a and b are the limits of integration and dx is the width
    of each bin.

    The integration error due to baseline noise is calculated as:

    .. math::
        error = \sigma_{vol} = \sigma \sqrt{n}

    if the noise_limits are set.

    Where:

    .. math::
        n =  \frac{|b-a|}{dx}+1

    sigma is the standard deviation of the baseline noise. n is the number
    of bins in the integration range.

    Parameters
    ----------
    data: array like
        1d array of intensity
    unit_conv: `fileiobase.unit_conversion` object
        unit_conversion object associated with data
    limits: array like
        With shape (2,) or (P, 2). Array with lower and upper integration
        limit. Or array with P rows of lower and upper integration limits.
    noise_limits: Optional[array like]
        With shape(2, ). Array with lower and upper limits to section of data
        with only noise. A larger range will likely yield a more accurate
        estimate. It is unwise to use the very end of the spectrum for the
        noise range.
    norm_to_range: Optional[int]
        If given, all values are normalized to the value determined within the
        range given by limits[norm_to_range, :]
    calibrate: Optional[float]
        If norm_to_range is given, the values are re-normalized so that the
        integral within limits[norm_to_range, :] equals some number.

    Returns
    -------
    array
        [value, ...] integration values

    if noise_limits is given:

    array
        [[value, error], ...] where error a one sigma estimate of the error
        only from the spectrum noise

    """
    # convert form limits in units to lower and upper index
    limits = np.array(limits)
    if limits.size == 2:
        limits = np.expand_dims(limits, axis=0)
    inds = [(unit_conv(x0, unit), unit_conv(x1, unit)) for (x0, x1) in limits]
    inds = np.array([sorted(ind) for ind in inds])

    # sum part of the integral
    sum_slice = np.array([np.sum(data[slice(*ind)]) for ind in inds])

    # dx part of the integral
    scale = unit_conv.ppm_scale()
    dx = abs(scale[1]-scale[0])

    # put together the integral
    values = sum_slice * dx

    if noise_limits is not None:

        # determine the standard deviation of the noise
        noise_inds = (unit_conv(noise_limits[0], unit),
                      unit_conv(noise_limits[1], unit))
        noise_inds = sorted(noise_inds)

        # the error (from noise) is  std * dx * limits range
        std = np.std(data[slice(*noise_inds)])
        errors = std * np.sqrt(inds[:, 1] - inds[:, 0])

        if norm_to_range is not None:
            # You must normalize the error before you normalize the value!
            errors = (errors / values[norm_to_range]) * calibrate
            values = (values / values[norm_to_range]) * calibrate
        return np.vstack((values, errors)).T

    else:

        if norm_to_range is not None:
            values = (values / values[norm_to_range]) * calibrate
        return values


def ndintegrate(data, unit_conv, limits, unit='ppm', noise_limits=None):
    r"""
    Integrate one nD data array within limits given in units. Data points must
    be equally spaced. Can only integrate one region per function call.

    The integration error due to baseline noise is calculated as:

    .. math::
        error = \sigma_{vol} = \sigma \sqrt{\prod_i^{d} n_{i}},

    if the noise_limits are set.

    Where:
    sigma is the standard deviation of the baseline noise. n is the number
    of bins in the integration range for each d dimensions.

    See integrate for more information.

    Parameters
    ----------
    data: array like
        1d array of intensity
    unit_convs: [`fileiobase.unit_conversion`, ] list
        list of unit_conversion object associated with each dim of data.
    limits: array like
        With shape (2,) or (d, 2). 1D Array with lower and upper integration
        limits for 1D . Or array with d rows of lower and upper integration
        limits for each dimension.
    noise_limits: Optional[array like]
        With shape(2, ). Array with lower and upper limits to section of data
        with only noise. A larger range will likely yield a more accurate
        estimate. It is unwise to use the very end of the spectrum for the
        noise range.

    Returns
    -------
    array
        [value, ...] integration values

    if noise_limits is given:

    array
        [[value, error], ...] where error a one sigma estimate of the error
        only from the spectrum noise
    """

    # determine the dimensionality of the data.
    d = np.ndim(data)

    try:
        iter(unit_conv)
    except TypeError:
        unit_conv = [unit_conv, ]

    if d != len(unit_conv):
        mesg = 'A unit_conversion object is needed for each dimension.'
        raise ValueError(mesg)

    limits = np.array(limits)
    if limits.ndim == 1:
        limits = np.expand_dims(limits, axis=0)

    if limits.shape[0] != d and limits.shape[1] != 2:
        mesg = 'A lower and upper limit is needed for each dimension.'
        raise ValueError(mesg)

    inds = [(uc(x[0], unit), uc(x[1], unit))
            for (uc, x) in zip(unit_conv, limits)]
    inds = [sorted(x) for x in inds]

    # the integrate_nd needs to be scaled by the bin width in ppm
    ppm_scales = [x.ppm_scale() for x in unit_conv]
    dx = np.prod(np.array([abs(x[1]-x[0]) for x in ppm_scales]))

    slice_sum = (data[tuple([slice(x[0], x[1]) for x in np.flipud(inds)])]).sum()

    value = slice_sum * dx

    if noise_limits is None:
        return value

    else:
        noise_limits = np.array(noise_limits)
        if noise_limits.size == 2:
            noise_limits = np.expand_dims(noise_limits, axis=0)

        if noise_limits.shape[0] != d and noise_limits.shape[1] != 2:
            mesg = 'If given, a noise limit is needed for each dimension.'
            raise ValueError(mesg)

        noise_inds = [(uc(x[0], unit), uc(x[1], unit))
                      for (uc, x) in zip(unit_conv, noise_limits)]
        noise_inds = [sorted(x) for x in noise_inds]

        # see docstring of integrate
        nm = np.prod(np.array([abs(x[1]-x[0]) for x in noise_inds]))
        std = np.std(data[[slice(x[0], x[1])for x in np.flipud(noise_inds)]])

        error = std * np.sqrt(nm)
        return np.hstack((value, error))
