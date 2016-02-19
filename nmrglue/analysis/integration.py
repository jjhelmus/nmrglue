import numpy as np


def integrate(data, unit_conv, limits, unit='ppm', noise_limits=None,
              norm_to_range=None, calibrate=1.0):
    """
    Integrate one 1D data array within limits given in units. Data points must
    be equally spaced.

    Functional form of integration is:

    .. math::
        value = \sum_a^b s(x_{i}) dx

    Where:
    s is the signal, a and b are the limits of integration and dx is the width
    of each bin.

    A simple error analysis is optionally performed as:

    ..math::
        error = \sigma_{vol} = \sigma \sqrt{n}

    Where:
    .. math::
        n = \frac{|b-a|}{dx}+1

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
