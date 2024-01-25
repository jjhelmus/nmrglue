""" Unit tests for nmrglue/analysis/integration.py module """

import numpy as np
from scipy.stats import multivariate_normal
from numpy.testing import assert_array_almost_equal

import nmrglue as ng
from nmrglue.analysis.integration import integrate, ndintegrate


# test helper functions
def _build_1d_ppm_scale():
    n = 500
    sw = 1/1000e-6
    obs = 100  # Mhz

    dx = (sw / n)/obs
    ppm_scale = np.arange(dx, dx * n, dx)

    # To match the apparent NMR convention the array has to be reversed.
    return ppm_scale[::-1]


# test functions
def test_1d_integrate():
    """ Test integration of synthetic 1D data with two peaks."""
    # generate test scale
    ppm_scale = _build_1d_ppm_scale()
    uc = ng.fileio.fileiobase.uc_from_freqscale(ppm_scale, 100)

    # generate test data
    data = (multivariate_normal.pdf(ppm_scale, mean=5, cov=0.01) +
            multivariate_normal.pdf(ppm_scale, mean=8, cov=0.01) * 2.)

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(ppm_scale, data)
    # plt.show()

    # Test with a single integral region
    assert integrate(data, uc, (4, 6)) - 1.0 <= np.finfo(float).eps
    assert integrate(data, uc, (7, 9)) - 2.0 <= np.finfo(float).eps

    # Test with multiple integral regions:
    limits = [(4, 6), (7, 9)]
    values_1 = integrate(data, uc, limits)
    assert_array_almost_equal(values_1,  [1., 2.])

    # Test with multiple integral regions and re-normalization of the values
    # to the value of the range between 7 ppm and 9 ppm
    values_2 = integrate(data, uc, limits, norm_to_range=1)
    assert_array_almost_equal(values_2,  [0.5, 1.])

    # Test with multiple integral regions and re-normalization of the values
    # to the value of the range between 7 ppm and 9 ppm and re-calibration
    # of the the value to 3.
    values_2 = integrate(data, uc, limits, norm_to_range=1, calibrate=3)
    assert_array_almost_equal(values_2,  [0.5 * 3, 3.])


def test_1d_integrate_with_noise():
    """ Test integration of synthetic 1D data with two peaks and noise"""

    # seed random number
    np.random.seed(0)

    # generate test scale
    ppm_scale = _build_1d_ppm_scale()
    uc = ng.fileio.fileiobase.uc_from_freqscale(ppm_scale, 100)

    # generate test data
    data = (multivariate_normal.pdf(ppm_scale, mean=5, cov=0.01) +
            multivariate_normal.pdf(ppm_scale, mean=8, cov=0.01) * 2)
    noise = np.random.randn(*data.shape) * 0.01
    data = data + noise

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(ppm_scale, data)
    # plt.show()

    # estimate of how much the noise can throw off the the integral
    limits = (4.0, 6.0)
    n = abs(uc(4.5, 'ppm')-uc(5.5, 'ppm'))
    max_error = np.std(noise) * np.sqrt(n)

    # Test with a single integral region sanity check.
    assert abs(integrate(data, uc, limits) - 1) <= max_error

    # Test renormalization of norms
    results = integrate(data, uc, ((4, 6), (7, 9)), noise_limits=(1, 2),
                        norm_to_range=1)

    # Test renormalization of values.
    assert abs(results[0, 0] - 0.5) <= max_error


def test_1d_ndintegrate():
    """ Test integration of synthetic with ndintegrate.
    """
    # generate test scale
    ppm_scale = _build_1d_ppm_scale()
    uc = ng.fileio.fileiobase.uc_from_freqscale(ppm_scale, 100)

    # generate test data
    data = (multivariate_normal.pdf(ppm_scale, mean=5, cov=0.01) +
            multivariate_normal.pdf(ppm_scale, mean=8, cov=0.01) * 2.)

    # Test with a single  integral region
    assert abs(ndintegrate(data, uc, (4, 6)) - 1.0) <= 1e-10
    assert abs(ndintegrate(data, uc, [(7, 9), ]) - 2.0) <= 1e-10


def test_2d_ndintegrate():
    """ Test integration of synthetic 2D data with ndintegrate. """
    # generate test scale
    ppm_scale = _build_1d_ppm_scale()
    uc = ng.fileio.fileiobase.uc_from_freqscale(ppm_scale, 100)

    # generate test data
    x, y = np.meshgrid(ppm_scale, ppm_scale)
    pos = np.empty(x.shape + (2,))
    pos[:, :, 0] = x
    pos[:, :, 1] = y
    rv = multivariate_normal([5, 8], [0.01, 0.01])
    data = rv.pdf(pos)

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.contour(x, y, data)
    # plt.show()

    limits = ((4, 6), (7, 9))

    assert abs(ndintegrate(data, [uc, uc], limits) - 1.0) <= 1e-10
