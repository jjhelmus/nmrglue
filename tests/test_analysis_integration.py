""" Unit tests for nmrglue/analysis/integration.py module """

import numpy as np
from scipy.stats import multivariate_normal
from numpy.testing import assert_array_almost_equal

import nmrglue as ng
from nmrglue.analysis.integration import integrate


# test helper functions
def _build_1d_ppm_scale():
    n = 500
    sw = 1/1000e-6
    obs = 100  # Mhz

    dx = (sw /n)/obs
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


def test_1d_integrate_withnoise():
    """ Test integration of synthetic 1D data with two peaks. and noise"""
    # generate test scale
    ppm_scale = _build_1d_ppm_scale()
    uc = ng.fileio.fileiobase.uc_from_freqscale(ppm_scale, 100)

    # generate test data
    data = (multivariate_normal.pdf(ppm_scale, mean=5, cov=0.01) +
            multivariate_normal.pdf(ppm_scale, mean=8, cov=0.01) * 2)
    noise = np.random.randn(*data.shape) * 0.1
    data = data + noise

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(ppm_scale, data)
    # plt.show()

    # estimate of how much the noise can throw off the the integral
    limits = (4.0, 6.0)
    dx = abs(ppm_scale[1]-ppm_scale[0])
    max_error = np.std(data) * dx * abs(limits[0]-limits[1]) * 3

    # Test with a single integral region sanity check.
    assert abs(integrate(data, uc, limits) - 1) <= max_error

    # Test renormalization of norms
    resutls= integrate(data, uc, ((4, 6), (7, 9)), noise_limits=(1, 2),
                       norm_to_range=1)

    # Test renormalization of values.
    assert abs(resutls[0, 0] - 0.5) <= max_error