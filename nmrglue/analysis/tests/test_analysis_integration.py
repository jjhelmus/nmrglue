""" Unit tests for nmrglue/analysis/integration.py module """

import numpy as np
import pytest
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


@pytest.mark.parametrize('width', [0.25, 2.0])
@pytest.mark.parametrize('norm_to_range', [None, 1])
@pytest.mark.parametrize('unit', ['ppm', 'hz'])
def test_integrate_noise_has_integral_units(width, norm_to_range, unit):
    """The integral and its noise estimate must use the same bin width."""
    data = np.arange(16, dtype=float) + 1
    uc = ng.fileiobase.unit_conversion(16, False, 16 * 100 * width, 100, 0)
    limits = [[uc.unit(a, unit), uc.unit(b, unit)]
              for a, b in [(4, 8), (8, 14)]]
    noise_limits = [uc.unit(0, unit), uc.unit(4, unit)]
    expected_values = np.array([data[4:8].sum(), data[8:14].sum()]) * width
    expected_errors = data[:4].std() * np.sqrt([4, 6]) * width
    if norm_to_range is not None:
        factor = 3.0 / expected_values[norm_to_range]
        expected_values *= factor
        expected_errors *= factor

    actual = integrate(data, uc, limits, unit=unit, noise_limits=noise_limits,
                       norm_to_range=norm_to_range, calibrate=3.0)

    np.testing.assert_allclose(actual[:, 0], expected_values, rtol=1e-12)
    np.testing.assert_allclose(actual[:, 1], expected_errors, rtol=1e-12)


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
