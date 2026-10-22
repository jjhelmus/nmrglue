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


@pytest.mark.parametrize('ndim', [1, 2, 3])
@pytest.mark.parametrize('unit', ['ppm', 'hz'])
@pytest.mark.parametrize('reverse', [False, True])
def test_ndintegrate_noise_uses_signal_bins_and_bin_volume(
        ndim, unit, reverse):
    """Propagate baseline noise through the same weighted sum as the volume."""
    shape = (6, 8, 10)[-ndim:]
    data = np.arange(np.prod(shape), dtype=float).reshape(shape) % 7 + 1
    widths = np.arange(1, ndim + 1) * 0.125
    ucs = [ng.fileiobase.unit_conversion(n, False, n * 100 * dx, 100, 0)
           for n, dx in zip(shape, widths)]
    signal_slices = tuple(slice(2, min(n - 1, 6)) for n in shape)
    noise_slices = tuple(slice(0, 2) for n in shape)

    def limits_for(slices):
        limits = [[uc.unit(sl.start, unit), uc.unit(sl.stop, unit)]
                  for uc, sl in zip(ucs, slices)][::-1]
        return [x[::-1] for x in limits] if reverse else limits

    signal = data[signal_slices]
    noise = data[noise_slices]
    expected = [signal.sum() * np.prod(widths),
                noise.std() * np.sqrt(signal.size) * np.prod(widths)]
    actual = ndintegrate(data, ucs[::-1], limits_for(signal_slices),
                         unit=unit, noise_limits=limits_for(noise_slices))

    np.testing.assert_allclose(actual, expected, rtol=1e-12)
    if ndim == 1:
        flat_noise_limits = limits_for(noise_slices)[0]
        actual = ndintegrate(data, ucs[0], limits_for(signal_slices),
                             unit=unit, noise_limits=flat_noise_limits)
        np.testing.assert_allclose(actual, expected, rtol=1e-12)
