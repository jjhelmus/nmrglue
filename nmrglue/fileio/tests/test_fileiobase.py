""" Unit tests for nmrglue/fileio/fileiobase.py module """

import copy
import os

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng
import pytest


# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
NMRPIPE_1D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_1d_freq.fid')


class DummyDataND(ng.fileiobase.data_nd):
    """Minimal data_nd implementation for testing base-class operations."""

    def __init__(self, order=None, fshape=(2, 3)):
        self.fshape = fshape
        if order is None:
            order = range(len(fshape))
        self.order = tuple(order)
        self.dtype = np.dtype("float64")
        self.__setdimandshape__()

    def __fcopy__(self, order):
        return DummyDataND(order, self.fshape)

    def __fgetitem__(self, slices):
        return np.arange(np.prod(self.fshape)).reshape(self.fshape)[slices]


def test_uc_from_freqscale():
    """
    Test that `uc_from_freqscale` gives equivalent results as `uc_from_udic`.
    """
    from nmrglue.fileio.fileiobase import uc_from_freqscale

    # read frequency test data
    dic, data = ng.pipe.read(NMRPIPE_1D_FREQ)

    # make udic and uc using uc_from_udic
    udic = ng.pipe.guess_udic(dic, data)
    uc = ng.fileiobase.uc_from_udic(udic)

    ppm_scale = uc.ppm_scale()
    uc_from_ppm = uc_from_freqscale(ppm_scale, udic[0]['obs'], 'ppm')
    new_ppm_scale = uc_from_ppm.ppm_scale()
    assert_array_equal(ppm_scale, new_ppm_scale)

    hz_scale = uc.hz_scale()
    uc_from_hz = uc_from_freqscale(hz_scale, udic[0]['obs'], 'hz')
    new_hz_scale = uc_from_hz.hz_scale()
    assert_array_equal(hz_scale, new_hz_scale)

    khz_scale = hz_scale * 1.0e3
    uc_from_khz = uc_from_freqscale(khz_scale, udic[0]['obs'], 'khz')
    new_khz_scale = uc_from_khz.hz_scale() * 1.0e3
    assert_array_equal(khz_scale, new_khz_scale)


# regression test for https://github.com/jjhelmus/nmrglue/issues/113
def test_uc_with_float_size():
    uc = ng.fileiobase.unit_conversion(
        size=64.0, cplx=False, sw=1.0, obs=1.0, car=1.0)
    scale = uc.ppm_scale()
    assert len(scale) == 64


def test_update_uc():
    uc = ng.fileiobase.unit_conversion(
        size=64.0, cplx=False, sw=1.0, obs=1.0, car=1.0)
    uc2 = ng.fileiobase.update_uc(uc, size=10, cplx=True, sw=2.0, car=5.2, obs=3.0)
    assert uc2._size == 10
    assert uc2._cplx is True
    assert abs(uc2._sw - 2.0) < 1e-5
    assert abs(uc2._car - 5.2) < 1e-5
    assert abs(uc2._obs - 3.0) < 1e-5


def test_data_nd_copy():
    data = DummyDataND()

    copied = copy.copy(data)

    assert copied is not data
    assert copied.order == data.order
    assert_array_equal(copied[:], data[:])


def test_data_nd_copy_and_transform_orders_are_independent():
    data = DummyDataND(fshape=(2, 3, 4))
    expected = np.arange(24).reshape(data.fshape)

    transformed = data.transpose(-1, -2, -3)
    copied_transformed = copy.copy(transformed)
    transformed_copy = copy.copy(data).swapaxes(-3, 2)

    assert data.order == (0, 1, 2)
    assert transformed.order == (2, 1, 0)
    assert copied_transformed.order == transformed.order
    assert transformed_copy.order == (2, 1, 0)
    assert_array_equal(copied_transformed[:], expected.transpose(2, 1, 0))
    assert_array_equal(transformed_copy[:], expected.swapaxes(-3, 2))


@pytest.mark.parametrize("axes", [(-1, 0), (0, -1)])
def test_data_nd_swapaxes_with_negative_axis(axes):
    data = DummyDataND()

    swapped = data.swapaxes(*axes)

    assert swapped.order == (1, 0)
    assert swapped.shape == (3, 2)
    assert_array_equal(swapped[:], np.arange(6).reshape(2, 3).swapaxes(0, 1))


def test_data_nd_swapaxes_3d_with_negative_axis():
    data = DummyDataND(fshape=(2, 3, 4))

    swapped = data.swapaxes(-3, 2)

    assert swapped.order == (2, 1, 0)
    expected = np.arange(24).reshape(2, 3, 4).swapaxes(-3, 2)
    assert_array_equal(swapped[:], expected)


@pytest.mark.parametrize("axes", [(-3, 0), (0, -3), (2, 0), (0, 2)])
def test_data_nd_swapaxes_rejects_invalid_axis(axes):
    with pytest.raises(ValueError):
        DummyDataND().swapaxes(*axes)


@pytest.mark.parametrize(
    "axes, expected_order",
    [
        ((0, 1, 2), (0, 1, 2)),
        ((2, 1, 0), (2, 1, 0)),
        ((-1, -2, -3), (2, 1, 0)),
        ((-3, 1, -1), (0, 1, 2)),
    ],
)
def test_data_nd_transpose_axes(axes, expected_order):
    data = DummyDataND(fshape=(2, 3, 4))

    transposed = data.transpose(*axes)

    assert transposed.order == expected_order
    expected = np.arange(24).reshape(data.fshape).transpose(axes)
    assert_array_equal(transposed[:], expected)


@pytest.mark.parametrize(
    "axes",
    [(-4, 0, 1), (0, 0, 1), (0, 1), (0, 1, 2, 3)],
)
def test_data_nd_transpose_rejects_invalid_axes(axes):
    with pytest.raises(ValueError):
        DummyDataND(fshape=(2, 3, 4)).transpose(*axes)
