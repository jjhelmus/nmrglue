""" Unit tests for nmrglue/fileio/glue.py module """

import numpy as np
from numpy.testing import assert_array_equal

from nmrglue.fileio import glue


class FakeDataset:
    """Minimal stand-in for an h5py Dataset: shape, dtype and __getitem__."""

    def __init__(self, array):
        self._array = np.asarray(array)
        self.shape = self._array.shape
        self.dtype = self._array.dtype

    def __getitem__(self, key):
        return self._array[key]


def test_glue_2d_contract():
    """glue_2d follows the data_nd numeric axis-order contract."""
    data = np.arange(12).reshape(3, 4)
    d = glue.glue_2d(FakeDataset(data))

    assert d.shape == (3, 4)
    assert d.ndim == 2
    assert_array_equal(d[:], data)
    assert_array_equal(np.array(list(d)), data)


def test_glue_2d_transpose_and_swapaxes():
    """glue_2d transpose and swapaxes return glue_2d with reordered axes."""
    data = np.arange(12).reshape(3, 4)
    d = glue.glue_2d(FakeDataset(data))

    transposed = d.transpose()
    assert isinstance(transposed, glue.glue_2d)
    assert transposed.shape == (4, 3)
    assert_array_equal(transposed[:], data.T)

    swapped = d.swapaxes(0, 1)
    assert isinstance(swapped, glue.glue_2d)
    assert swapped.shape == (4, 3)
    assert_array_equal(swapped[:], data.T)


def test_glue_3d_contract():
    """glue_3d follows the data_nd numeric axis-order contract."""
    data = np.arange(24).reshape(2, 3, 4)
    d = glue.glue_3d(FakeDataset(data))

    assert d.shape == (2, 3, 4)
    assert d.ndim == 3
    assert_array_equal(d[:], data)

    transposed = d.transpose()
    assert isinstance(transposed, glue.glue_3d)
    assert transposed.shape == (4, 3, 2)
    assert_array_equal(transposed[:], data.transpose(2, 1, 0))
