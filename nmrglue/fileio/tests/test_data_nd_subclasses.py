"""Test ``data_nd`` operations on concrete file-format subclasses."""

import copy
import os

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from nmrglue.fileio import bruker, fileiobase, pipe, rnmrtk


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
NMRPIPE_2D_FREQ = os.path.join(DATA_DIR, "nmrpipe_2d_freq.ft2")


def assert_view(view, expected, expected_order):
    """Check a lazy view's metadata, full values, and representative slices."""
    assert isinstance(view, fileiobase.data_nd)
    assert not isinstance(view, np.ndarray)
    assert tuple(view.order) == tuple(expected_order)
    assert view.shape == expected.shape
    assert view.dtype == expected.dtype
    # data_nd historically squeezes length-one axes from indexed results.
    assert_array_equal(view[:], np.squeeze(expected))
    assert_array_equal(view[0], np.squeeze(expected[0]))
    assert_array_equal(view[..., 0], np.squeeze(expected[..., 0]))


def assert_data_nd_contract(data, expected):
    """Check copies and axis transforms against a materialized NumPy oracle."""
    original_order = tuple(data.order)
    original_shape = data.shape

    copied = copy.copy(data)
    assert type(copied) is type(data)
    assert copied is not data
    assert_view(copied, expected, original_order)

    transposed = data.transpose(-1, -2)
    assert type(transposed) is type(data)
    assert_view(transposed, expected.transpose(-1, -2), original_order[::-1])

    swapped = data.swapaxes(-1, 0)
    assert type(swapped) is type(data)
    assert_view(swapped, expected.swapaxes(-1, 0), original_order[::-1])

    copied_transformed = copy.copy(swapped)
    assert type(copied_transformed) is type(data)
    assert_view(
        copied_transformed,
        expected.swapaxes(-1, 0),
        original_order[::-1],
    )

    transformed_copy = copy.copy(data).swapaxes(-1, 0)
    assert type(transformed_copy) is type(data)
    assert_view(
        transformed_copy,
        expected.swapaxes(-1, 0),
        original_order[::-1],
    )
    assert tuple(data.order) == original_order
    assert data.shape == original_shape


def test_nmrpipe_data_nd_contract():
    _, data = pipe.read_lowmem(NMRPIPE_2D_FREQ)
    _, expected = pipe.read(NMRPIPE_2D_FREQ)

    assert_data_nd_contract(data, expected)


@pytest.mark.parametrize(
    "isfloat, storage_dtype",
    [(False, "<i4"), (True, "<f8")],
)
def test_bruker_data_nd_contract(tmp_path, isfloat, storage_dtype):
    expected = np.arange(6, dtype=np.int32).reshape(1, 6)
    filename = tmp_path / "bruker.bin"
    filename.write_bytes(expected.astype(storage_dtype).tobytes())
    data = bruker.bruker_nd(
        str(filename), expected.shape, False, False, isfloat=isfloat
    )

    assert_data_nd_contract(data, expected)
    assert copy.copy(data).isfloat is isfloat
    assert data.swapaxes(-1, 0).isfloat is isfloat


def test_rnmrtk_data_nd_contract(tmp_path):
    expected = np.arange(6, dtype=np.float32).reshape(2, 3)
    filename = tmp_path / "rnmrtk.sec"
    filename.write_bytes(expected.astype("<f4").tobytes())
    data = rnmrtk.rnmrtk_nd(
        str(filename), expected.shape, False, False
    )

    assert_data_nd_contract(data, expected)
