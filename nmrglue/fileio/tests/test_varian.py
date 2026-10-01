"""Tests for the Varian file reader and writer."""

import copy
import struct

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from nmrglue.fileio import varian


FILE_HEADER = struct.Struct(">6lhhl")
BLOCK_HEADER_SIZE = struct.calcsize(">4hl4f")
HEADER_NAMES = (
    "nblocks",
    "ntraces",
    "np",
    "ebytes",
    "tbytes",
    "bbytes",
    "vers_id",
    "status",
    "nbheaders",
)


def make_file_dic(shape, *, np_points=None, nblocks=None, ntraces=1,
                  ebytes=4, tbytes=None, bbytes=None):
    """Return a minimal int32 Varian file dictionary."""
    if np_points is None:
        np_points = shape[-1] * 2
    if nblocks is None:
        nblocks = int(np.prod(shape[:-1])) if len(shape) > 1 else 1
    if tbytes is None:
        tbytes = np_points * ebytes
    if bbytes is None:
        bbytes = BLOCK_HEADER_SIZE + ntraces * tbytes
    return {
        "np": np_points,
        "nblocks": nblocks,
        "ntraces": ntraces,
        "ebytes": ebytes,
        "tbytes": tbytes,
        "bbytes": bbytes,
        "vers_id": 0,
        "status": 5,
        "nbheaders": 1,
        "S_FLOAT": 0,
        "S_32": 1,
    }


def make_data(shape):
    """Return deterministic complex data representable as int32."""
    values = np.arange(np.prod(shape)).reshape(shape)
    return (values + 1j * (values + 1)).astype(np.complex128)


def read_raw_file_header(path):
    """Read the nine fields of a Varian file header independently."""
    with path.open("rb") as fileobj:
        values = FILE_HEADER.unpack(fileobj.read(FILE_HEADER.size))
    return dict(zip(HEADER_NAMES, values))


@pytest.mark.parametrize("shape", [(4,), (2, 3), (2, 3, 4)])
def test_write_fid_lowmem_corrects_structural_header(tmp_path, shape):
    """Correct all header fields derived from the physical data layout."""
    data = make_data(shape)
    dic = make_file_dic(
        shape,
        np_points=2,
        nblocks=1,
        tbytes=8,
        bbytes=36,
    )
    path = tmp_path / "test.fid"

    with pytest.warns(UserWarning):
        varian.write_fid_lowmem(path, dic, data)

    expected_nblocks = int(np.prod(shape[:-1])) if len(shape) > 1 else 1
    expected = {
        "np": shape[-1] * 2,
        "nblocks": expected_nblocks,
        "ntraces": 1,
        "ebytes": 4,
        "tbytes": shape[-1] * 2 * 4,
        "bbytes": BLOCK_HEADER_SIZE + shape[-1] * 2 * 4,
    }
    header = read_raw_file_header(path)
    assert {key: header[key] for key in expected} == expected
    assert {key: dic[key] for key in expected} == expected
    assert path.stat().st_size == FILE_HEADER.size + (
        expected["nblocks"] * expected["bbytes"]
    )

    _, normal = varian.read_fid(path, shape=shape, torder="f")
    assert normal.dtype == np.dtype("complex128")
    assert_array_equal(normal, data)

    if len(shape) == 1:
        # The historical low-memory reader represents an unspecified 1D
        # shape as a single-block 2D object.
        _, lowmem = varian.read_fid_lowmem(path, torder="f")
        lowmem_shape = (1, shape[0])
        first_expected = last_expected = data
    else:
        _, lowmem = varian.read_fid_lowmem(path, shape=shape, torder="f")
        lowmem_shape = shape
        first_expected = data[0]
        last_expected = data[-1]
    assert lowmem.shape == lowmem_shape
    assert lowmem.dtype == np.dtype("complex128")
    assert_array_equal(lowmem[:], data)
    assert_array_equal(lowmem[0], first_expected)
    assert_array_equal(lowmem[-1], last_expected)


@pytest.mark.parametrize("shape", [(4,), (2, 3), (2, 3, 4)])
def test_write_fid_lowmem_matches_regular_writer(tmp_path, shape):
    """Consistent inputs produce the same bytes through both writers."""
    data = make_data(shape)
    dic = make_file_dic(shape)
    regular_path = tmp_path / "regular.fid"
    lowmem_path = tmp_path / "lowmem.fid"

    varian.write_fid(regular_path, copy.deepcopy(dic), data, torder="f")
    varian.write_fid_lowmem(
        lowmem_path,
        copy.deepcopy(dic),
        data,
        torder="f",
    )

    assert lowmem_path.read_bytes() == regular_path.read_bytes()


def test_write_fid_lowmem_corrects_trace_and_element_sizes(tmp_path):
    """Header element and trace sizes follow the bytes actually written."""
    shape = (2, 3)
    data = make_data(shape)
    dic = make_file_dic(shape, ntraces=2, ebytes=2)
    path = tmp_path / "test.fid"

    with pytest.warns(UserWarning):
        varian.write_fid_lowmem(path, dic, data)

    header = read_raw_file_header(path)
    assert header["ntraces"] == 1
    assert header["ebytes"] == 4
    assert header["tbytes"] == 24
    assert header["bbytes"] == 52
    assert path.stat().st_size == FILE_HEADER.size + 2 * 52
    _, lowmem = varian.read_fid_lowmem(path, shape=shape)
    assert_array_equal(lowmem[:], data)


def test_write_fid_lowmem_can_preserve_header_dimensions(tmp_path):
    """Passing correct=False preserves all supplied structural fields."""
    shape = (2, 3)
    dic = make_file_dic(
        shape,
        np_points=2,
        nblocks=1,
        ntraces=2,
        ebytes=2,
        tbytes=4,
        bbytes=36,
    )
    original = copy.deepcopy(dic)
    path = tmp_path / "test.fid"

    with pytest.warns(UserWarning):
        varian.write_fid_lowmem(
            path,
            dic,
            make_data(shape),
            correct=False,
        )

    assert dic == original
    header = read_raw_file_header(path)
    for key in HEADER_NAMES:
        assert header[key] == original[key]
