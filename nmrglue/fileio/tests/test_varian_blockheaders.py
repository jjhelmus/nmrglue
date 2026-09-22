"""Tests for reading Varian blocks containing multiple traces."""

import io

import numpy as np
from numpy.testing import assert_array_equal

from nmrglue.fileio import varian


def test_get_block_ntraces_with_blockheader():
    """Read multi-trace data when block-header details are requested."""
    blockheader = bytes(28)
    trace_data = np.array([1, 2, 3, 4], dtype=">i4").tobytes()

    dic, data = varian.get_block_ntraces(
        io.BytesIO(blockheader + trace_data),
        ntraces=2,
        pts=2,
        nbheaders=1,
        dt=np.dtype(">i4"),
        read_blockhead=True,
    )

    assert dic["scale"] == 0
    assert_array_equal(data, [[1, 2], [3, 4]])
