import os

import numpy as np
from numpy.testing import assert_array_almost_equal

import nmrglue as ng

from setup import DATA_DIR


def test_tecmag_load_time_domain():
    ref1, data = ng.tecmag.read(os.path.join(
        DATA_DIR, 'tecmag', 'LiCl_ref1.tnt'))

    real, imag, usec = np.loadtxt(os.path.join(
        DATA_DIR, 'tecmag', 'LiCl_ref1.txt'), skiprows=3, unpack=True)

    assert_array_almost_equal(data.real.squeeze(), real, decimal=3)
    assert_array_almost_equal(data.imag.squeeze(), imag, decimal=3)
    assert_array_almost_equal(np.arange(ref1['npts'][0]) *
                              ref1['dwell'][0] * 1e6, usec, decimal=3)
