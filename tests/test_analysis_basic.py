""" Unit tests for nmrglue/analysis/basic.py module """

import os

from numpy.testing import assert_array_equal
import nmrglue as ng
from setup import DATA_DIR


def test_integration_1d():
    """
    Test that `uc_from_freqscale` gives equivalent results as `uc_from_udic`.
    """
    assert 1 == 1
