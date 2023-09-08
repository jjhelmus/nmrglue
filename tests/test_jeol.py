""" Tests for the fileio.jeol submodule """


import tempfile
import os
import shutil
from pathlib import Path

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng

from setup import DATA_DIR


def test_parse():
    pass