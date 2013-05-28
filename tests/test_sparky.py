""" Tests for the fileio.sparky submodule """

import tempfile
import os
import os.path

from numpy.testing import assert_array_equal
import nmrglue as ng

from setup import DATA_DIR

# subroutines
def dic_similar(dic1, dic2):
    """ Check is Sparky parameter dictionaries are the same. """
    if dic1.keys() != dic2.keys():
        print "Not same keys!"
    assert dic1.keys() == dic2.keys()
    for key in dic1.keys():
        if dic1[key] != dic2[key]:
            print key, dic1[key], dic2[key]
        assert dic1[key] == dic2[key]
    return True

def write_readback(dic, data):
    """ Write and readback of Sparky data. """
    # write out and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write(tf, dic, data)
    rdic, rdata = ng.sparky.read(tf)
    assert_array_equal(data, rdata)
    assert dic_similar(dic, rdic)
    os.remove(tf)

def lowmem_write_readback(dic, data):
    """ Lowmemory write and readback of Sparky data."""
    # write out and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write_lowmem(tf, dic, data)
    rdic, rdata = ng.sparky.read_lowmem(tf)
    tup = tuple(range(data.ndim))
    assert_array_equal(data[tup], rdata[tup])
    assert dic_similar(dic, rdic)
    os.remove(tf)

# tests
def test_2d():
    """ reading/writing of 2D sparky file """
    dic, data = ng.sparky.read(os.path.join(DATA_DIR, "sparky_2d", "data.ucsf"))
    assert data.shape == (2048, 4096)
    assert round(data[0, 1], 2) == 1601.83
    assert round(data[15, 20], 2) == 4281.06
    write_readback(dic, data)

def test_2d_lowmem():
    """ lowmemory reading/writing of 2D sparky file """
    dic, data = ng.sparky.read_lowmem(os.path.join(DATA_DIR, "sparky_2d", "data.ucsf")
    assert data.shape == (2048, 4096)
    assert round(data[0, 1], 2) == 1601.83
    assert round(data[15, 20], 2) == 4281.06
    lowmem_write_readback(dic, data)

def test_3d():
    """ reading/writing of 3D sparky file """
    dic, data = ng.sparky.read(os.path.join(DATA_DIR, "sparky_3d", "data.ucsf")
    assert data.shape == (128, 128, 4096)
    assert round(data[0, 1, 2], 2) == 25980.13
    assert round(data[11, 15, 20], 2) == -15256.97
    write_readback(dic, data)

def test_3d_lowmem():
    """ lowmemory reading/writing of 3D sparky file """
    dic, data = ng.sparky.read_lowmem(os.path.join(DATA_DIR, "sparky_3d", "data.ucsf")
    assert data.shape == (128, 128, 4096)
    assert round(data[0, 1, 2], 2) == 25980.13
    assert round(data[11, 15, 20], 2) == -15256.97
    lowmem_write_readback(dic, data)
