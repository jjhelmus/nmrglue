""" Tests for the fileio.varian submodule """

import tempfile
import os
import shutil

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng
from nose.plugins.attrib import attr

from setup import DATA_DIR


# subroutines
def write_readback(dic, data):
    """ write out and readback a Agilent/varian file directory. """
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, dic, data)
    rdic, rdata = ng.varian.read(td)
    shutil.rmtree(td)
    assert_array_equal(data, rdata)
    assert dic == rdic


def lowmem_write_readback(dic, data):
    """ lowmemory write out and readback a Agilent/Varian directory. """
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td, dic, data)
    rdic, rdata = ng.varian.read_lowmem(td)
    # check value [0,1,...]
    s = tuple(range(data.ndim))
    assert data[s] == rdata[s]
    assert dic == rdic
    shutil.rmtree(td)


def write_fid_readback(dic, data, shape, torder):
    """ Writeout and readback a Agilent/Varian fid file """
    tf = tempfile.mktemp(dir=".")
    ng.varian.write_fid(tf, dic, data, torder=torder)
    rdic, rdata = ng.varian.read_fid(tf, shape=shape, torder=torder)
    os.remove(tf)
    assert_array_equal(data, rdata)
    assert dic == rdic


def lowmem_fid_write_readback(dic, data, shape, torder):
    """ Lowmemory writeout and readback a Agilent/Varian fid file """
    tf = tempfile.mktemp(dir=".")
    ng.varian.write_fid_lowmem(tf, dic, data, torder=torder)
    rdic, rdata = ng.varian.read_fid_lowmem(tf, shape=shape, torder=torder)
    # check value [0,1,...]
    s = tuple(range(data.ndim))
    assert data[s] == rdata[s]
    assert dic == rdic
    os.remove(tf)


# tests


@attr(speed='fast')
def test_1d():
    """ reading/writing of 1D Varian file """
    dic, data = ng.varian.read(os.path.join(DATA_DIR, "agilent_1d"))
    assert data.shape == (1500, )
    assert np.abs(data[0].real - 91899.24) <= 0.01
    assert np.abs(data[0].imag - -1964.70) <= 0.01
    assert np.abs(data[1].real - 168844.25) <= 0.01
    assert np.abs(data[1].imag - -49503.41) <= 0.01
    write_readback(dic, data)


@attr(speed='fast')
def test_2d():
    """ reading/writing of 2D Varian file """
    dic, data = ng.varian.read(os.path.join(DATA_DIR, "agilent_2d"))
    assert data.shape == (332, 1500)
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - -223.20) <= 0.01
    assert np.abs(data[10, 18].real - 17.93) <= 0.01
    assert np.abs(data[10, 18].imag - -67.20) <= 0.01
    write_readback(dic, data)


@attr(speed='fast')
def test_2d_lowmem():
    """ low memory reading/writing of 2D Varian file """
    dic, data = ng.varian.read_lowmem(os.path.join(DATA_DIR, "agilent_2d"))
    assert data.shape == (332, 1500)
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - -223.20) <= 0.01
    assert np.abs(data[10, 18].real - 17.93) <= 0.01
    assert np.abs(data[10, 18].imag - -67.20) <= 0.01
    lowmem_write_readback(dic, data)


@attr(speed='fast')
def test_2d_tppi():
    """ reading/writing of 2D Varian file with TPPI encoding """
    dic, data = ng.varian.read(os.path.join(DATA_DIR, "agilent_2d_tppi"))
    assert data.shape == (600, 1400)
    assert np.abs(data[0, 1].real - -4589.29) <= 0.01
    assert np.abs(data[0, 1].imag - -1691.82) <= 0.01
    assert np.abs(data[10, 18].real - -166.62) <= 0.01
    assert np.abs(data[10, 18].imag - -594.73) <= 0.01
    write_readback(dic, data)


@attr(speed='fast')
def test_2d_tppi_lowmem():
    """ low memory reading/writing of 2D Varian file with TPPI encoding """
    dic, data = ng.varian.read_lowmem(os.path.join(DATA_DIR,
                                      "agilent_2d_tppi"))
    assert data.shape == (600, 1400)
    assert np.abs(data[0, 1].real - -4589.29) <= 0.01
    assert np.abs(data[0, 1].imag - -1691.82) <= 0.01
    assert np.abs(data[10, 18].real - -166.62) <= 0.01
    assert np.abs(data[10, 18].imag - -594.73) <= 0.01
    lowmem_write_readback(dic, data)


@attr(speed='slow')
def test_3d():
    """ reading/writing of 3D Varian file """
    dic, data = ng.varian.read(os.path.join(DATA_DIR, "agilent_3d"))
    assert data.shape == (128, 88, 1250)
    assert np.abs(data[0, 1, 2].real - 7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 11, 18].real - -9.36) <= 0.01
    assert np.abs(data[10, 11, 18].imag - -7.75) <= 0.01
    write_readback(dic, data)


@attr(speed='slow')
def test_3d_lowmem():
    """ low memory reading/writing of 3D Varian file """
    dic, data = ng.varian.read_lowmem(os.path.join(DATA_DIR, "agilent_3d"))
    assert data.shape == (128, 88, 1250)
    assert np.abs(data[0, 1, 2].real - 7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 11, 18].real - -9.36) <= 0.01
    assert np.abs(data[10, 11, 18].imag - -7.75) <= 0.01
    lowmem_write_readback(dic, data)


@attr(speed='slow')
def test_4d():
    """ reading/writing of 4D Varian fid file """
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic, data = ng.varian.read_fid(os.path.join(DATA_DIR, "agilent_4d/fid"),
                                   shape=(8, 12, 16, 1400), torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert np.abs(data[0, 1, 2, 3].real - 395.11) <= 0.01
    assert np.abs(data[0, 1, 2, 3].imag - 52.72) <= 0.01
    assert np.abs(data[3, 10, 11, 18].real - 51.81) <= 0.01
    assert np.abs(data[3, 10, 11, 18].imag - 16.01) <= 0.01
    write_fid_readback(dic, data, (8, 12, 16, 1400), 'r')


@attr(speed='slow')
def test_4d_lowmem():
    """ low memory reading/writing of 4D Varian fid file """
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic, data = ng.varian.read_fid_lowmem(
        os.path.join(DATA_DIR, "agilent_4d/fid"),
        shape=(8, 12, 16, 1400), torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert np.abs(data[0, 1, 2, 3].real - 395.11) <= 0.01
    assert np.abs(data[0, 1, 2, 3].imag - 52.72) <= 0.01
    assert np.abs(data[3, 10, 11, 18].real - 51.81) <= 0.01
    assert np.abs(data[3, 10, 11, 18].imag - 16.01) <= 0.01
    lowmem_fid_write_readback(dic, data, (8, 12, 16, 1400), 'r')
