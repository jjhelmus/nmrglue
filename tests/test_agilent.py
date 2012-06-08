""" Tests for the fileio.varian submodule """

import tempfile
import os
import shutil

from numpy.testing import assert_array_equal
import nmrglue as ng

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

# test
def test_1d():
    """ reading/writing of 1D Varian file """    
    dic, data = ng.varian.read(DATA_DIR + "agilent_1d")
    assert data.shape == (1500, )
    assert round(data[0].real, 2) == 91899.24
    assert round(data[0].imag, 2) == -1964.70
    assert round(data[1].real, 2) == 168844.25
    assert round(data[1].imag, 2) == -49503.41
    write_readback(dic, data)

def test_2d():
    """ reading/writing of 2D Varian file """    
    dic, data = ng.varian.read(DATA_DIR + "agilent_2d")
    assert data.shape == (332, 1500) 
    assert round(data[0, 1].real, 2) == 360.07 
    assert round(data[0, 1].imag, 2) == -223.20
    assert round(data[10, 18].real, 2) == 17.93
    assert round(data[10, 18].imag, 2) == -67.20
    write_readback(dic, data)

def test_2d_lowmem():
    """ low memory reading/writing of 2D Varian file """    
    dic, data = ng.varian.read_lowmem(DATA_DIR + "agilent_2d")
    assert data.shape == (332, 1500) 
    assert round(data[0, 1].real, 2) == 360.07 
    assert round(data[0, 1].imag, 2) == -223.20
    assert round(data[10, 18].real, 2) == 17.93
    assert round(data[10, 18].imag, 2) == -67.20
    lowmem_write_readback(dic, data)

def test_2d_tppi():
    """ reading/writing of 2D Varian file with TPPI encoding """    
    dic, data = ng.varian.read(DATA_DIR + "agilent_2d_tppi")
    assert data.shape == (600, 1400) 
    assert round(data[0, 1].real, 2) == -4589.29
    assert round(data[0, 1].imag, 2) == -1691.82
    assert round(data[10, 18].real, 2) == -166.62
    assert round(data[10, 18].imag, 2) == -594.73
    write_readback(dic, data)

def test_2d_tppi_lowmem():
    """ low memory reading/writing of 2D Varian file with TPPI encoding """    
    dic, data = ng.varian.read_lowmem(DATA_DIR + "agilent_2d_tppi")
    assert data.shape == (600, 1400) 
    assert round(data[0, 1].real, 2) == -4589.29
    assert round(data[0, 1].imag, 2) == -1691.82
    assert round(data[10, 18].real, 2) == -166.62
    assert round(data[10, 18].imag, 2) == -594.73
    lowmem_write_readback(dic, data)

def test_3d():
    """ reading/writing of 3D Varian file """    
    dic, data = ng.varian.read(DATA_DIR + "agilent_3d")
    assert data.shape == (128, 88, 1250)
    assert round(data[0, 1, 2].real, 2) == 7.98
    assert round(data[0, 1, 2].imag, 2) == 33.82
    assert round(data[10, 11, 18].real, 2) == -9.36
    assert round(data[10, 11, 18].imag, 2) == -7.75
    write_readback(dic, data)

def test_3d_lowmem():
    """ low memory reading/writing of 3D Varian file """    
    dic, data = ng.varian.read_lowmem(DATA_DIR + "agilent_3d") 
    assert data.shape == (128, 88, 1250)
    assert round(data[0, 1, 2].real, 2) == 7.98
    assert round(data[0, 1, 2].imag, 2) == 33.82
    assert round(data[10, 11, 18].real, 2) == -9.36
    assert round(data[10, 11, 18].imag, 2) == -7.75
    lowmem_write_readback(dic, data)

def test_4d():
    """ reading/writing of 4D Varian fid file """    
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic, data = ng.varian.read_fid(DATA_DIR + "agilent_4d/fid",
                shape=(8, 12, 16, 1400), torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert round(data[0, 1, 2, 3].real, 2) == 395.11 
    assert round(data[0, 1, 2, 3].imag, 2) == 52.72
    assert round(data[3, 10, 11, 18].real, 2) == 51.81
    assert round(data[3, 10, 11, 18].imag, 2) == 16.01
    write_fid_readback(dic, data, (8, 12, 16, 1400),'r')

def test_4d_lowmem():
    """ low memory reading/writing of 4D Varian fid file """    
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic, data = ng.varian.read_fid_lowmem(DATA_DIR + "agilent_4d/fid",
                shape=(8, 12, 16, 1400),torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert round(data[0, 1, 2, 3].real, 2) == 395.11 
    assert round(data[0, 1, 2, 3].imag, 2) == 52.72
    assert round(data[3, 10, 11, 18].real, 2) == 51.81
    assert round(data[3, 10, 11, 18].imag, 2) == 16.01
    lowmem_fid_write_readback(dic, data, (8, 12, 16, 1400),'r')
