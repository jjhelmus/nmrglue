import nmrglue as ng
import numpy as np

import tempfile
import os
import shutil
import glob

from numpy.testing import assert_array_equal

# subroutines

def write_readback(dic,data):
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td,dic,data)
    rdic,rdata = ng.varian.read(td)
    shutil.rmtree(td)
    assert_array_equal(data,rdata)
    assert dic == rdic

def lowmem_write_readback(dic,data):
    # lowmemory write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td,dic,data)
    rdic,rdata = ng.varian.read_lowmem(td)
    # check value [0,1,...]
    s = tuple(range(data.ndim))
    assert data[s] == rdata[s]
    assert dic == rdic
    shutil.rmtree(td)

def write_fid_readback(dic,data,shape,torder):
    # write out and readback a fid file
    tf = tempfile.mktemp(dir=".")
    ng.varian.write_fid(tf,dic,data,torder=torder)
    rdic,rdata = ng.varian.read_fid(tf,shape=shape,torder=torder)
    os.remove(tf)
    assert_array_equal(data,rdata)
    assert dic == rdic

def lowmem_fid_write_readback(dic,data,shape,torder):
    # lowmemory write out and readback
    tf = tempfile.mktemp(dir=".")
    ng.varian.write_fid_lowmem(tf,dic,data,torder=torder)
    rdic,rdata = ng.varian.read_fid_lowmem(tf,shape=shape,torder=torder)
    # check value [0,1,...]
    s = tuple(range(data.ndim))
    assert data[s] == rdata[s]
    assert dic == rdic
    os.remove(tf)


# test
def test_1d():
    """ reading/writing of 1D Varian file """    
    dic,data = ng.varian.read("common_data/1d_varian/")
    assert data.shape == (1500,)
    assert round(data[0].real,2) == 91899.24
    assert round(data[0].imag,2) == -1964.70
    assert round(data[1].real,2) == 168844.25
    assert round(data[1].imag,2) == -49503.41
    write_readback(dic,data)


def test_2d():
    """ reading/writing of 2D Varian file """    
    dic,data = ng.varian.read("common_data/2d_varian/")
    assert data.shape == (332, 1500) 
    assert round(data[0,1].real,2) == 360.07 
    assert round(data[0,1].imag,2) == -223.20
    assert round(data[10,18].real,2) == 17.93
    assert round(data[10,18].imag,2) == -67.20
    write_readback(dic,data)


def test_2d_lowmem():
    """ low memory reading/writing of 2D Varian file """    
    dic,data = ng.varian.read_lowmem("common_data/2d_varian/")
    assert data.shape == (332, 1500) 
    assert round(data[0,1].real,2) == 360.07 
    assert round(data[0,1].imag,2) == -223.20
    assert round(data[10,18].real,2) == 17.93
    assert round(data[10,18].imag,2) == -67.20
    lowmem_write_readback(dic,data)


def test_2d_tppi():
    """ reading/writing of 2D Varian file with TPPI encoding """    
    dic,data = ng.varian.read("common_data/2d_varian_tppi/")
    assert data.shape == (600, 1400) 
    assert round(data[0,1].real,2) == -4589.29
    assert round(data[0,1].imag,2) == -1691.82
    assert round(data[10,18].real,2) == -166.62
    assert round(data[10,18].imag,2) == -594.73
    write_readback(dic,data)


def test_2d_tppi_lowmem():
    """ low memory reading/writing of 2D Varian file with TPPI encoding """    
    dic,data = ng.varian.read_lowmem("common_data/2d_varian_tppi/")
    assert data.shape == (600, 1400) 
    assert round(data[0,1].real,2) == -4589.29
    assert round(data[0,1].imag,2) == -1691.82
    assert round(data[10,18].real,2) == -166.62
    assert round(data[10,18].imag,2) == -594.73
    lowmem_write_readback(dic,data)


def test_3d():
    """ reading/writing of 3D Varian file """    
    dic,data = ng.varian.read("common_data/3d_varian")
    assert data.shape == (128, 88, 1250)
    assert round(data[0,1,2].real,2) == 7.98
    assert round(data[0,1,2].imag,2) == 33.82
    assert round(data[10,11,18].real,2) == -9.36
    assert round(data[10,11,18].imag,2) == -7.75
    write_readback(dic,data)


def test_3d_lowmem():
    """ low memory reading/writing of 3D Varian file """    
    dic,data = ng.varian.read_lowmem("common_data/3d_varian") 
    assert data.shape == (128, 88, 1250)
    assert round(data[0,1,2].real,2) == 7.98
    assert round(data[0,1,2].imag,2) == 33.82
    assert round(data[10,11,18].real,2) == -9.36
    assert round(data[10,11,18].imag,2) == -7.75
    lowmem_write_readback(dic,data)


def test_4d():
    """ reading/writing of 4D Varian fid file """    
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic,data = ng.varian.read_fid("common_data/4d_varian/fid",
                shape=(8, 12, 16, 1400),torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert round(data[0,1,2,3].real,2) == 395.11 
    assert round(data[0,1,2,3].imag,2) == 52.72
    assert round(data[3,10,11,18].real,2) == 51.81
    assert round(data[3,10,11,18].imag,2) == 16.01
    write_fid_readback(dic,data,(8, 12, 16, 1400),'r')


def test_4d_lowmem():
    """ low memory reading/writing of 4D Varian fid file """    
    # since this is a fake 4D with no procpar we need to explicitly
    # provide the shape and trace ordering parameters
    dic,data = ng.varian.read_fid_lowmem("common_data/4d_varian/fid",
                shape=(8, 12, 16, 1400),torder='r')
    assert data.shape == (8, 12, 16, 1400)
    assert round(data[0,1,2,3].real,2) == 395.11 
    assert round(data[0,1,2,3].imag,2) == 52.72
    assert round(data[3,10,11,18].real,2) == 51.81
    assert round(data[3,10,11,18].imag,2) == 16.01
    lowmem_fid_write_readback(dic,data,(8, 12, 16, 1400),'r')




