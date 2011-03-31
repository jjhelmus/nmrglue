import nmrglue as ng
import numpy as np

import tempfile
import os
import shutil
import glob

from numpy.testing import assert_array_equal

# subroutines
def dic_similar(dic1,dic2):
    if dic1.keys() != dic2.keys():
        print "Not same keys!"
    assert dic1.keys() == dic2.keys()
    for key in dic1.keys():
        if dic1[key] != dic2[key]:
            print key,dic1[key],dic2[key]
        assert dic1[key] == dic2[key]
    return True

def write_readback(dic,data):
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td,dic,data)
    rdic,rdata = ng.bruker.read(td)
    assert_array_equal(data,rdata)
    assert dic_similar(dic,rdic)
    shutil.rmtree(td)

def lowmem_write_readback(dic,data):
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td,dic,data)
    rdic,rdata = ng.bruker.read_lowmem(td)
    tup = tuple(range(data.ndim))
    assert_array_equal(data[tup],rdata[tup])
    assert dic_similar(dic,rdic)
    shutil.rmtree(td)


# tests
def test_jcamp1():
    """ reading/writing of JCAMP file 1"""
    dic = ng.bruker.read_jcamp("common_data/1d_bruker/acqus")
    assert dic['LFILTER'] == 200
    assert len(dic['PRECHAN']) == 16
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_jcamp(dic,tf)
    ndic = ng.bruker.read_jcamp(tf)
    assert dic_similar(dic,ndic)
    os.remove(tf)
    
def test_jcamp2():
    """ reading/writing of JCAMP file 2"""
    dic = ng.bruker.read_jcamp("common_data/2d_bruker/acqu2s")
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_jcamp(dic,tf)
    ndic = ng.bruker.read_jcamp(tf)
    assert dic_similar(dic,ndic)
    os.remove(tf)

def test_pprog():
    """ reading/writing of pulse program"""
    dic = ng.bruker.read_pprog("common_data/3d_bruker/pulseprogram")
    assert dic['var']['LOOPC'] == '58'
    assert len(dic['incr']) == 4
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_pprog(tf,dic)
    ndic = ng.bruker.read_pprog(tf)
    assert dic_similar(dic,ndic)
    os.remove(tf)

def test_1d():
    """ reading/writing of 1D bruker data"""
    dic,data = ng.bruker.read("common_data/1d_bruker/")
    assert dic['FILE_SIZE'] == 16384
    assert data.shape == (2048,)
    assert round(data[20].real,2) == -282.0
    assert round(data[20].imag,2) == 14.0
    assert round(data[91].real,2) == -840842.0
    assert round(data[91].imag,2) == -1116591.0
    write_readback(dic,data)


def test_2d():
    """ reading/writing of 2D bruker data"""
    dic,data = ng.bruker.read("common_data/2d_bruker/")
    assert dic['FILE_SIZE'] == 3686400
    assert data.shape == (600, 768)
    assert round(data[0,40].real,2) == 28.0
    assert round(data[0,40].imag,2) == -286.0
    assert round(data[13,91].real,2) == -7279.0
    assert round(data[13,91].imag,2) == -17680.0
    write_readback(dic,data)

def test_2d_lowmem():
    """ lowmemory reading/writing of 2D bruker data"""
    dic,data = ng.bruker.read_lowmem("common_data/2d_bruker/")
    assert dic['FILE_SIZE'] == 3686400
    assert data.shape == (600, 768)
    assert round(data[0,40].real,2) == 28.0
    assert round(data[0,40].imag,2) == -286.0
    assert round(data[13,91].real,2) == -7279.0
    assert round(data[13,91].imag,2) == -17680.0
    lowmem_write_readback(dic,data)


def test_3d():
    """ reading/writing of 3D bruker data"""
    dic,data = ng.bruker.read("common_data/3d_bruker/")
    assert dic['FILE_SIZE'] == 91226112
    assert data.shape == (116,128,768)
    assert round(data[0,0,40].real,2) == 18.0
    assert round(data[0,0,40].imag,2) == -66.0
    assert round(data[5,13,91].real,2) == 1138.0 
    assert round(data[5,13,91].imag,2) == 3482.0
    write_readback(dic,data)


def test_3d_lowmem():
    """ low memory reading/writing of 3D bruker data"""
    dic,data = ng.bruker.read_lowmem("common_data/3d_bruker/")
    assert dic['FILE_SIZE'] == 91226112
    assert data.shape == (116,128,768)
    assert round(data[0,0,40].real,2) == 18.0
    assert round(data[0,0,40].imag,2) == -66.0
    assert round(data[5,13,91].real,2) == 1138.0 
    assert round(data[5,13,91].imag,2) == 3482.0
    lowmem_write_readback(dic,data)
