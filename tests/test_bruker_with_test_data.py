""" Tests for the fileio.bruker submodule """


import tempfile
import os
import shutil
from pathlib import Path

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng

from setup import DATA_DIR


# subroutines
def dic_similar(dic1, dic2):
    """ Compared two Bruker parameter dictionaries"""
    if dic1.keys() != dic2.keys():
        print("Not same keys!")
        print(dic1.keys())
        print(dic2.keys())
    assert dic1.keys() == dic2.keys()
    for key in dic1.keys():
        if dic1[key] != dic2[key]:
            print(key, dic1[key], dic2[key])
        assert dic1[key] == dic2[key]
    return True


def write_readback(dic, data):
    """ Write out and readback a Bruker directory """
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, dic, data, write_procs=True)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(data, rdata)
    assert dic_similar(dic, rdic)
    shutil.rmtree(td)


def write_readback_pdata(dic, data, pdata_folder=False):
    """ Write out and readback a Bruker processed dataset """
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_pdata(td, dic, data, write_procs=True,
                          pdata_folder=pdata_folder)
    if pdata_folder:
        rdic, rdata = ng.bruker.read_pdata(os.path.join(td, 'pdata',
                                           str(pdata_folder)), scale_data=False)
    else:
        rdic, rdata = ng.bruker.read_pdata(td, scale_data=False)
    assert_array_equal(data, rdata)
    assert dic_similar(dic, rdic)
    shutil.rmtree(td)


def lowmem_write_readback(dic, data):
    """ Lowmemory write out and readback of a Bruker directory """
    # write out and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td, dic, data)
    rdic, rdata = ng.bruker.read_lowmem(td, read_procs=False)
    tup = tuple(range(data.ndim))
    assert_array_equal(data[tup], rdata[tup])
    assert dic_similar(dic, rdic)
    shutil.rmtree(td)


# tests
def test_jcamp1():
    """ reading/writing of JCAMP file 1"""
    dic = ng.bruker.read_jcamp(os.path.join(DATA_DIR, "bruker_1d", "acqus"))
    assert dic['LFILTER'] == 200
    assert len(dic['PRECHAN']) == 16
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_jcamp(dic, tf)
    ndic = ng.bruker.read_jcamp(tf)
    assert dic_similar(dic, ndic)
    os.remove(tf)


def test_jcamp2():
    """ reading/writing of JCAMP file 2"""
    dic = ng.bruker.read_jcamp(os.path.join(DATA_DIR, "bruker_2d", "acqu2s"))
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_jcamp(dic, tf)
    ndic = ng.bruker.read_jcamp(tf)
    assert dic_similar(dic, ndic)
    os.remove(tf)


def test_pprog():
    """ reading/writing of pulse program"""
    dic = ng.bruker.read_pprog(os.path.join(DATA_DIR, "bruker_3d",
                                            "pulseprogram"))
    assert dic['var']['LOOPC'] == '58'
    assert len(dic['incr']) == 4
    tf = tempfile.mktemp(dir='.')
    ng.bruker.write_pprog(tf, dic)
    ndic = ng.bruker.read_pprog(tf)
    assert dic_similar(dic, ndic)
    os.remove(tf)


def test_1d():
    """ reading/writing of 1D bruker data"""
    dic, data = ng.bruker.read(os.path.join(DATA_DIR, "bruker_1d"))
    pdic, pdata = ng.bruker.read(Path(DATA_DIR) / "bruker_1d")
    assert dic['FILE_SIZE'] == pdic["FILE_SIZE"]== 16384
    assert data.shape == (2048, )
    assert_array_equal(data,pdata)
    assert np.abs(data[20].real - -282.0) <= 0.01
    assert np.abs(data[20].imag - 14.0) <= 0.01
    assert np.abs(data[91].real - -840842.0) <= 0.01
    assert np.abs(data[91].imag - -1116591.0) <= 0.01
    write_readback(dic, data)


def test_2d():
    """ reading/writing of 2D bruker data"""
    dic, data = ng.bruker.read(os.path.join(DATA_DIR, "bruker_2d"))
    assert dic['FILE_SIZE'] == 3686400
    assert data.shape == (600, 768)
    assert np.abs(data[0, 40].real - 28.0) <= 0.01
    assert np.abs(data[0, 40].imag - -286.0) <= 0.01
    assert np.abs(data[13, 91].real - -7279.0) <= 0.01
    assert np.abs(data[13, 91].imag - -17680.0) <= 0.01
    write_readback(dic, data)


def test_2d_lowmem():
    """ lowmemory reading/writing of 2D bruker data"""
    dic, data = ng.bruker.read_lowmem(os.path.join(DATA_DIR, "bruker_2d"),
                                      read_procs=False)
    assert dic['FILE_SIZE'] == 3686400
    assert data.shape == (600, 768)
    assert np.abs(data[0, 40].real - 28.0) <= 0.01
    assert np.abs(data[0, 40].imag - -286.0) <= 0.01
    assert np.abs(data[13, 91].real - -7279.0) <= 0.01
    assert np.abs(data[13, 91].imag - -17680.0) <= 0.01
    lowmem_write_readback(dic, data)


def test_3d():
    """ reading/writing of 3D bruker data"""
    dic, data = ng.bruker.read(os.path.join(DATA_DIR, "bruker_3d"))
    assert dic['FILE_SIZE'] == 91226112
    assert data.shape == (116, 128, 768)
    assert np.abs(data[0, 0, 40].real - 18.0) <= 0.01
    assert np.abs(data[0, 0, 40].imag - -66.0) <= 0.01
    assert np.abs(data[5, 13, 91].real - 1138.0) <= 0.01
    assert np.abs(data[5, 13, 91].imag - 3482.0) <= 0.01
    write_readback(dic, data)


def test_3d_lowmem():
    """ low memory reading/writing of 3D bruker data"""
    dic, data = ng.bruker.read_lowmem(os.path.join(DATA_DIR, "bruker_3d"),
                                      read_procs=False)
    assert dic['FILE_SIZE'] == 91226112
    assert data.shape == (116, 128, 768)
    assert np.abs(data[0, 0, 40].real - 18.0) <= 0.01
    assert np.abs(data[0, 0, 40].imag - -66.0) <= 0.01
    assert np.abs(data[5, 13, 91].real - 1138.0) <= 0.01
    assert np.abs(data[5, 13, 91].imag - 3482.0) <= 0.01
    lowmem_write_readback(dic, data)


def test_read_pdata_1d():
    """ read processed 1D data """
    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'bruker_1d',
                                                  'pdata', '1'))
    assert dic['procs']['OFFSET'] == 13.03153
    assert dic['procs']['SF'] == 600.13
    assert dic['procs']['FT_mod'] == 6
    assert data[9] - 189610.5 <= 0.001
    assert data[644] - 398782.375 <= 0.001
    assert data[1144] - 288069.375 <= 0.001
    assert data[1486] - 281011.875 <= 0.001
    assert data[1708] - 170066.875 <= 0.001


def test_read_pdata_2d():
    """ read processed 2d data """
    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'bruker_2d',
                                                  'pdata', '1'))
    assert dic['procs']['OFFSET'] == 11.60683
    assert dic['procs']['SF'] == 800.13
    assert dic['proc2s']['OFFSET'] == 143.1681
    assert dic['proc2s']['SF'] == 81.076469
    assert data[2, 217] - 291066.5 <= 0.001
    assert data[10, 271] - 140808.375 <= 0.001
    assert data[24, 219] - 197628.75 <= 0.001
    assert data[405, 189] - 134437.75 <= 0.001
    assert data[507, 258] - 221842.125 <= 0.001


def test_write_pdata_1d():
    """ writing of processed 1D bruker data """
    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'bruker_1d',
                                     'pdata', '1'), read_acqus=False)
    write_readback_pdata(dic=dic, data=data)
    write_readback_pdata(dic=dic, data=data, pdata_folder=90)


def test_write_pdata_2d():
    """ writing of processed 2D bruker data """
    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'bruker_2d',
                                     'pdata', '1'), read_acqus=False)
    write_readback_pdata(dic=dic, data=data)
    write_readback_pdata(dic=dic, data=data, pdata_folder=90)


def test_read_nuslist():
    """ reading nuslist """
    with open("tmp_nuslist", "w") as f:
        f.write("""0 0\n10 20\n50 21\n9 8\n7 8\n20 20""")

    nuslist = ng.bruker.read_nuslist(fname="tmp_nuslist")
    assert nuslist == [(0, 0), (10, 20), (50, 21), (9, 8), (7, 8), (20, 20)]

    os.remove("tmp_nuslist")


def test_read_vdlist():
    """ reading vdlist """
    with open("tmp_vdlist", "w") as f:
        f.write("""1n\n10n\n50u\n20u\n30m\n50m\n1\n2\n""")

    vdlist = ng.bruker.read_vdlist(".", fname="tmp_vdlist")
    true_vdlist = [1e-9, 10e-9, 50e-6, 20e-6, 30e-3, 50e-3, 1.0, 2.0] 
    for i, j in zip(vdlist, true_vdlist):
        assert i -j < 1e-10

    os.remove("tmp_vdlist")