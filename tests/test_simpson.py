
import nmrglue.fileio.simpson as simpson
import numpy as np
from numpy.testing import assert_allclose, assert_raises
import os.path

from setup import DATA_DIR

DD_1D = os.path.join(DATA_DIR, 'simpson_1d')
DD_2D = os.path.join(DATA_DIR, 'simpson_2d')

def test_1d_time():
    """ reading 1D time domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_1D, '1d_text.fid'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_1D, '1d_bin.fid'))
    xreim_units, xreim_data = simpson.read(os.path.join(DD_1D, '1d_ftext.fid'))
    rd, rawbin_data = simpson.read(os.path.join(DD_1D, '1d_rawbin.fid'), spe=False, ndim=1)
  
    # check data in text file
    assert text_data.shape == (4096, )
    assert text_data.dtype == 'complex64'
    assert round(text_data[0].real, 2) == 2.0
    assert round(text_data[0].imag, 2) == 0.0
    assert round(text_data[1].real, 2) == 1.78
    assert round(text_data[1].imag, 2) == -0.01

    # data in all files should be close
    assert np.allclose(rawbin_data, text_data)
    assert np.allclose(rawbin_data, bin_data)
    assert np.allclose(rawbin_data, xreim_data)


def test_1d_freq():
    """ reading 1D freq domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_1D, '1d_text.spe'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_1D, '1d_bin.spe'))
    xreim_units, xreim_data = simpson.read(os.path.join(DD_1D, '1d_ftext.spe'))
    rd, rawbin_data = simpson.read(os.path.join(DD_1D, '1d_rawbin.spe'), spe=True, ndim=1)
    
    # check data in text file
    assert text_data.shape == (4096, )
    assert text_data.dtype == 'complex64'
    assert round(text_data[2048].real, 2) == 40.34
    assert round(text_data[2048].imag, 2) == -1.51
    assert round(text_data[2049].real, 2) == 39.58
    assert round(text_data[2049].imag, 2) == -3.97


    # data in all file should be close
    assert np.allclose(rawbin_data, text_data)
    assert np.allclose(rawbin_data, bin_data)
    assert np.allclose(rawbin_data, xreim_data)

def test_2d_time():
    """ reading 2D time domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_2D, '2d_text.fid'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_2D, '2d.fid'))
    xyreim_units, xyreim_data = simpson.read(os.path.join(DD_2D, '2d_ftext.fid'))
    rd, rawbin_data = simpson.read(os.path.join(DD_2D, '2d_raw.fid'), NP=128, NI=48, ndim=2, 
                                    spe=False)
   
    # check data in text file
    assert text_data.shape == (48, 128)
    assert text_data.dtype == 'complex64'
    assert round(text_data[0,0].real, 2) == 1.00
    assert round(text_data[0,0].imag, 2) == 0.03
    assert round(text_data[0,1].real, 2) == 0.75
    assert round(text_data[0,1].imag, 2) == 0.59
    assert round(text_data[1,0].real, 2) == 0.89
    assert round(text_data[1,0].imag, 2) == 0.03

    # data in all files should be close
    assert np.allclose(rawbin_data, text_data)
    assert np.allclose(rawbin_data, bin_data)
    assert np.allclose(rawbin_data, xyreim_data)

def test_2d_freq():
    """ reading 2D freq domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_2D, '2d_text.spe'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_2D, '2d.spe'))
    xyreim_units, xyreim_data = simpson.read(os.path.join(DD_2D, '2d_ftext.spe'))
    rd, rawbin_data = simpson.read(os.path.join(DD_2D, '2d_raw.spe'), ndim=2, NP=256, 
                                    NI=512, spe=True)
    
    # check data in text file
    assert text_data.shape == (512, 256)
    assert text_data.dtype == 'complex64'
    assert round(text_data[4, 150].real, 2) == 0.29
    assert round(text_data[4, 150].imag, 2) == 0.34
    assert round(text_data[4, 151].real, 2) == 0.13
    assert round(text_data[4, 151].imag, 2) == 0.16
    assert round(text_data[5, 150].real, 2) == 0.41
    assert round(text_data[5, 150].imag, 2) == 0.14
 
    # data in text, bin and xyreim files should all be close
    assert np.allclose(text_data, bin_data)
    assert np.allclose(text_data, xyreim_data)
    
    # rawbin should be close except for first point along each vector
    assert np.allclose(rawbin_data[:, 1:], text_data[:, 1:])

def test_exceptions_read():
    """ raising exceptions due to missing read parameters """
    
    # missing spe parameter
    assert_raises(ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'))

    # missing ndim parameter
    assert_raises(ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'), spe=False)

    # missing NP/NI parameter
    assert_raises(ValueError, simpson.read, os.path.join(DD_2D, '2d_raw.fid'), spe=False,
            ndim=2)

    # bad ftype
    assert_raises(ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'), ftype='a')
