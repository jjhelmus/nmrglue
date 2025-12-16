import nmrglue.fileio.simpson as simpson
import numpy as np
from numpy.testing import assert_allclose, assert_raises
import os.path
import pytest

from setup import DATA_DIR

DD_1D = os.path.join(DATA_DIR, 'simpson_1d')
DD_2D = os.path.join(DATA_DIR, 'simpson_2d')


def test_1d_time():
    """ reading 1D time domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_1D, '1d_text.fid'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_1D, '1d_bin.fid'))
    xreim_units, xreim_data = simpson.read(os.path.join(DD_1D, '1d_ftext.fid'))
    rd, rawbin_data = simpson.read(
        os.path.join(DD_1D, '1d_rawbin.fid'), spe=False, ndim=1)

    # check data in text file
    assert text_data.shape == (1, 4096)
    assert text_data.dtype == 'complex64'
    assert np.abs(text_data[0][0].real - 2.0) <= 0.01
    assert np.abs(text_data[0][0].imag - 0.0) <= 0.01
    assert np.abs(text_data[0][1].real - 1.78) <= 0.01
    assert np.abs(text_data[0][1].imag - -0.01) <= 0.01

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
    rd, rawbin_data = simpson.read(
        os.path.join(DD_1D, '1d_rawbin.spe'), spe=True, ndim=1)

    # check data in text file
    assert text_data.shape == (1, 4096)
    assert text_data.dtype == 'complex64'
    assert np.abs(text_data[0][2048].real - 40.34) <= 0.01
    assert np.abs(text_data[0][2048].imag - -1.51) <= 0.01
    assert np.abs(text_data[0][2049].real - 39.58) <= 0.01
    assert np.abs(text_data[0][2049].imag - -3.97) <= 0.01

    # data in all file should be close
    assert np.allclose(rawbin_data, text_data)
    assert np.allclose(rawbin_data, bin_data)
    assert np.allclose(rawbin_data, xreim_data)


def test_2d_time():
    """ reading 2D time domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_2D, '2d_text.fid'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_2D, '2d.fid'))
    xyreim_units, xyreim_data = simpson.read(
        os.path.join(DD_2D, '2d_ftext.fid'))
    rd, rawbin_data = simpson.read(
        os.path.join(DD_2D, '2d_raw.fid'), NP=128, NI=48, ndim=2, spe=False)

    # check data in text file
    assert text_data.shape == (48, 128)
    assert text_data.dtype == 'complex64'
    assert np.abs(text_data[0, 0].real - 1.00) <= 0.01
    assert np.abs(text_data[0, 0].imag - 0.03) <= 0.01
    assert np.abs(text_data[0, 1].real - 0.75) <= 0.01
    assert np.abs(text_data[0, 1].imag - 0.59) <= 0.01
    assert np.abs(text_data[1, 0].real - 0.89) <= 0.01
    assert np.abs(text_data[1, 0].imag - 0.03) <= 0.01

    # data in all files should be close
    assert np.allclose(rawbin_data, text_data)
    assert np.allclose(rawbin_data, bin_data)
    assert np.allclose(rawbin_data, xyreim_data)

@pytest.mark.slow
def test_2d_freq():
    """ reading 2D freq domain files """
    # read the text, binary, xreim, and rawbin data
    text_dic, text_data = simpson.read(os.path.join(DD_2D, '2d_text.spe'))
    bin_dic, bin_data = simpson.read(os.path.join(DD_2D, '2d.spe'))
    xyreim_units, xyreim_data = simpson.read(
        os.path.join(DD_2D, '2d_ftext.spe'))
    rd, rawbin_data = simpson.read(
        os.path.join(DD_2D, '2d_raw.spe'), ndim=2, NP=256, NI=512, spe=True)

    # check data in text file
    assert text_data.shape == (512, 256)
    assert text_data.dtype == 'complex64'
    assert np.abs(text_data[4, 150].real - 0.29) <= 0.01
    assert np.abs(text_data[4, 150].imag - 0.34) <= 0.01
    assert np.abs(text_data[4, 151].real - 0.13) <= 0.01
    assert np.abs(text_data[4, 151].imag - 0.16) <= 0.01
    assert np.abs(text_data[5, 150].real - 0.41) <= 0.01
    assert np.abs(text_data[5, 150].imag - 0.14) <= 0.01

    # data in text, bin and xyreim files should all be close
    assert np.allclose(text_data, bin_data)
    assert np.allclose(text_data, xyreim_data)

    # rawbin should be close except for first point along each vector
    assert np.allclose(rawbin_data[:, 1:], text_data[:, 1:])


def test_exceptions_read():
    """ raising exceptions due to missing read parameters """

    # missing spe parameter
    assert_raises(
        ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'))

    # missing ndim parameter
    assert_raises(
        ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'),
        spe=False)

    # missing NP/NI parameter
    assert_raises(
        ValueError, simpson.read, os.path.join(DD_2D, '2d_raw.fid'),
        spe=False, ndim=2)

    # bad ftype
    assert_raises(
        ValueError, simpson.read, os.path.join(DD_1D, '1d_rawbin.fid'),
        ftype='a')
