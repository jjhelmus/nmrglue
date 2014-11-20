""" Tests for the fileio.rnmrtk submodule """

from __future__ import print_function

import tempfile
import os

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng
from nose.plugins.attrib import attr

from setup import DATA_DIR


# subroutines
def write_readback(dic, data):
    """ write out and readback a RNMRTK file. """
    tf = tempfile.mktemp(suffix='.sec', dir='.')
    ng.rnmrtk.write(tf, dic, data)
    rdic, rdata = ng.rnmrtk.read(tf)
    os.remove(tf)
    os.remove(tf.replace('.sec', '.par'))
    assert_array_equal(data, rdata)
    assert dic == rdic


def lowmem_write_readback(dic, data):
    """ lowmemory write out and readback a RNMTRK file. """
    tf = tempfile.mktemp(suffix='.sec', dir='.')
    ng.rnmrtk.write_lowmem(tf, dic, data)
    rdic, rdata = ng.rnmrtk.read_lowmem(tf)
    # check value [0,1,...]
    s = tuple(range(data.ndim))
    assert data[s] == rdata[s]
    assert dic == rdic
    print(tf)
    os.remove(tf)
    os.remove(tf.replace('.sec', '.par'))


# tests


@attr(speed='fast')
def test_1d_time():
    """ reading/writing of 1D RNMRTK time domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, 'rnmrtk_1d', 'time_1d.sec'))
    assert data.shape == (1500, )
    assert np.abs(data[0].real - 91899.24) <= 0.01
    assert np.abs(data[0].imag - -1964.70) <= 0.01
    assert np.abs(data[1].real - 168844.25) <= 0.01
    assert np.abs(data[1].imag - -49503.41) <= 0.01
    assert dic['sw'][0] == 50000.0
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 99.0
    write_readback(dic, data)


@attr(speed='fast')
def test_1d_freq():
    """ reading/writing of 1D RNMRTK frequency domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, 'rnmrtk_1d', 'freq_1d.sec'))
    assert data.shape == (4096, )
    assert np.abs(data[0] - -1726.76) <= 0.01
    assert np.abs(data[1] - -1702.76) <= 0.01
    assert dic['sw'][0] == 50000.0
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 99.0
    write_readback(dic, data)


@attr(speed='fast')
def test_2d_time():
    """ reading/writing of 2D RNMRTK time domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, "rnmrtk_2d", "time_2d.sec"))
    assert data.shape == (332, 1500)
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - -223.20) <= 0.01
    assert np.abs(data[10, 18].real - 17.93) <= 0.01
    assert np.abs(data[10, 18].imag - -67.20) <= 0.01
    assert dic['sw'][1] == 50000.0
    assert dic['sf'][1] == 125.69
    assert dic['ppm'][1] == 55.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 50.65
    assert dic['ppm'][0] == 120.0
    write_readback(dic, data)


@attr(speed='fast')
def test_2d_freq():
    """ reading/writing of 2D RNMRTK frequency domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, "rnmrtk_2d", "freq_2d.sec"))
    assert data.shape == (2048, 4096)
    assert np.abs(data[0, 1] - -.19) <= 0.01
    assert np.abs(data[10, 18] - 0.88) <= 0.01
    assert dic['sw'][1] == 50000.0
    assert dic['sf'][1] == 125.69
    assert dic['ppm'][1] == 55.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 50.65
    assert dic['ppm'][0] == 120.0
    write_readback(dic, data)


@attr(speed='fast')
def test_2d_time_lowmem():
    """ low memory reading/writing of 2D RNMRTK time domain file """
    dic, data = ng.rnmrtk.read_lowmem(
        os.path.join(DATA_DIR, "rnmrtk_2d", "time_2d.sec"))
    assert data.shape == (332, 1500)
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - -223.20) <= 0.01
    assert np.abs(data[10, 18].real - 17.93) <= 0.01
    assert np.abs(data[10, 18].imag - -67.20) <= 0.01
    assert dic['sw'][1] == 50000.0
    assert dic['sf'][1] == 125.69
    assert dic['ppm'][1] == 55.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 50.65
    assert dic['ppm'][0] == 120.0
    lowmem_write_readback(dic, data)


@attr(speed='fast')
def test_2d_freq_lowmem():
    """ low memory reading/writing of 2D RNMRTK frequency domain file """
    dic, data = ng.rnmrtk.read_lowmem(
        os.path.join(DATA_DIR, "rnmrtk_2d", "freq_2d.sec"))
    assert data.shape == (2048, 4096)
    assert np.abs(data[0, 1] - -.19) <= 0.01
    assert np.abs(data[10, 18] - 0.88) <= 0.01
    assert dic['sw'][1] == 50000.0
    assert dic['sf'][1] == 125.69
    assert dic['ppm'][1] == 55.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 50.65
    assert dic['ppm'][0] == 120.0
    lowmem_write_readback(dic, data)


@attr(speed='slow')
def test_3d_time():
    """ reading/writing of 3D RNMRTK time domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, "rnmrtk_3d", "time_3d.sec"))
    assert data.shape == (128, 88, 1250)
    assert np.abs(data[0, 1, 2].real - 7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 11, 18].real - -9.36) <= 0.01
    assert np.abs(data[10, 11, 18].imag - -7.75) <= 0.01
    assert dic['sw'][2] == 50000.0
    assert dic['sf'][2] == 125.68
    assert dic['ppm'][2] == 56.0
    assert dic['sw'][1] == 2777.778
    assert dic['sf'][1] == 50.65
    assert dic['ppm'][1] == 120.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 56.0
    write_readback(dic, data)


@attr(speed='slow')
def test_3d_freq():
    """ reading/writing of 3D RNMRTK frequency domain file """
    dic, data = ng.rnmrtk.read(
        os.path.join(DATA_DIR, "rnmrtk_3d", "freq_3d.sec"))
    assert data.shape == (128, 128, 4096)
    assert np.abs(data[0, 1, 2] - 3.23) <= 0.01
    assert np.abs(data[10, 11, 18] - 1.16) <= 0.01
    assert dic['sw'][2] == 50000.0
    assert dic['sf'][2] == 125.68
    assert dic['ppm'][2] == 56.0
    assert dic['sw'][1] == 2777.778
    assert dic['sf'][1] == 50.65
    assert dic['ppm'][1] == 120.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 56.0
    write_readback(dic, data)


@attr(speed='slow')
def test_3d_time_lowmem():
    """ low memory reading/writing of 3D RNMRTK time domain file """
    dic, data = ng.rnmrtk.read_lowmem(
        os.path.join(DATA_DIR, "rnmrtk_3d", "time_3d.sec"))
    assert data.shape == (128, 88, 1250)
    assert np.abs(data[0, 1, 2].real - 7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 11, 18].real - -9.36) <= 0.01
    assert np.abs(data[10, 11, 18].imag - -7.75) <= 0.01
    assert dic['sw'][2] == 50000.0
    assert dic['sf'][2] == 125.68
    assert dic['ppm'][2] == 56.0
    assert dic['sw'][1] == 2777.778
    assert dic['sf'][1] == 50.65
    assert dic['ppm'][1] == 120.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 56.0
    lowmem_write_readback(dic, data)


@attr(speed='slow')
def test_3d_freq_lowmem():
    """ low memory reading/writing of 3D RNMRTK frequency domain file """
    dic, data = ng.rnmrtk.read_lowmem(
        os.path.join(DATA_DIR, "rnmrtk_3d", "freq_3d.sec"))
    assert data.shape == (128, 128, 4096)
    assert np.abs(data[0, 1, 2] - 3.23) <= 0.01
    assert np.abs(data[10, 11, 18] - 1.16) <= 0.01
    assert dic['sw'][2] == 50000.0
    assert dic['sf'][2] == 125.68
    assert dic['ppm'][2] == 56.0
    assert dic['sw'][1] == 2777.778
    assert dic['sf'][1] == 50.65
    assert dic['ppm'][1] == 120.0
    assert dic['sw'][0] == 5555.556
    assert dic['sf'][0] == 125.68
    assert dic['ppm'][0] == 56.0
    lowmem_write_readback(dic, data)


@attr(speed='fast')
def test_3d_transpose():
    """ reading/writing of transposed 3D RNMRTK time domain file """

    dir_3d = os.path.join(DATA_DIR, 'rnmrtk_3d')
    # T1 T2 T3 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t1_t2_t3.sec"))
    assert data.shape == (128, 88, 36)
    assert np.abs(data[2, 6, 4].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)

    # T1 T3 T2 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t1_t3_t2.sec"))
    assert data.shape == (128, 72, 44)
    assert np.abs(data[2, 8, 3].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)

    # T2 T1 T3 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t2_t1_t3.sec"))
    assert data.shape == (88, 128, 36)
    assert np.abs(data[6, 2, 4].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)

    # T2 T3 T1 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t2_t3_t1.sec"))
    assert data.shape == (88, 72, 64)
    assert np.abs(data[6, 8, 1].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)

    # T3 T1 T2 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t3_t1_t2.sec"))
    assert data.shape == (72, 128, 44)
    assert np.abs(data[8, 2, 3].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)

    # T3 T2 T1 ordering
    dic, data = ng.rnmrtk.read(os.path.join(dir_3d, "time_3d_t3_t2_t1.sec"))
    assert data.shape == (72, 88, 64)
    assert np.abs(data[8, 6, 1].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    write_readback(dic, data)


@attr(speed='slow')
def test_3d_transpose_lowmem():
    """ low mem. reading/writing of transposed 3D RNMRTK time domain file """

    dir_3d = os.path.join(DATA_DIR, 'rnmrtk_3d')
    # T1 T2 T3 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t1_t2_t3.sec"))
    assert data.shape == (128, 88, 36)
    assert np.abs(data[2, 6, 4].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)

    # T1 T3 T2 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t1_t3_t2.sec"))
    assert data.shape == (128, 72, 44)
    assert np.abs(data[2, 8, 3].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)

    # T2 T1 T3 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t2_t1_t3.sec"))
    assert data.shape == (88, 128, 36)
    assert np.abs(data[6, 2, 4].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)

    # T2 T3 T1 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t2_t3_t1.sec"))
    assert data.shape == (88, 72, 64)
    assert np.abs(data[6, 8, 1].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)

    # T3 T1 T2 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t3_t1_t2.sec"))
    assert data.shape == (72, 128, 44)
    assert np.abs(data[8, 2, 3].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)

    # T3 T2 T1 ordering
    dic, data = ng.rnmrtk.read_lowmem(os.path.join(
        dir_3d, "time_3d_t3_t2_t1.sec"))
    assert data.shape == (72, 88, 64)
    assert np.abs(data[8, 6, 1].real - -1.82) <= 0.01
    assert dic['npts'] == [64, 44, 36]
    lowmem_write_readback(dic, data)
