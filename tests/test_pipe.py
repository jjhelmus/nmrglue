""" Tests for the fileio.pipe submodule """

from __future__ import print_function
from nmrglue.fileio.pipe import fdata2dic

import tempfile
import os
import glob
import io

import numpy as np
from numpy.testing import assert_array_equal
import nmrglue as ng

from setup import DATA_DIR

# subroutines

def write_readback(dic, data):
    """ Write out a NMRPipe file and read back in. """
    # write out and read back
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, dic, data)
    rdic, rdata = ng.pipe.read(tf)
    os.remove(tf)
    assert_array_equal(data, rdata)
    assert dic == rdic


def write_readback_3d(dic, data):
    """ Write out a 3D NMRPipe file and read back in. """
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, dic, data)
    rdic, rdata = ng.pipe.read(tf)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)
    assert_array_equal(data, rdata)
    assert dic == rdic


def lowmem_write_readback_3d(dic, data):
    """ Lowmemory write out and readback of 3D NMRPipe file. """
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, dic, data)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    # check first values are the same
    s = tuple([0] * data.ndim)
    assert data[s] == rdata[s]
    assert dic == rdic
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)


def lowmem_write_readback_4d(dic, data):
    """ Lowmemory write out and readback of 4D NMRPipe file. """
    # write out and read back
    tf = tempfile.mktemp(dir=".") + "%02d%03d"
    ng.pipe.write_lowmem(tf, dic, data)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    # check first values are the same
    s = tuple([0] * data.ndim)
    assert data[s] == rdata[s]
    assert dic == rdic
    for f in glob.glob(tf[:-8] + "*"):
        os.remove(f)


def lowmem_write_readback(dic, data):
    """ Lowmemory write out and readback of NMRPipe file. """
    # lowmemory write out and read back
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, dic, data)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    # check first values are the same
    s = tuple([0] * data.ndim)
    assert data[s] == rdata[s]
    assert dic == rdic
    os.remove(tf)


def check_ppm_limits(dic, data, dim, limits):
    """ Check PPM Limits """
    uc0 = ng.pipe.make_uc(dic, data, dim=dim)
    climits = [round(i, 2) for i in uc0.ppm_limits()]
    print(limits)
    print(climits)
    assert limits == climits

def read_with_bytes_or_buffer(filename):
    """ Check reading pipe files from filename, io.BytesIO or bytes buffer """
    dic, data = ng.pipe.read(filename)
    # read bytes
    with open(filename, "rb") as binary_stream:
        data_bytes = binary_stream.read()
        bdic, bdata = ng.pipe.read(data_bytes)
    # read from io.BytesIO
    with open(filename, "rb") as binary_stream:
        bdic2, bdata2 = ng.pipe.read(binary_stream)
    assert dic == bdic == bdic2
    assert_array_equal(data, bdata)
    assert_array_equal(data, bdata2)

# tests
def test_get_fdata_bytes():
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid") 
    with open(data_path, "rb") as binary_stream:
        bytes_stream = binary_stream.read()
        bdata = ng.fileio.pipe.get_fdata(bytes_stream)
    data = ng.fileio.pipe.get_fdata(data_path)
    assert_array_equal(data, bdata)


def test_get_data_bytes():
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid") 
    with open(data_path, "rb") as binary_stream:
        bytes_stream = binary_stream.read()
        bdata = ng.fileio.pipe.get_data(bytes_stream)
    data = ng.fileio.pipe.get_data(data_path)
    assert_array_equal(data, bdata)


def test_fdata2dic_bytes():
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid") 
    with open(data_path, "rb") as binary_stream:
        bytes_stream = binary_stream.read()
        bfdata = ng.fileio.pipe.get_fdata(bytes_stream)
    fdata = ng.fileio.pipe.get_fdata(data_path)
    dic = fdata2dic(fdata)
    bdic = fdata2dic(bfdata)
    assert dic == bdic


def test_fdata_data_bytes():
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid") 
    with open(data_path, "rb") as binary_stream:
        bytes_stream = binary_stream.read()
        bdic_data, bdata = ng.fileio.pipe.get_fdata_data(bytes_stream)
    dic_data, data = ng.fileio.pipe.get_fdata_data(data_path)
    assert_array_equal(dic_data, bdic_data)
    assert_array_equal(data, bdata)


def test_read_bytes_1d_time():
    """ reading NMRPipe data from io.BytesIO binary stream or bytes buffer"""
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid")
    read_with_bytes_or_buffer(data_path)


def test_read_bytes_1d_freq():
    """ reading NMRPipe data from io.BytesIO binary stream or bytes buffer"""
    data_path = os.path.join(DATA_DIR, "nmrpipe_1d", "test.ft")
    read_with_bytes_or_buffer(data_path)


def test_read_bytes_2d_time():
    """ reading NMRPipe data from io.BytesIO binary stream """
    data_path = os.path.join(DATA_DIR, "nmrpipe_2d", "test.fid")
    read_with_bytes_or_buffer(data_path)


def test_read_bytes_2d_freq():
    """ reading NMRPipe data from io.BytesIO binary stream """
    data_path = os.path.join(DATA_DIR, "nmrpipe_2d", "test.ft2")
    read_with_bytes_or_buffer(data_path)


def test_read_bytes_3d_time():
    """ reading NMRPipe data from io.BytesIO binary stream """
    data_path = os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.fid")
    read_with_bytes_or_buffer(data_path)


def test_read_bytes_3d_freq():
    """ reading NMRPipe data from io.BytesIO binary stream """
    data_path = os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.ft3")
    read_with_bytes_or_buffer(data_path)


def test_1d_time():
    """reading/writing of 1D NMRPipe time domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_1d", "test.fid"))
    assert data.shape == (1500, )
    assert data.dtype == 'complex64'
    assert np.abs(data[0].real - 91899.24) <= 0.01
    assert np.abs(data[0].imag - 1964.70) <= 0.01
    assert np.abs(data[1].real - 168844.25) <= 0.01
    assert np.abs(data[1].imag - 49503.41) <= 0.01
    write_readback(dic, data)


def test_1d_freq():
    """reading/writing of 1D NMRPipe freq domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_1d", "test.ft"))
    assert data.shape == (4096, )
    assert data.dtype == 'float32'
    assert np.abs(data[0] - -63790.) <= 1.
    assert np.abs(data[1] - -63159.9) <= 0.1
    assert np.abs(data[100] - -29308.) <= 1.
    write_readback(dic, data)
    check_ppm_limits(dic, data, 0, [297.92, -99.82])


def test_1d_cut():
    """reading/writing of 1D NMRPipe with EXTracted region"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_1d", "test_cut.ft"))
    assert data.shape == (2766, )
    assert data.dtype == 'float32'
    assert np.abs(data[0] - -12123.7) <= 0.1
    assert np.abs(data[1] - -8979) <= 1.
    assert np.abs(data[100] - -7625.3) <= 0.1
    write_readback(dic, data)
    check_ppm_limits(dic, data, 0, [278.59, 10.03])


def test_2d_time():
    """reading/writing of 2D NMRPipe time domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test.fid"))
    assert data.shape == (332, 1500)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - 223.20) <= 0.01
    assert np.abs(data[10, 22].real - -26.76) <= 0.01
    assert np.abs(data[10, 22].imag - 42.67) <= 0.01
    write_readback(dic, data)


def test_2d_freq():
    """reading/writing of 2D NMRPipe freq domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test.ft2"))
    assert data.shape == (2048, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1] - 1601.83) <= 0.01
    assert np.abs(data[10, 22] - 3079.44) <= 0.01
    write_readback(dic, data)
    check_ppm_limits(dic, data, 0, [174.84, 65.21])
    check_ppm_limits(dic, data, 1, [253.90, -143.80])


def test_2d_time_tran():
    """reading/writing of TransPosed 2D NMRPipe time domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test_tp.fid"))
    assert data.shape == (3000, 166)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1].real - 99.46) <= 0.01
    assert np.abs(data[0, 1].imag - -80.63) <= 0.01
    assert np.abs(data[10, 22].real - 34.66) <= 0.01
    assert np.abs(data[10, 22].imag - 35.00) <= 0.01
    write_readback(dic, data)


def test_2d_freq_tran():
    """reading/writing of TransPosed 2D NMRPipe freq domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test_tp.ft2"))
    assert data.shape == (4096, 2048)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1] - -1525.10) <= 0.01
    assert np.abs(data[10, 22] - 1731.94) <= 0.01
    write_readback(dic, data)
    check_ppm_limits(dic, data, 0, [253.90, -143.80])
    check_ppm_limits(dic, data, 1, [174.84, 65.21])


def test_2d_time_lowmem():
    """lowmemory reading/writing of 2D NMRPipe time domain file"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test.fid"))
    assert data.shape == (332, 1500)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1].real - 360.07) <= 0.01
    assert np.abs(data[0, 1].imag - 223.20) <= 0.01
    assert np.abs(data[10, 22].real - -26.76) <= 0.01
    assert np.abs(data[10, 22].imag - 42.67) <= 0.01
    lowmem_write_readback(dic, data)


def test_2d_freq_lowmem():
    """lowmemory reading/writing of 2D NMRPipe freq domain file"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_2d", "test.ft2"))
    assert data.shape == (2048, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1] - 1601.83) <= 0.01
    assert np.abs(data[10, 22] - 3079.44) <= 0.01
    lowmem_write_readback(dic, data)
    check_ppm_limits(dic, data, 0, [174.84, 65.21])
    check_ppm_limits(dic, data, 1, [253.90, -143.80])


def test_3d_time():
    """reading/writing of 3D NMRPipe time domain data"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "data", "test%03d.fid"))
    sdic, sdata = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "data", "test001.fid"))
    assert data.shape == (128, 88, 1250)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2].real - -7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 22, 5].real - 15.71) <= 0.01
    assert np.abs(data[10, 22, 5].imag - 15.1) <= 0.01

    # and the first slice
    assert sdata.shape == (88, 1250)
    assert sdata.dtype == 'complex64'
    assert np.abs(sdata[1, 2].real - -7.98) <= 0.01
    assert np.abs(sdata[1, 2].imag - 33.82) <= 0.01
    assert np.abs(sdata[22, 5].real - 22.65) <= 0.01
    assert np.abs(sdata[22, 5].imag - 13.65) <= 0.01

    # slice/data matching
    assert_array_equal(data[0], sdata)

    write_readback_3d(dic, data)


def test_3d_freq():
    """reading/writing of 3D NMRPipe freq domain data"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))
    sdic, sdata = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test001.ft3"))
    assert data.shape == (128, 128, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2] - 25980.13) <= 0.01
    assert np.abs(data[10, 22, 5] - 1561.09) <= 0.01
    check_ppm_limits(dic, data, 0, [78.10, 34.24])
    check_ppm_limits(dic, data, 1, [147.42, 93.01])
    check_ppm_limits(dic, data, 2, [254.92, -142.83])

    # and the first slice
    assert sdata.shape == (128, 4096)
    assert sdata.dtype == 'float32'
    assert np.abs(sdata[1, 2] - 25980) <= 1.
    assert np.abs(sdata[22, 5] - -8336) <= 1.
    check_ppm_limits(sdic, sdata, 0, [147.42, 93.01])
    check_ppm_limits(sdic, sdata, 1, [254.92, -142.83])

    # slice/data matching
    assert_array_equal(data[0], sdata)

    write_readback_3d(dic, data)


def test_3d_time_lowmem():
    """lowmemory reading/writing of 3D NMRPipe time domain data"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "data", "test%03d.fid"))
    assert data.shape == (128, 88, 1250)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2].real - -7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 22, 5].real - 15.71) <= 0.01
    assert np.abs(data[10, 22, 5].imag - 15.1) <= 0.01
    lowmem_write_readback_3d(dic, data)


def test_3d_freq_lowmem():
    """lowmemory reading/writing of 3D NMRPipe freq domain data"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))
    assert data.shape == (128, 128, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2] - 25980.13) <= 0.01
    assert np.abs(data[10, 22, 5] - 1561.09) <= 0.01
    check_ppm_limits(dic, data, 0, [78.10, 34.24])
    check_ppm_limits(dic, data, 1, [147.42, 93.01])
    check_ppm_limits(dic, data, 2, [254.92, -142.83])
    lowmem_write_readback_3d(dic, data)


def test_3d_stream_time():
    """reading/writing of 3D NMRPipe data stream time domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.fid"))
    assert data.shape == (128, 88, 1250)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2].real - -7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 22, 5].real - 15.71) <= 0.01
    assert np.abs(data[10, 22, 5].imag - 15.1) <= 0.01
    write_readback(dic, data)


def test_3d_stream_freq():
    """reading/writing of 3D NMRPipe data stream freq domain file"""
    dic, data = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.ft3"))
    assert data.shape == (128, 128, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2] - 25980.13) <= 0.01
    assert np.abs(data[10, 22, 5] - 1561.09) <= 0.01
    check_ppm_limits(dic, data, 0, [78.10, 34.24])
    check_ppm_limits(dic, data, 1, [147.42, 93.01])
    check_ppm_limits(dic, data, 2, [254.92, -142.83])
    write_readback(dic, data)


def test_3d_stream_time_lowmem():
    """lowmemory reading/writing of 3D NMRPipe data stream time domain file"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.fid"))
    assert data.shape == (128, 88, 1250)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2].real - -7.98) <= 0.01
    assert np.abs(data[0, 1, 2].imag - 33.82) <= 0.01
    assert np.abs(data[10, 22, 5].real - 15.71) <= 0.01
    assert np.abs(data[10, 22, 5].imag - 15.1) <= 0.01
    lowmem_write_readback(dic, data)


def test_3d_stream_freq_lowmem():
    """lowmemory reading/writing of 3D NMRPipe data stream freq domain file"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "full3D.ft3"))
    assert data.shape == (128, 128, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2] - 25980.13) <= 0.01
    assert np.abs(data[10, 22, 5] - 1561.09) <= 0.01
    check_ppm_limits(dic, data, 0, [78.10, 34.24])
    check_ppm_limits(dic, data, 1, [147.42, 93.01])
    check_ppm_limits(dic, data, 2, [254.92, -142.83])
    lowmem_write_readback(dic, data)


def test_3d_slicing():
    """lowmemory 3D slicing"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))
    fdic, fdata = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))

    # null indexing
    assert_array_equal(fdata[1:0, 0, 0], data[1:0, 0, 0])
    assert_array_equal(fdata[50:0, 20:19, 12:12], data[50:0, 20:19, 12:12])
    assert_array_equal(fdata[64:64], data[64:64])

    # single element
    assert_array_equal(fdata[0, 0, 0], data[0, 0, 0])

    # incomplete indexing
    assert_array_equal(fdata[5], data[5])
    assert_array_equal(fdata[6, 7], data[6, 7])

    # ellipsis
    assert_array_equal(fdata[..., 7], data[..., 7])
    assert_array_equal(fdata[..., 2, 100], data[..., 2, 100])
    assert_array_equal(fdata[..., 2:10, :], data[..., 2:10, :])

    # simple slicing
    assert_array_equal(fdata[0, 1:4, 0], data[0, 1:4, 0])
    assert_array_equal(fdata[0:10, 0:12, 0:165], data[0:10, 0:12, 0:165])
    assert_array_equal(fdata[5:12, 5:8, 100:200], data[5:12, 5:8, 100:200])

    # missing endpoint
    assert_array_equal(fdata[50:, 0, :10], data[50:, 0, :10])

    # overindexed
    assert_array_equal(fdata[0:90, 50:100, 10:20], data[0:90, 50:100, 10:20])

    # extended slicing
    assert_array_equal(fdata[2:8:2, 0:20:7, :50:5], data[2:8:2, 0:20:7, :50:5])

    # negative extended slicing
    assert_array_equal(fdata[0, 0, 100:50:-1], data[0, 0, 100:50:-1])
    assert_array_equal(fdata[20:5:-3, 90:10:-7, 0], data[20:5:-3, 90:10:-7, 0])

    # overindexing with negative extended slicing
    assert_array_equal(fdata[78:1:-1, 0, 0], data[78:1:-1, 0, 0])
    assert_array_equal(fdata[0, 90:0:-1, 0], data[0, 90:0:-1, 0])

    # negative extended slicing with implicit ends
    assert_array_equal(fdata[0::-1, 0, 0], data[0::-1, 0, 0])
    assert_array_equal(fdata[0, :10:-2, 0], data[0, :10:-2, 0])
    assert_array_equal(fdata[:10:-3, 0, 0], data[:10:-3, 0, 0])
    assert_array_equal(fdata[0, 10::-1, 0], data[0, 10::-1, 0])
    assert_array_equal(fdata[::-1, ::-10, ::-7], data[::-1, ::-10, ::-7])

    # negative indexing
    assert_array_equal(fdata[-1, -2, -3], data[-1, -2, -3])
    assert_array_equal(fdata[-5:, -2:, -3], data[-5:, -2:, -3])
    assert_array_equal(fdata[-5:-1, -5:-2, -3], data[-5:-1, -5:-2, -3])
    assert_array_equal(fdata[-900:-1, -70:, 0], data[-900:-1, -70:, 0])
    assert_array_equal(fdata[-1:-100:-3, -1::-1, 0],
                       data[-1:-100:-3, -1::-1, 0])


def test_3d_tranpose():
    """lowmemory 3D axis transposing and swaps"""
    dic, data = ng.pipe.read_lowmem(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))
    fdic, fdata = ng.pipe.read(
        os.path.join(DATA_DIR, "nmrpipe_3d", "ft", "test%03d.ft3"))

    assert_array_equal(fdata.transpose()[0, 1, 2], data.transpose()[0, 1, 2])
    assert_array_equal(fdata.transpose((2, 0, 1))[0, 1, 2],
                       data.transpose((2, 0, 1))[0, 1, 2])
    assert_array_equal(fdata.swapaxes(0, 1)[0, 1, 2],
                       data.swapaxes(0, 1)[0, 1, 2])
    assert_array_equal(fdata.swapaxes(2, 0)[0, 1, 2],
                       data.swapaxes(2, 0)[0, 1, 2])


def test_4d_single_index_time():
    """reading/writing of 4D single-index NMRPipe time domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "time_1index",
                         "test%03d.fid")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "time_1index",
                         "test018.fid")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 12, 16, 1400)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2, 3].real - -395.11) <= 0.01
    assert np.abs(data[0, 1, 2, 3].imag - 52.72) <= 0.01
    assert np.abs(data[5, 9, 11, 987].real - -35.09) <= 0.01
    assert np.abs(data[5, 9, 11, 987].imag - 33.07) <= 0.01

    # check the slice
    assert sdata.shape == (16, 1400)
    assert sdata.dtype == 'complex64'
    assert np.abs(sdata[1, 2].real - 75.93) <= 0.01
    assert np.abs(sdata[1, 2].imag - 5.55) <= 0.01
    assert np.abs(sdata[7, 800].real - -8.93) <= 0.01
    assert np.abs(sdata[7, 800].imag - -10.24) <= 0.01

    # slice/data matching
    assert_array_equal(data[1, 5], sdata)

    lowmem_write_readback_3d(dic, data)


def test_4d_two_index_time():
    """reading/writing of 4D double-index NMRPipe time domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "time_2index",
                         "test%02d%03d.fid")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "time_2index",
                         "test02006.fid")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 12, 16, 1400)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2, 3].real - -395.11) <= 0.01
    assert np.abs(data[0, 1, 2, 3].imag - 52.72) <= 0.01
    assert np.abs(data[5, 9, 11, 987].real - -35.09) <= 0.01
    assert np.abs(data[5, 9, 11, 987].imag - 33.07) <= 0.01

    # check the slice
    assert sdata.shape == (16, 1400)
    assert sdata.dtype == 'complex64'
    assert np.abs(sdata[1, 2].real - 75.93) <= 0.01
    assert np.abs(sdata[1, 2].imag - 5.55) <= 0.01
    assert np.abs(sdata[7, 800].real - -8.93) <= 0.01
    assert np.abs(sdata[7, 800].imag - -10.24) <= 0.01

    # slice/data matching
    assert_array_equal(data[1, 5], sdata)

    lowmem_write_readback_4d(dic, data)


def test_4d_stream_time():
    """reading/writing of 4D data stream NMRPipe time domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "full4D.fid")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "time_2index",
                         "test02006.fid")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 12, 16, 1400)
    assert data.dtype == 'complex64'
    assert np.abs(data[0, 1, 2, 3].real - -395.11) <= 0.01
    assert np.abs(data[0, 1, 2, 3].imag - 52.72) <= 0.01
    assert np.abs(data[5, 9, 11, 987].real - -35.09) <= 0.01
    assert np.abs(data[5, 9, 11, 987].imag - 33.07) <= 0.01

    # check the slice
    assert sdata.shape == (16, 1400)
    assert sdata.dtype == 'complex64'
    assert np.abs(sdata[1, 2].real - 75.93) <= 0.01
    assert np.abs(sdata[1, 2].imag - 5.55) <= 0.01
    assert np.abs(sdata[7, 800].real - -8.93) <= 0.01
    assert np.abs(sdata[7, 800].imag - -10.24) <= 0.01

    # slice/data matching
    assert_array_equal(data[1, 5], sdata)

    lowmem_write_readback(dic, data)


def test_4d_single_index_freq():
    """reading/writing of 4D single-index NMRPipe freq domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "ft_1index", "test%03d.ft4")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "ft_1index", "test053.ft4")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 16, 16, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2, 3] - -2703.98) <= 0.01
    assert np.abs(data[5, 9, 11, 891] - 5212.07) <= 0.01
    check_ppm_limits(dic, data, 0, [321.03, -65.77])
    check_ppm_limits(dic, data, 1, [321.03, -93.40])
    check_ppm_limits(dic, data, 2, [232.62, -16.04])
    check_ppm_limits(dic, data, 3, [298.92, -98.83])

    # check the slice
    assert sdata.shape == (16, 4096)
    assert sdata.dtype == 'float32'
    assert np.abs(sdata[1, 2] - 602.70) <= 0.01
    assert np.abs(sdata[12, 900] - 2717.60) <= 0.01
    check_ppm_limits(sdic, sdata, 0, [232.62, -16.04])
    check_ppm_limits(sdic, sdata, 1, [298.92, -98.83])

    # slice/data matching
    assert_array_equal(data[3, 4], sdata)

    lowmem_write_readback_3d(dic, data)


def test_4d_two_index_freq():
    """reading/writing of 4D double-index NMRPipe freq domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "ft_2index",
                         "test%02d%03d.ft4")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "ft_2index", "test04005.ft4")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 16, 16, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2, 3] - -2703.98) <= 0.01
    assert np.abs(data[5, 9, 11, 891] - 5212.07) <= 0.01
    check_ppm_limits(dic, data, 0, [321.03, -65.77])
    check_ppm_limits(dic, data, 1, [321.03, -93.40])
    check_ppm_limits(dic, data, 2, [232.62, -16.04])
    check_ppm_limits(dic, data, 3, [298.92, -98.83])

    # check the slice
    assert sdata.shape == (16, 4096)
    assert sdata.dtype == 'float32'
    assert np.abs(sdata[1, 2] - 602.70) <= 0.01
    assert np.abs(sdata[12, 900] - 2717.60) <= 0.01
    check_ppm_limits(sdic, sdata, 0, [232.62, -16.04])
    check_ppm_limits(sdic, sdata, 1, [298.92, -98.83])

    # slice/data matching
    assert_array_equal(data[3, 4], sdata)

    lowmem_write_readback_4d(dic, data)


def test_4d_stream_index_freq():
    """reading/writing of 4D data stream NMRPipe freq domain file"""
    fmask = os.path.join(DATA_DIR, "nmrpipe_4d", "full4D.ft4")
    dic, data = ng.pipe.read_lowmem(fmask)

    fname = os.path.join(DATA_DIR, "nmrpipe_4d", "ft_2index", "test04005.ft4")
    sdic, sdata = ng.pipe.read(fname)

    assert data.shape == (8, 16, 16, 4096)
    assert data.dtype == 'float32'
    assert np.abs(data[0, 1, 2, 3] - -2703.98) <= 0.01
    assert np.abs(data[5, 9, 11, 891] - 5212.07) <= 0.01
    check_ppm_limits(dic, data, 0, [321.03, -65.77])
    check_ppm_limits(dic, data, 1, [321.03, -93.40])
    check_ppm_limits(dic, data, 2, [232.62, -16.04])
    check_ppm_limits(dic, data, 3, [298.92, -98.83])

    # check the slice
    assert sdata.shape == (16, 4096)
    assert sdata.dtype == 'float32'
    assert np.abs(sdata[1, 2] - 602.70) <= 0.01
    assert np.abs(sdata[12, 900] - 2717.60) <= 0.01
    check_ppm_limits(sdic, sdata, 0, [232.62, -16.04])
    check_ppm_limits(sdic, sdata, 1, [298.92, -98.83])

    # slice/data matching
    assert_array_equal(data[3, 4], sdata)

    lowmem_write_readback(dic, data)
