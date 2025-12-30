""" Unit tests for nmrglue/fileio/pipe.py module """

import os
import tempfile
import glob

from numpy.testing import assert_array_equal
import nmrglue as ng
import pytest

# NMRPipe files being tested, these are created by the script:
# create_test_data_nmrpipe.sh
# they are checked into the repository, so NMRPipe is not required to
# run the tests
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

NMRPIPE_1D_TIME = os.path.join(DATA_DIR, 'nmrpipe_1d_time.fid')
NMRPIPE_1D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_1d_freq.fid')
NMRPIPE_1D_EXT = os.path.join(DATA_DIR, 'nmrpipe_1d_ext.fid')

NMRPIPE_2D_TIME = os.path.join(DATA_DIR, 'nmrpipe_2d_time.fid')
NMRPIPE_2D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_2d_freq.ft2')
NMRPIPE_2D_TIME_TP = os.path.join(DATA_DIR, 'nmrpipe_2d_time_tp.fid')
NMRPIPE_2D_FREQ_TP = os.path.join(DATA_DIR, 'nmrpipe_2d_freq_tp.ft2')

NMRPIPE_3D_TIME = os.path.join(DATA_DIR, 'nmrpipe_3d_time.dir',
                               'nmrpipe_3d_time_%03d.fid')
NMRPIPE_3D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_3d_freq.dir',
                               'nmrpipe_3d_freq_%03d.ft3')
NMRPIPE_3D_TIME_STREAM = os.path.join(DATA_DIR, 'nmrpipe_3d_time.fid')
NMRPIPE_3D_FREQ_STREAM = os.path.join(DATA_DIR, 'nmrpipe_3d_freq.ft3')

NMRPIPE_4D_TIME_2 = os.path.join(DATA_DIR, 'nmrpipe_4d_time_2.dir',
                                 'nmrpipe_4d_time_%03d_%03d.fid')
NMRPIPE_4D_TIME_1 = os.path.join(DATA_DIR, 'nmrpipe_4d_time_1.dir',
                                 'nmrpipe_4d_time_%03d.fid')
NMRPIPE_4D_TIME_STREAM = os.path.join(DATA_DIR, 'nmrpipe_4d_time.fid')

NMRPIPE_4D_FREQ_1 = os.path.join(DATA_DIR, 'nmrpipe_4d_freq_1.dir',
                                 'nmrpipe_4d_freq_%03d.ft4')
NMRPIPE_4D_FREQ_2 = os.path.join(DATA_DIR, 'nmrpipe_4d_freq_2.dir',
                                 'nmrpipe_4d_freq_%03d_%03d.ft4')
NMRPIPE_4D_FREQ_STREAM = os.path.join(DATA_DIR, 'nmrpipe_4d_freq.ft4')
NMRPIPE_TABLE = os.path.join(DATA_DIR, 'test.tab')


def check_simple_roundtrip(dic, data, n_percents=0, lowmem=False):
    """ Check write/read roundtrip for simple data. """
    try:
        # create one or more temporary file name(s)
        base_fname = tempfile.mktemp(prefix='nmrgluetest', dir='.')
        tmpfname = base_fname + '_%03d' * n_percents

        # write data to it, read back and check that same
        if lowmem:
            ng.pipe.write_lowmem(tmpfname, dic, data)
            rdic, rdata = ng.pipe.read_lowmem(tmpfname)
            s = tuple([0] * data.ndim)
            assert data[s] == rdata[s]
        else:
            ng.pipe.write(tmpfname, dic, data)
            rdic, rdata = ng.pipe.read(tmpfname)
            assert_array_equal(data, rdata)
        assert dic == rdic

    finally:    # remove the temporary files
        for f in glob.glob(base_fname + "*"):
            os.remove(f)


def check_ppm_limits(dic, data, dim, limits):
    """ Check PPM Limits """
    uc0 = ng.pipe.make_uc(dic, data, dim=dim)
    l0, l1 = uc0.ppm_limits()
    assert round(l0, 2) == limits[0]
    assert round(l1, 2) == limits[1]


def test_1d_time():
    """ read/write of 1D NMRPipe time domain file. """
    dic, data = ng.pipe.read(NMRPIPE_1D_TIME)
    assert data.shape == (16, )
    assert data.dtype == 'complex64'
    assert data[0].real == 1.
    assert data[0].imag == -1.
    assert data[1].real == 2.
    assert data[1].imag == -2.
    check_simple_roundtrip(dic, data)


def test_1d_freq():
    """ read/write of 1D NMRPipe frequency domain file. """
    dic, data = ng.pipe.read(NMRPIPE_1D_FREQ)
    assert data.shape == (16, )
    assert data.dtype == 'float32'
    assert data[0].real == 1.
    assert data[1].real == 2.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [149.00, 55.25])


def test_1d_ext():
    """ read/write of 1D NMRPipe file with an EXTacted region. """
    dic, data = ng.pipe.read(NMRPIPE_1D_EXT)
    assert data.shape == (8, )
    assert data.dtype == 'float32'
    assert data[0].real == 1.
    assert data[1].real == 2.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [130.25, 86.50])


def test_2d_time():
    """ read/write of 2D NMRPipe time domain file. """
    dic, data = ng.pipe.read(NMRPIPE_2D_TIME)
    assert data.shape == (4, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == -1.
    assert data[0, 1].real == 2.
    assert data[0, 1].imag == -2.
    assert data[1, 0].real == 1
    assert data[1, 0].imag == -1
    check_simple_roundtrip(dic, data)


def test_2d_freq():
    """ read/write of 2D NMRPipe frequency domain file. """
    dic, data = ng.pipe.read(NMRPIPE_2D_FREQ)
    assert data.shape == (2, 8)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 2.
    assert data[1, 0] == 1.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [179.00, 99.00])
    check_ppm_limits(dic, data, 1, [54.70, -32.80])


def test_2d_time_tp():
    """ read/write of 2D NMRPipe time domain file transposed. """
    dic, data = ng.pipe.read(NMRPIPE_2D_TIME_TP)
    assert data.shape == (16, 2)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == 1.
    assert data[0, 1].real == 1.
    assert data[0, 1].imag == 1.
    assert data[1, 0].real == -1.
    assert data[1, 0].imag == -1.
    assert data[2, 0].real == 2.
    assert data[2, 0].imag == 2.
    check_simple_roundtrip(dic, data)


def test_2d_freq_tp():
    """ read/write of 2D NMRPipe frequency domain file transposed. """
    dic, data = ng.pipe.read(NMRPIPE_2D_FREQ_TP)
    assert data.shape == (8, 2)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 1.
    assert data[1, 0] == 2.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [54.70, -32.80])
    check_ppm_limits(dic, data, 1, [179.00, 99.00])


def test_2d_time_lowmem():
    """ lowmem read/write of 2D NMRPipe time domain file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_2D_TIME)
    assert data.shape == (4, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == -1.
    assert data[0, 1].real == 2.
    assert data[0, 1].imag == -2.
    assert data[1, 0].real == 1
    assert data[1, 0].imag == -1
    check_simple_roundtrip(dic, data, lowmem=True)


def test_2d_freq_lowmem():
    """ lowmem read/write of 2D NMRPipe frequency domain file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_2D_FREQ)
    assert data.shape == (2, 8)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 2.
    assert data[1, 0] == 1.
    check_simple_roundtrip(dic, data, lowmem=True)
    check_ppm_limits(dic, data, 0, [179.00, 99.00])
    check_ppm_limits(dic, data, 1, [54.70, -32.80])


def test_3d_time():
    """ read/write of 3D NMRPipe time domain file. """
    dic, data = ng.pipe.read(NMRPIPE_3D_TIME)
    assert data.shape == (4, 6, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0].real == 1.
    assert data[0, 0, 0].imag == -1.
    assert data[0, 0, 1].real == 2.
    assert data[0, 0, 1].imag == -2.
    assert data[0, 1, 0].real == 1
    assert data[0, 1, 0].imag == -1
    assert data[1, 0, 0].real == 1
    assert data[1, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, 1)


def test_3d_time_slice():
    """ read/write of 3D NMRPipe time domain first slice. """
    dic, data = ng.pipe.read(NMRPIPE_3D_TIME % (1))
    assert data.shape == (6, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == -1.
    assert data[0, 1].real == 2.
    assert data[0, 1].imag == -2.
    assert data[1, 0].real == 1
    assert data[1, 0].imag == -1
    check_simple_roundtrip(dic, data)


def test_3d_freq():
    """ read/write of 3D NMRPipe frequency domain file. """
    dic, data = ng.pipe.read(NMRPIPE_3D_FREQ)
    assert data.shape == (2, 3, 8)
    assert data.dtype == 'float32'
    assert data[0, 0, 0] == 1.
    assert data[0, 0, 1] == 1.
    assert data[0, 1, 0] == 1.
    assert data[1, 0, 0] == 2.
    check_simple_roundtrip(dic, data, 1)
    check_ppm_limits(dic, data, 0, [220.00, 120.00])
    check_ppm_limits(dic, data, 1, [152.33, 45.67])
    check_ppm_limits(dic, data, 2, [54.70, -32.80])


def test_3d_freq_slice():
    """ read/write of 3D NMRPipe frequency domain first slice. """
    dic, data = ng.pipe.read(NMRPIPE_3D_FREQ % (1))
    assert data.shape == (3, 8)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 1.
    assert data[1, 0] == 1.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [152.33, 45.67])
    check_ppm_limits(dic, data, 1, [54.70, -32.80])


def test_3d_time_lowmem():
    """ lowmem read/write of 3D NMRPipe time domain file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_3D_TIME)
    assert data.shape == (4, 6, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0].real == 1.
    assert data[0, 0, 0].imag == -1.
    assert data[0, 0, 1].real == 2.
    assert data[0, 0, 1].imag == -2.
    assert data[0, 1, 0].real == 1
    assert data[0, 1, 0].imag == -1
    assert data[1, 0, 0].real == 1
    assert data[1, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, 1, lowmem=True)


def test_3d_freq_lowmem():
    """ lowmem read/write of 3D NMRPipe frequency domain file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_3D_FREQ)
    assert data.shape == (2, 3, 8)
    assert data.dtype == 'float32'
    assert data[0, 0, 0] == 1.
    assert data[0, 0, 1] == 1.
    assert data[0, 1, 0] == 1.
    assert data[1, 0, 0] == 2.
    check_simple_roundtrip(dic, data, 1, lowmem=True)
    check_ppm_limits(dic, data, 0, [220.00, 120.00])
    check_ppm_limits(dic, data, 1, [152.33, 45.67])
    check_ppm_limits(dic, data, 2, [54.70, -32.80])


def test_3d_time_stream():
    """ read/write of 3D NMRPipe time domain stream (single file). """
    dic, data = ng.pipe.read(NMRPIPE_3D_TIME_STREAM)
    assert data.shape == (4, 6, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0].real == 1.
    assert data[0, 0, 0].imag == -1.
    assert data[0, 0, 1].real == 2.
    assert data[0, 0, 1].imag == -2.
    assert data[0, 1, 0].real == 1
    assert data[0, 1, 0].imag == -1
    assert data[1, 0, 0].real == 1
    assert data[1, 0, 0].imag == -1
    check_simple_roundtrip(dic, data)


def test_3d_time_stream_lowmem():
    """ lowmem read/write of 3D NMRPipe time domain stream. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_3D_TIME_STREAM)
    assert data.shape == (4, 6, 8)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0].real == 1.
    assert data[0, 0, 0].imag == -1.
    assert data[0, 0, 1].real == 2.
    assert data[0, 0, 1].imag == -2.
    assert data[0, 1, 0].real == 1
    assert data[0, 1, 0].imag == -1
    assert data[1, 0, 0].real == 1
    assert data[1, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, lowmem=True)


def test_3d_freq_stream():
    """ read/write of 3D NMRPipe frequency domain stream. """
    dic, data = ng.pipe.read(NMRPIPE_3D_FREQ_STREAM)
    assert data.shape == (2, 3, 8)
    assert data.dtype == 'float32'
    assert data[0, 0, 0] == 1.
    assert data[0, 0, 1] == 1.
    assert data[0, 1, 0] == 1.
    assert data[1, 0, 0] == 2.
    check_simple_roundtrip(dic, data)
    check_ppm_limits(dic, data, 0, [220.00, 120.00])
    check_ppm_limits(dic, data, 1, [152.33, 45.67])
    check_ppm_limits(dic, data, 2, [54.70, -32.80])


def test_3d_freq_stream_lowmem():
    """ lowmem read/write of 3D NMRPipe frequency domain stream. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_3D_FREQ_STREAM)
    assert data.shape == (2, 3, 8)
    assert data.dtype == 'float32'
    assert data[0, 0, 0] == 1.
    assert data[0, 0, 1] == 1.
    assert data[0, 1, 0] == 1.
    assert data[1, 0, 0] == 2.
    check_simple_roundtrip(dic, data, lowmem=True)
    check_ppm_limits(dic, data, 0, [220.00, 120.00])
    check_ppm_limits(dic, data, 1, [152.33, 45.67])
    check_ppm_limits(dic, data, 2, [54.70, -32.80])


# All 4D tests are lowmem reads and writes

def test_4d_time_double_index():
    """ read/write of 4D NMRPipe time domain file (two index). """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_TIME_2)
    assert data.shape == (4, 6, 8, 5)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0, 0].real == 1.
    assert data[0, 0, 0, 0].imag == -1.
    assert data[0, 0, 0, 1].real == 2.
    assert data[0, 0, 0, 1].imag == -2.
    assert data[0, 0, 1, 0].real == 1
    assert data[0, 0, 1, 0].imag == -1
    assert data[0, 1, 0, 0].real == 1
    assert data[0, 1, 0, 0].imag == -1
    assert data[1, 0, 0, 0].real == 1
    assert data[1, 0, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, 2, lowmem=True)


def test_4d_time_single_index():
    """ read/write of 4D NMRPipe time domain file (single index). """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_TIME_1)
    assert data.shape == (4, 6, 8, 5)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0, 0].real == 1.
    assert data[0, 0, 0, 0].imag == -1.
    assert data[0, 0, 0, 1].real == 2.
    assert data[0, 0, 0, 1].imag == -2.
    assert data[0, 0, 1, 0].real == 1
    assert data[0, 0, 1, 0].imag == -1
    assert data[0, 1, 0, 0].real == 1
    assert data[0, 1, 0, 0].imag == -1
    assert data[1, 0, 0, 0].real == 1
    assert data[1, 0, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, 1, lowmem=True)


def test_4d_time_stream():
    """ read/write of 4D NMRPipe stream file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_TIME_STREAM)
    assert data.shape == (4, 6, 8, 5)
    assert data.dtype == 'complex64'
    assert data[0, 0, 0, 0].real == 1.
    assert data[0, 0, 0, 0].imag == -1.
    assert data[0, 0, 0, 1].real == 2.
    assert data[0, 0, 0, 1].imag == -2.
    assert data[0, 0, 1, 0].real == 1
    assert data[0, 0, 1, 0].imag == -1
    assert data[0, 1, 0, 0].real == 1
    assert data[0, 1, 0, 0].imag == -1
    assert data[1, 0, 0, 0].real == 1
    assert data[1, 0, 0, 0].imag == -1
    check_simple_roundtrip(dic, data, lowmem=True)


def test_4d_time_double_index_slice():
    """ read/write of 4D NMRPipe time domain (two index) first slice. """
    dic, data = ng.pipe.read(NMRPIPE_4D_TIME_2 % (1, 1))
    assert data.shape == (8, 5)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == -1.
    assert data[0, 1].real == 2.
    assert data[0, 1].imag == -2.
    assert data[1, 0].real == 1
    assert data[1, 0].imag == -1
    check_simple_roundtrip(dic, data)


def test_4d_time_single_index_slice():
    """ read/write of 4D NMRPipe time domain (single index) first slice. """
    dic, data = ng.pipe.read(NMRPIPE_4D_TIME_1 % (1))
    assert data.shape == (8, 5)
    assert data.dtype == 'complex64'
    assert data[0, 0].real == 1.
    assert data[0, 0].imag == -1.
    assert data[0, 1].real == 2.
    assert data[0, 1].imag == -2.
    assert data[1, 0].real == 1
    assert data[1, 0].imag == -1
    check_simple_roundtrip(dic, data)


def test_4d_freq_single_index():
    """ read/write of 4D NMRPipe frequency domain file (single index). """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_FREQ_1)
    assert data.shape == (2, 3, 4, 5)
    assert data.dtype == 'float32'
    assert data[0, 0, 0, 0] == 1.
    assert data[0, 0, 0, 1] == 1.
    assert data[0, 0, 1, 0] == 1.
    assert data[0, 1, 0, 0] == 1.
    assert data[1, 0, 0, 0] == 2.
    check_ppm_limits(dic, data, 0, [180.00, 80.00])
    check_ppm_limits(dic, data, 1, [186.67, 53.33])
    check_ppm_limits(dic, data, 2, [179.00, 59.00])
    check_ppm_limits(dic, data, 3, [44.70, -35.30])
    check_simple_roundtrip(dic, data, 1, lowmem=True)


def test_4d_freq_double_index():
    """ read/write of 4D NMRPipe frequency domain file (double index). """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_FREQ_2)
    assert data.shape == (2, 3, 4, 5)
    assert data.dtype == 'float32'
    assert data[0, 0, 0, 0] == 1.
    assert data[0, 0, 0, 1] == 1.
    assert data[0, 0, 1, 0] == 1.
    assert data[0, 1, 0, 0] == 1.
    assert data[1, 0, 0, 0] == 2.
    check_ppm_limits(dic, data, 0, [180.00, 80.00])
    check_ppm_limits(dic, data, 1, [186.67, 53.33])
    check_ppm_limits(dic, data, 2, [179.00, 59.00])
    check_ppm_limits(dic, data, 3, [44.70, -35.30])
    check_simple_roundtrip(dic, data, 2, lowmem=True)


def test_4d_freq_stream():
    """ read/write of 4D NMRPipe frequency domain stream file. """
    dic, data = ng.pipe.read_lowmem(NMRPIPE_4D_FREQ_STREAM)
    assert data.shape == (2, 3, 4, 5)
    assert data.dtype == 'float32'
    assert data[0, 0, 0, 0] == 1.
    assert data[0, 0, 0, 1] == 1.
    assert data[0, 0, 1, 0] == 1.
    assert data[0, 1, 0, 0] == 1.
    assert data[1, 0, 0, 0] == 2.
    check_ppm_limits(dic, data, 0, [180.00, 80.00])
    check_ppm_limits(dic, data, 1, [186.67, 53.33])
    check_ppm_limits(dic, data, 2, [179.00, 59.00])
    check_ppm_limits(dic, data, 3, [44.70, -35.30])
    check_simple_roundtrip(dic, data, lowmem=True)


def test_4d_freq_double_index_slice():
    """ read/write of 4D NMRPipe freq domain (two index) first slice. """
    dic, data = ng.pipe.read(NMRPIPE_4D_FREQ_2 % (1, 1))
    assert data.shape == (4, 5)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 1.
    assert data[1, 0] == 1.
    check_ppm_limits(dic, data, 0, [179.00, 59.00])
    check_ppm_limits(dic, data, 1, [44.70, -35.30])
    check_simple_roundtrip(dic, data)


def test_4d_freq_single_index_slice():
    """ read/write of 4D NMRPipe freq domain (single index) first slice. """
    dic, data = ng.pipe.read(NMRPIPE_4D_FREQ_1 % (1))
    assert data.shape == (4, 5)
    assert data.dtype == 'float32'
    assert data[0, 0] == 1.
    assert data[0, 1] == 1.
    assert data[1, 0] == 1.
    check_ppm_limits(dic, data, 0, [179.00, 59.00])
    check_ppm_limits(dic, data, 1, [44.70, -35.30])
    check_simple_roundtrip(dic, data)


def test_make_uc():
    """ ng.pipe.make_uc """
    dic, data = ng.pipe.read(NMRPIPE_2D_FREQ)

    # H1 dimension
    uc = ng.pipe.make_uc(dic, data)
    assert uc._size == 8
    assert uc._cplx is False
    assert uc._sw == 50000.0
    assert uc._obs == 500.0
    assert round(uc._car) == 2350.

    # H1 dimension
    uc = ng.pipe.make_uc(dic, data, 1)
    assert uc._size == 8
    assert uc._cplx is False
    assert uc._sw == 50000.0
    assert uc._obs == 500.0
    assert round(uc._car) == 2350.

    # C13 dimension
    uc = ng.pipe.make_uc(dic, data, 0)
    assert uc._size == 2
    assert uc._cplx is False
    assert uc._sw == 20000.0
    assert uc._obs == 125.0
    assert round(uc._car) == 12375.


def test_guess_udic():
    """ ng.pipe.guess_udic """
    dic, data = ng.pipe.read(NMRPIPE_2D_FREQ)
    udic = ng.pipe.guess_udic(dic, data)
    assert udic[0]['car'] == 12375.0
    assert udic[0]['complex'] is False
    assert udic[0]['encoding'] == 'states'
    assert udic[0]['freq'] is True
    assert udic[0]['label'] == 'C13'
    assert udic[0]['obs'] == 125.0
    assert udic[0]['size'] == 2
    assert udic[0]['sw'] == 20000.0
    assert udic[0]['time'] is False
    assert round(udic[1]['car']) == 2350
    assert udic[1]['complex'] is False
    assert udic[1]['encoding'] == 'states'
    assert udic[1]['freq'] is True
    assert udic[1]['label'] == 'H1'
    assert udic[1]['obs'] == 500.0
    assert udic[1]['size'] == 8
    assert udic[1]['sw'] == 50000.0
    assert udic[1]['time'] is False
    assert udic['ndim'] == 2

    dic, data = ng.pipe.read(NMRPIPE_2D_TIME)
    udic = ng.pipe.guess_udic(dic, data)
    assert udic[0]['car'] == 12375.0
    assert udic[0]['complex'] is True
    assert udic[0]['encoding'] == 'states'
    assert udic[0]['freq'] is False
    assert udic[0]['label'] == 'C13'
    assert udic[0]['obs'] == 125.0
    assert udic[0]['size'] == 4
    assert udic[0]['sw'] == 20000.0
    assert udic[0]['time'] is True
    assert round(udic[1]['car']) == 2350
    assert udic[1]['complex'] is True
    assert udic[1]['encoding'] == 'states'
    assert udic[1]['freq'] is False
    assert udic[1]['label'] == 'H1'
    assert udic[1]['obs'] == 500.0
    assert udic[1]['size'] == 8
    assert udic[1]['sw'] == 50000.0
    assert udic[1]['time'] is True
    assert udic['ndim'] == 2

def test_read_table():
    comments, tbl_format, tbl = ng.pipe.read_table(NMRPIPE_TABLE)

    assert len(comments) == 13
    assert len(tbl_format) == 25
    assert len(tbl) == 5
    assert len(tbl[0]) == 25

    assert tbl_format[0] == '%5d'
    assert tbl_format[1] == '%9.3f'

    assert tbl['INDEX'][0] == 1
    assert tbl['INDEX'][-1] == 5

    assert tbl['X1'][0] == 1328
    assert tbl['X1'][-1] == 1161

def test_write_table():
    ref_comments, ref_tbl_format, ref_tbl = ng.pipe.read_table(NMRPIPE_TABLE)
    try:
        tbl_fname = tempfile.mktemp(dir='.')
        ng.pipe.write_table(tbl_fname, ref_comments, ref_tbl_format, ref_tbl)
        comments, tbl_format, tbl = ng.pipe.read_table(tbl_fname)
        assert all(ref_tbl == tbl)
        assert ref_comments == comments
        assert ref_tbl_format == tbl_format
    finally:
        os.remove(tbl_fname)
