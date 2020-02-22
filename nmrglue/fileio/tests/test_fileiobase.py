""" Unit tests for nmrglue/fileio/fileiobase.py module """

import os

from numpy.testing import assert_array_equal
import nmrglue as ng


# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
NMRPIPE_1D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_1d_freq.fid')


def test_uc_from_freqscale():
    """
    Test that `uc_from_freqscale` gives equivalent results as `uc_from_udic`.
    """
    from nmrglue.fileio.fileiobase import uc_from_freqscale

    # read fequency test data
    dic, data = ng.pipe.read(NMRPIPE_1D_FREQ)

    # make udic and uc using uc_from_udic
    udic = ng.pipe.guess_udic(dic, data)
    uc = ng.fileiobase.uc_from_udic(udic)

    ppm_scale = uc.ppm_scale()
    uc_from_ppm = uc_from_freqscale(ppm_scale, udic[0]['obs'], 'ppm')
    new_ppm_scale = uc_from_ppm.ppm_scale()
    assert_array_equal(ppm_scale, new_ppm_scale)

    hz_scale = uc.hz_scale()
    uc_from_hz = uc_from_freqscale(hz_scale, udic[0]['obs'], 'hz')
    new_hz_scale = uc_from_hz.hz_scale()
    assert_array_equal(hz_scale, new_hz_scale)

    khz_scale = hz_scale * 1.0e3
    uc_from_khz = uc_from_freqscale(khz_scale, udic[0]['obs'], 'khz')
    new_khz_scale = uc_from_khz.hz_scale() * 1.0e3
    assert_array_equal(khz_scale, new_khz_scale)


# regression test for https://github.com/jjhelmus/nmrglue/issues/113
def test_uc_with_float_size():
    uc = ng.fileiobase.unit_conversion(
        size=64.0, cplx=False, sw=1.0, obs=1.0, car=1.0)
    scale = uc.ppm_scale()
    assert len(scale) == 64
