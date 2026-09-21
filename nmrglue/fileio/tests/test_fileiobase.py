""" Unit tests for nmrglue/fileio/fileiobase.py module """

import os

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import nmrglue as ng
import pytest


# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
NMRPIPE_1D_FREQ = os.path.join(DATA_DIR, 'nmrpipe_1d_freq.fid')


def test_uc_from_freqscale():
    """
    Test that `uc_from_freqscale` gives equivalent results as `uc_from_udic`.
    """
    from nmrglue.fileio.fileiobase import uc_from_freqscale

    # read frequency test data
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


def test_uc_with_reference_frequency():
    """A reference frequency is used for ppm and Hz conversions."""
    uc = ng.fileiobase.unit_conversion(
        size=4, cplx=False, sw=400.0, obs=100.0, car=100.0, ref=80.0)

    assert_allclose(uc.ppm_scale(), [3.75, 2.5, 1.25, 0.0])
    assert_allclose(uc.hz_scale(), [300.0, 200.0, 100.0, 0.0])
    assert_allclose(uc.f(uc.ppm_scale(), 'ppm'), np.arange(4))


def test_uc_reference_defaults_to_observation_frequency():
    """Omitting ref preserves the historical unit conversion behavior."""
    legacy = ng.fileiobase.unit_conversion(
        size=4, cplx=False, sw=400.0, obs=100.0, car=100.0)
    explicit = ng.fileiobase.unit_conversion(
        size=4, cplx=False, sw=400.0, obs=100.0, car=100.0, ref=100.0)

    assert_array_equal(legacy.ppm_scale(), explicit.ppm_scale())
    assert_array_equal(legacy.hz_scale(), explicit.hz_scale())


def test_reference_frequency_is_used_by_exporters():
    """Formats with one frequency field preserve the referenced ppm scale."""
    udic = ng.fileiobase.create_blank_udic(1)
    udic[0].update(size=8, obs=100.0, ref=80.0, car=100.0)

    pipe_dic = ng.pipe.create_dic(udic)
    sparky_dic = ng.sparky.create_dic(udic)
    rnmrtk_dic = ng.rnmrtk.create_dic(udic)

    assert pipe_dic['FDF2OBS'] == 80.0
    assert pipe_dic['FDF2CAR'] == 1.25
    assert sparky_dic['w1']['spectrometer_freq'] == 80.0
    assert sparky_dic['w1']['xmtr_freq'] == 1.25
    assert rnmrtk_dic['sf'] == [80.0]
    assert rnmrtk_dic['ppm'] == [1.25]


def test_glue_reference_frequency_roundtrip(tmp_path):
    """Glue round-trips an optional reference frequency."""
    pytest.importorskip('h5py')
    from nmrglue.fileio import glue

    data = np.zeros(8)
    referenced = ng.fileiobase.create_blank_udic(1)
    referenced[0].update(size=8, obs=100.0, ref=80.0, car=100.0)
    legacy = ng.fileiobase.create_blank_udic(1)
    legacy[0].update(size=8, obs=100.0, car=100.0)

    for name, udic in [('referenced.glue', referenced),
                       ('legacy.glue', legacy)]:
        path = tmp_path / name
        glue.write(path, udic, data)
        read_dic, _ = glue.read(path)
        assert read_dic[0].get('ref') == udic[0].get('ref')


def test_spinsolve_extracts_reference_frequency():
    """The Spinsolve $SF parameter is retained as a reference frequency."""
    params = {
        '$BF1': ['43.5'],
        '$SF': ['43.4998'],
        '$SW': ['100.0'],
        '.OBSERVENUCLEUS': ['^1H'],
    }

    adic = ng.spinsolve.get_udic_from_jcamp_dict(params)

    assert adic['obs'] == 43.5
    assert adic['ref'] == 43.4998
    assert_allclose(adic['car'], 200.0)


def test_update_uc():
    uc = ng.fileiobase.unit_conversion(
        size=64.0, cplx=False, sw=1.0, obs=1.0, car=1.0)
    uc2 = ng.fileiobase.update_uc(uc, size=10, cplx=True, sw=2.0, car=5.2, obs=3.0)
    assert uc2._size == 10
    assert uc2._cplx is True
    assert abs(uc2._sw - 2.0) < 1e-5
    assert abs(uc2._car - 5.2) < 1e-5
    assert abs(uc2._obs - 3.0) < 1e-5

    referenced = ng.fileiobase.unit_conversion(
        size=64, cplx=False, sw=1.0, obs=1.0, car=1.0, ref=0.8)
    updated = ng.fileiobase.update_uc(referenced, obs=2.0)
    assert updated._obs == 2.0
    assert updated._ref == 0.8
