"""Tests for the fileio.jeol submodule"""

import os

import numpy as np
import nmrglue as ng
import pytest

from setup import DATA_DIR
JEOLDATA = os.path.join(DATA_DIR, "jeol")


def test_1d_real():
    """data with a real dimension"""

    pass


def test_1d_complex_1():
    """data with a complex dimension"""

    basename = "cyclosporine_Carbon-1-1"
    jdic, jdata = ng.jeol.read(os.path.join(JEOLDATA, f"{basename}.jdf"))
    pdic, pdata = ng.pipe.read(os.path.join(JEOLDATA, f"{basename}.fid"))
    assert np.allclose(jdata, pdata, rtol=1e-7)

def test_1d_complex_1_udic():
    basename = "cyclosporine_Carbon-1-1"
    dic, data = ng.jeol.read(os.path.join(JEOLDATA, f"{basename}.jdf"))
    udic = ng.jeol.guess_udic(dic, data)
    assert udic["ndim"] == 1
    assert udic[0]["sw"] - 39556.962025316454 < 0.0001
    assert udic[0]["obs"] - 125.76529768332652 < 0.0001
    assert udic[0]["car"] - 12576.529768332653 < 0.0001
    assert udic[0]["size"] == 65536
    assert udic[0]["label"] == "Carbon13"
    assert udic[0]["encoding"] == "complex"

def test_1d_complex_2():
    """data with a complex dimension"""

    basename = "cyclosporine_Proton-1-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_1d_complex_3():
    """data with a complex dimension"""

    basename = "2mM sucrose_WATERGATE-3-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_1d_complex_4():
    """data with a complex dimension"""

    basename = "2-ethyl-1-indanone 1_PROTON-21-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_2d_rr():
    """data with two real dimensions"""

    pass

@pytest.mark.slow
def test_2d_cc_1():
    """data with two complex dimensions"""

    basename = "cyclosporine_tocsy-1-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)

@pytest.mark.slow
def test_2d_cc_1_udic():
    basename = "cyclosporine_tocsy-1-1"
    dic, data = ng.jeol.read(os.path.join(JEOLDATA, f"{basename}.jdf"))
    udic = ng.jeol.guess_udic(dic, data)
    assert udic["ndim"] == 2
    assert udic[0]["sw"] - 6262.5250501002 < 1e-4
    assert udic[0]["obs"] - 500.15991520961256 < 1e-4
    assert udic[0]["car"] - 2500.7995760480626 < 1e-4
    assert udic[0]["size"] == 2048
    assert udic[0]["label"] == "Proton"
    assert udic[0]["encoding"] == "complex"
    assert udic[1]["sw"] - 5002.000800320128 < 1e-4
    assert udic[1]["obs"] - 500.15991520961256 < 1e-4
    assert udic[1]["car"] - 2500.7995760480626 < 1e-4
    assert udic[1]["size"] == 4096
    assert udic[1]["label"] == "Proton"
    assert udic[1]["encoding"] == "complex"

@pytest.mark.slow
def test_2d_cc_2():
    """data with two complex dimensions"""

    basename = "cyclosporine_HSQC-1-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_2d_cc_nus():
    """NUS data with two complex dimensions"""

    basename = "2-ethyl-1-indanone 1_HSQC_NUS-3-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_2d_rc_nus():
    """NUS data with two real_complex dimensions"""
    basename = "2-ethyl-1-indanone 1_HMBC_NUS-2-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)

@pytest.mark.slow
def test_2d_rc():
    """data with two real_complex dimensions"""

    basename = "cyclosporine_dqf_cosy_pfg-1-1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)


def test_2d_cr():
    """data with one complex dimension and one real dimension (example: 1d optimization)"""

    basename = "jeol-2d_1"
    jdic, jdata = ng.jeol.read(f"{JEOLDATA}/{basename}.jdf")
    pdic, pdata = ng.pipe.read(f"{JEOLDATA}/{basename}.fid")
    assert np.allclose(jdata, pdata, rtol=1e-7)