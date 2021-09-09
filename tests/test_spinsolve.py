""" Tests for the fileio.spinsolve submodule """


import nmrglue as ng
from pathlib import Path

from setup import DATA_DIR


def test_acqu():
    """ read nmr_fid.dx """
    dic, data = ng.spinsolve.read(Path(DATA_DIR) / "spinsolve" / "ethanol", "nmr_fid.dx")
    assert dic["acqu"]["Sample"] == "EtOH"
    assert dic["acqu"]["Solvent"] == "None"


def test_jcamp_dx():
    """ read nmr_fid.dx """
    dic, data = ng.spinsolve.read(Path(DATA_DIR) / "spinsolve" / "ethanol", "nmr_fid.dx")
    assert data.size == 32768
    assert data.shape == (32768,)
    assert "Magritek Spinsolve" in dic["dx"]["_comments"][0]


def test_data1d():
    """ read nmr_fid.dx """
    dic, data = ng.spinsolve.read(Path(DATA_DIR) / "spinsolve" / "ethanol", "data.1d")
    assert dic["spectrum"]["xDim"] == 32768
    assert len(dic["spectrum"]["xaxis"]) == 32768
    assert data.size == 32768
    assert data.shape == (32768,)


def test_guess_acqu():
    """ guess_udic based on acqu dictionary """
    dic, data = ng.spinsolve.read(Path(DATA_DIR) / "spinsolve" / "ethanol", "nmr_fid.dx")
    udic = ng.spinsolve.guess_udic(dic, data)
    assert udic[0]["sw"] == 5000
    assert 43.49 < udic[0]["obs"] < 43.50
    assert 206 < udic[0]["car"] < 207
    assert udic[0]["size"] == 32768
    assert udic[0]["label"] == "1H"


def test_guess_jcamp_dx():
    """ guess_udic based on dx dictionary """
    dic, data = ng.spinsolve.read(Path(DATA_DIR) / "spinsolve" / "ethanol", "nmr_fid.dx")

    # Drop acqu dict that would be used as default
    dic["acqu"] = {}

    udic = ng.spinsolve.guess_udic(dic, data)
    assert 4999 < udic[0]["sw"] < 5001
    assert 43.49 < udic[0]["obs"] < 43.50
    assert 206 < udic[0]["car"] < 207
    assert udic[0]["size"] == 32768
    assert udic[0]["label"] == "1H"
