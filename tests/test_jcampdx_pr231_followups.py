"""Self-contained JCAMP-DX follow-up tests for PR #231."""

import numpy as np
import pytest

import nmrglue as ng


@pytest.mark.parametrize("components", ["RI", "R", "I"])
def test_read_as_complex(tmp_path, components):
    content = (
        "##TITLE=Synthetic FID\n"
        "##JCAMP-DX=6.0\n"
        "##DATA TYPE=NMR FID\n"
        "##DATA CLASS=NTUPLES\n"
        "##NTUPLES=NMR FID\n"
        "##SYMBOL=X,R,I\n"
        "##FACTOR=1,1,1\n"
        "##FIRST=0,1,5\n"
        "##LAST=0.75,4,8\n"
        "##UNITS=SECONDS,ARBITRARY UNITS,ARBITRARY UNITS\n"
    )
    if "R" in components:
        content += "##DATA TABLE=(X++(R..R)),XYDATA\n0 1 2 3 4\n"
    if "I" in components:
        content += "##DATA TABLE=(X++(I..I)),XYDATA\n0 5 6 7 8\n"
    if components != "RI":
        content = content.replace("##FACTOR=1,1,1\n", "")
    content += "##END=\n"
    path = tmp_path / "fid.dx"
    path.write_text(content)

    dic, default = ng.jcampdx.read(path)
    explicit_dic, explicit_default = ng.jcampdx.read(path, as_complex=False)
    complex_dic, converted = ng.jcampdx.read(path, as_complex=True)
    assert dic == explicit_dic == complex_dic
    if components == "RI":
        assert isinstance(default, list)
        for expected, actual in zip(default, explicit_default):
            np.testing.assert_array_equal(actual, expected)
        assert converted.dtype == np.complex128
        np.testing.assert_array_equal(converted, [1+5j, 2+6j, 3+7j, 4+8j])
        assert ng.jcampdx.guess_udic(complex_dic, converted)[0]["complex"]
    elif components == "R":
        np.testing.assert_array_equal(default, [1, 2, 3, 4])
        np.testing.assert_array_equal(explicit_default, default)
        np.testing.assert_array_equal(converted, default)
    else:
        assert default[0] is explicit_default[0] is converted[0] is None
        np.testing.assert_array_equal(default[1], [5, 6, 7, 8])
        np.testing.assert_array_equal(explicit_default[1], default[1])
        np.testing.assert_array_equal(converted[1], default[1])


def test_read_as_complex_spectrum(tmp_path):
    path = tmp_path / "spectrum.dx"
    path.write_text(
        "##TITLE=Synthetic spectrum\n"
        "##JCAMP-DX=5.0\n"
        "##DATA TYPE=NMR SPECTRUM\n"
        "##XYDATA=(X++(Y..Y))\n0 10 20 30\n"
        "##END=\n"
    )
    dic, default = ng.jcampdx.read(path)
    converted_dic, converted = ng.jcampdx.read(path, as_complex=True)
    assert dic == converted_dic
    assert converted.dtype == default.dtype
    np.testing.assert_array_equal(converted, [10, 20, 30])
    np.testing.assert_array_equal(converted, default)


@pytest.mark.parametrize("components", ["RI", "R", "I"])
def test_guess_udic_available_component(components):
    dic = {
        "DATATYPE": ["NMR FID"],
        "FIRSTX": ["0"],
        "LASTX": ["0.75"],
        "XUNITS": ["SECONDS"],
    }
    data = [np.arange(4.) if "R" in components else None,
            np.arange(4.) if "I" in components else None]
    udic = ng.jcampdx.guess_udic(dic, data)[0]
    assert udic["size"] == 4
    assert udic["sw"] == pytest.approx(4.)
    assert udic["time"] is True
    assert udic["complex"] is False


@pytest.mark.parametrize("data", [None, [], [None, None]])
def test_guess_udic_no_available_component(data):
    dic = {"FIRSTX": ["0"], "LASTX": ["1"], "XUNITS": ["HZ"]}
    with pytest.warns(UserWarning, match="No data, cannot set udic size"):
        udic = ng.jcampdx.guess_udic(dic, data)
    assert udic[0]["size"] == ng.fileiobase.create_blank_udic(1)[0]["size"]


@pytest.mark.parametrize("datatype, ntuples, is_fid", [
    (None, "NMR FID", True),
    (None, " nmr fid ", True),
    (None, "NMR SPECTRUM", False),
    (None, None, False),
    ("NMR SPECTRUM", "NMR FID", False),
    ("NMR FID", "NMR SPECTRUM", True),
])
def test_guess_udic_fid_detection(datatype, ntuples, is_fid):
    dic = {"FIRSTX": ["0"], "LASTX": ["0.75"], "XUNITS": ["SECONDS"]}
    if datatype is not None:
        dic["DATATYPE"] = [datatype]
    if ntuples is not None:
        dic["NTUPLES"] = [ntuples]
    udic = ng.jcampdx.guess_udic(dic, np.arange(4.))[0]
    assert udic["time"] is is_fid
    assert udic["freq"] is (not is_fid)
    assert udic["sw"] == pytest.approx(4. if is_fid else .75)


def test_read_fid_without_datatype(tmp_path):
    path = tmp_path / "fid.dx"
    path.write_text(
        "##TITLE=Synthetic FID without DATA TYPE\n"
        "##JCAMP-DX=6.0\n"
        "##DATA CLASS=NTUPLES\n"
        "##NTUPLES=NMR FID\n"
        "##SYMBOL=X,R,I\n"
        "##FACTOR=1,1,1\n"
        "##FIRST=0,1,5\n"
        "##LAST=0.75,4,8\n"
        "##UNITS=SECONDS,ARBITRARY UNITS,ARBITRARY UNITS\n"
        "##DATA TABLE=(X++(R..R)),XYDATA\n0 1 2 3 4\n"
        "##DATA TABLE=(X++(I..I)),XYDATA\n0 5 6 7 8\n"
        "##END=\n"
    )
    dic, data = ng.jcampdx.read(path)
    assert "DATATYPE" not in dic
    udic = ng.jcampdx.guess_udic(dic, data)[0]
    assert udic["time"] is True
    assert udic["freq"] is False
    assert udic["size"] == 4
    assert udic["sw"] == pytest.approx(4.)
