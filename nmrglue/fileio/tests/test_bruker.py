""" Unit tests for nmrglue/fileio/bruker.py module """

import os
import shutil
import locale
import tempfile
import warnings

import numpy as np
import nmrglue as ng
import pytest

# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'bruker_test_data')


def _make_2d_dic(indirect_dic):
    """Return minimal Bruker parameters for a 2D dataset."""
    return {
        'acqus': {
            'AQ_mod': 3,
            'NUC1': '1H',
            'O1': 0.0,
            'SFO1': 400.0,
            'SW_h': 10000.0,
        },
        'procs': {'SF': 400.0},
        **indirect_dic,
    }


@pytest.mark.parametrize('fnmode', [1, 2])
def test_guess_udic_magnitude_indirect_dimension(fnmode):
    """QF and QSEC indirect dimensions are real, not complex."""
    dic = _make_2d_dic({
        'acqu2s': {
            'FnMODE': fnmode,
            'NUC1': '1H',
            'SFO1': 400.0,
            'SW': 10.0,
        },
        'proc2s': {'SF': 400.0},
    })
    data = np.zeros((40, 128), dtype=np.complex128)

    udic = ng.bruker.guess_udic(dic, data)
    converter = ng.convert.converter()
    converter.from_bruker(dic, data)
    pipe_dic, pipe_data = converter.to_pipe()
    _, transposed = ng.pipe_proc.tp(pipe_dic.copy(), pipe_data)

    assert udic[0]['encoding'] == 'magnitude'
    assert udic[0]['complex'] is False
    assert pipe_dic['FDF2QUADFLAG'] == 0.0
    assert pipe_dic['FDF1QUADFLAG'] == 1.0
    assert transposed.shape == (128, 40)


@pytest.mark.parametrize('mc2', [0, 1])
def test_guess_udic_magnitude_indirect_dimension_from_proc_params(mc2):
    """Processed QF and QSEC indirect dimensions are real, not complex."""
    dic = _make_2d_dic({
        'proc2s': {
            'AXNUC': '1H',
            'MC2': mc2,
            'OFFSET': 0.0,
            'SF': 400.0,
            'SW_p': 4000.0,
        },
    })
    data = np.zeros((40, 128), dtype=np.complex128)

    udic = ng.bruker.guess_udic(dic, data)

    assert udic[0]['encoding'] == 'magnitude'
    assert udic[0]['complex'] is False

def test_read_pdata():
    """Reading processed bruker 1D data"""

    # read processed data
    specdir = os.path.join(DATA_DIR, '1', 'pdata', '1')
    dic, data = ng.bruker.read_pdata(specdir)

    # check that the data is correct
    assert np.all(data == [-36275840.0, -34775104.0])

    # check dictionaries are correct
    assert len(dic.keys()) == 2
    assert 'procs' in dic.keys()
    assert 'acqus' in dic.keys()
    assert len(dic['procs'].keys()) == 131
    assert len(dic['acqus'].keys()) == 298

    # specifying proc(s) files
    # (fired this error in 0.8 version:
    # UnboundLocalError: local variable 'pdata_path' referenced before
    # assignment)

    proc = os.path.join(specdir, 'procs')
    dic, data = ng.bruker.read_pdata(specdir, procs_files = [proc])
    assert len(dic['procs'].keys()) == 131


def test_reorder_submatrix():
    """reordering submatrix back and forth"""

    # make a dummy matrix
    data = np.arange(16, dtype='float64').reshape(4, 4)

    # correctly reordered matrix
    rdata = np.array([[ 0,  1,  4,  5],
                      [ 2,  3,  6,  7],
                      [ 8,  9, 12, 13],
                      [10, 11, 14, 15]], dtype='float64')

    # reorder from the submatrix form
    r1data = ng.bruker.reorder_submatrix(data, shape=(4, 4),
                                         submatrix_shape=(2, 2), reverse=False)

    # reorder to the submatrix form
    r2data = ng.bruker.reorder_submatrix(r1data, shape=(4, 4),
                                         submatrix_shape=(2, 2), reverse=True)
    # checks
    assert np.all(rdata == r1data)
    assert np.all(data == r2data)


def test_write_pdata():
    """ Writing a processed Bruker dataset """

    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, '1', 'pdata', '1'))

    # write to a temperory file
    td = tempfile.mkdtemp('.')
    ng.bruker.write_pdata(td, dic, data, write_procs=True, pdata_folder=10)

    assert os.path.isdir(os.path.join(td, 'pdata', '10'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', 'procs'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', 'proc'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', '1r'))

    rdic, rdata = ng.bruker.read_pdata(os.path.join(td, 'pdata', '10'))

    assert np.all(data == rdata)
    assert rdic['procs'].keys() == dic['procs'].keys()
    shutil.rmtree(td)


def _real_acqus_bytes():
    """Bytes of the real acqus file shipped with the test data."""
    with open(os.path.join(DATA_DIR, '1', 'acqus'), 'rb') as f:
        return f.read()


def _write_temp(content):
    fd, temp_path = tempfile.mkstemp()
    with os.fdopen(fd, 'wb') as f:
        f.write(content)
    return temp_path


def test_read_jcamp_cp1252():
    """cp1252-encoded acqus (e.g. degree sign) decodes correctly"""
    # real acqus with a realistic cp1252 (non utf-8) parameter added,
    # as written by instruments configured with a western locale
    content = _real_acqus_bytes().replace(
        b"##END=", "##$SOLVENT= <CDCl3 at 25\u00b0C>\n##END=".encode("cp1252"))
    temp_path = _write_temp(content)
    try:
        dic = ng.bruker.read_jcamp(temp_path)
        # correct character, no U+FFFD corruption
        assert dic["SOLVENT"] == "CDCl3 at 25\u00b0C"
        # remainder of the real file still parsed
        assert dic["LOCKED"] is True
    finally:
        os.remove(temp_path)


def test_read_jcamp_undecodable_bytes(monkeypatch):
    """bytes invalid in both utf-8 and cp1252 do not crash the reader"""
    # 0x81 is undefined in cp1252 and invalid utf-8; latin-1 fallback.
    # A utf-8 locale adds nothing, so the result does not depend on the
    # machine running the test.
    monkeypatch.setattr(locale, "getpreferredencoding", lambda *a: "UTF-8")
    content = _real_acqus_bytes().replace(
        b"##END=", b"##$BAD= <\x81>\n##END=")
    temp_path = _write_temp(content)
    try:
        with pytest.warns(UserWarning, match="latin-1"):
            dic = ng.bruker.read_jcamp(temp_path)
        assert dic["BAD"] == "\x81"  # latin-1 maps byte to same codepoint
        assert dic["LOCKED"] is True  # rest of file intact
    finally:
        os.remove(temp_path)


def test_read_jcamp_explicit_encoding():
    """an explicit encoding is tried first"""
    # 0xb0 is a degree sign in cp1252 but an infinity sign in mac-roman;
    # with an explicit encoding the caller's choice must win
    content = _real_acqus_bytes().replace(
        b"##END=", b"##$TEMPUNIT= <\xb0C>\n##END=")
    temp_path = _write_temp(content)
    try:
        dic = ng.bruker.read_jcamp(temp_path, encoding="mac-roman")
        assert dic["TEMPUNIT"] == "\u221eC"
        dic = ng.bruker.read_jcamp(temp_path)
        assert dic["TEMPUNIT"] == "\u00b0C"
    finally:
        os.remove(temp_path)


def test_read_jcamp_utf8_bom():
    """a byte order mark does not hide the first record"""
    content = b"\xef\xbb\xbf" + _real_acqus_bytes()
    temp_path = _write_temp(content)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dic = ng.bruker.read_jcamp(temp_path)
        # left in place, the BOM makes ##TITLE= unrecognisable and it is
        # discarded as an extraneous line
        assert not [w for w in caught if "Extraneous line" in str(w.message)]
        assert dic["_coreheader"][0].startswith("##TITLE=")
        assert dic["LOCKED"] is True  # remainder of the real file still parsed
    finally:
        os.remove(temp_path)


def test_read_jcamp_line_endings():
    """CRLF and CR-only line endings read as LF does"""
    lf = _real_acqus_bytes().replace(b"\r\n", b"\n")
    expected = None
    for eol in (b"\n", b"\r\n", b"\r"):
        temp_path = _write_temp(lf.replace(b"\n", eol))
        try:
            dic = ng.bruker.read_jcamp(temp_path)
        finally:
            os.remove(temp_path)
        if expected is None:
            expected = dic
        assert dic == expected
    assert expected["LOCKED"] is True


def test_read_jcamp_locale_encoding(monkeypatch):
    """the locale's encoding is tried before falling back to latin-1"""
    # 0x81 fails utf-8 and cp1252 but is a letter in cp1251
    content = _real_acqus_bytes().replace(
        b"##END=", b"##$NAME= <\x81>\n##END=")
    temp_path = _write_temp(content)
    try:
        monkeypatch.setattr(locale, "getpreferredencoding",
                            lambda *a: "cp1251")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            dic = ng.bruker.read_jcamp(temp_path)
        assert dic["NAME"] == "\u0403"
        assert dic["LOCKED"] is True

        # a latin-1 locale is a choice, not a fallback: no warning
        monkeypatch.setattr(locale, "getpreferredencoding",
                            lambda *a: "ISO-8859-1")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            dic = ng.bruker.read_jcamp(temp_path)
        assert dic["NAME"] == "\x81"

        # a locale encoding Python does not know is skipped
        monkeypatch.setattr(locale, "getpreferredencoding",
                            lambda *a: "no-such-codec")
        with pytest.warns(UserWarning, match="latin-1"):
            dic = ng.bruker.read_jcamp(temp_path)
        assert dic["NAME"] == "\x81"
    finally:
        os.remove(temp_path)
