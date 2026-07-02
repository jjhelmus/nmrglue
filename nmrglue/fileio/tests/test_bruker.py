""" Unit tests for nmrglue/fileio/bruker.py module """

import os
import warnings

import numpy as np
import pytest
import nmrglue as ng
import shutil
import tempfile

# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'bruker_test_data')

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


def test_read_jcamp_undecodable_bytes():
    """bytes invalid in both utf-8 and cp1252 do not crash the reader"""
    # 0x81 is undefined in cp1252 and invalid utf-8; latin-1 fallback
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
