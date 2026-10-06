"""
tests/test_tecmag.py — Tests for nmrglue.fileio.tecmag

Covers two layers:

1. Integration tests — read a real .tnt file from DATA_DIR and verify
   the parsed data against a known reference.

2. Unit tests — synthesise minimal valid .tnt byte sequences to exercise
   read(), guess_tables(), and the internal helpers _read_pascal_pair()
   and _convert_si() without requiring external test data.
"""
import os
import struct

import numpy as np
from numpy.testing import assert_allclose, assert_array_almost_equal
import pytest

import nmrglue as ng
from nmrglue.fileio.tecmag import (
    TNTTMAG, TNTTMG2, guess_tables, _read_pascal_pair, _convert_si,
)

from setup import DATA_DIR


# =============================================================================
# Helpers — synthesise minimal valid .tnt byte sequences
# =============================================================================

def _pascal_string(s):
    """Build a length-prefixed Pascal string (4-byte LE length + bytes)."""
    if isinstance(s, str):
        s = s.encode('ascii')
    return struct.pack('<i', len(s)) + s


def _tlv_section(tag, payload, flag=0):
    """Build a TLV record: 4-byte tag + 4-byte flag + 4-byte length + payload."""
    tag_bytes = tag.encode('ascii') if isinstance(tag, str) else tag
    assert len(tag_bytes) == 4, "Tag must be exactly 4 bytes"
    return tag_bytes + struct.pack('<II', flag, len(payload)) + payload


def _build_minimal_tnt(table_pairs=(), n_points=8, n_records=4):
    """Synthesise a minimal valid .tnt file with TMAG/PSEQ/DATA/TMG2 sections.

    Parameters
    ----------
    table_pairs : iterable of (name, content) str tuples
        Pascal-string pairs to embed in the PSEQ section.
    n_points, n_records : int
        Shape of the DATA section (complex64 buffer).
    """
    magic = b'TNT1.008'

    # TMAG — patch actual_npts and key spectral params; rest zeroed
    tmag_arr = np.frombuffer(bytes(TNTTMAG.itemsize), dtype=TNTTMAG).copy()
    tmag_arr[0]['actual_npts'][0] = n_points
    tmag_arr[0]['actual_npts'][1] = n_records
    tmag_arr[0]['actual_npts'][2] = 1
    tmag_arr[0]['actual_npts'][3] = 1
    tmag_arr[0]['sw'][0]      = 20000.0
    tmag_arr[0]['ob_freq'][0] = 21.31
    tmag_payload = tmag_arr.tobytes()

    # DATA — zeros, complex64
    data_payload = np.zeros(n_points * n_records, dtype='<c8').tobytes()

    # PSEQ — Pascal-string table pairs with padding
    pseq_payload = b'\x00' * 32
    for name, content in table_pairs:
        pseq_payload += _pascal_string(name)
        pseq_payload += _pascal_string(content)
        pseq_payload += b'\x00' * 8

    # TMG2 — zeroed processing parameters
    tmg2_payload = bytes(TNTTMG2.itemsize)

    sections = (
        _tlv_section('TMAG', tmag_payload) +
        _tlv_section('DATA', data_payload) +
        _tlv_section('PSEQ', pseq_payload) +
        _tlv_section('TMG2', tmg2_payload)
    )
    return magic + sections


def _write_tnt(tmp_path, table_pairs=(), n_points=8, n_records=4):
    """Write a synthetic .tnt into pytest's tmp_path and return the path str.

    Using tmp_path (pytest fixture) instead of tempfile.mkstemp avoids
    PermissionError on Windows when np.memmap holds the file open after read().
    pytest deletes the directory after the test function exits, by which point
    all references to the memmap have been released.
    """
    path = tmp_path / "test.tnt"
    path.write_bytes(_build_minimal_tnt(table_pairs, n_points, n_records))
    return str(path)


# =============================================================================
# Integration tests — real .tnt data
# =============================================================================

_LICL_TNT = os.path.join(DATA_DIR, 'tecmag', 'LiCl_ref1.tnt')
_LICL_TXT = os.path.join(DATA_DIR, 'tecmag', 'LiCl_ref1.txt')


@pytest.mark.slow
@pytest.mark.skipif(
    not os.path.exists(_LICL_TNT),
    reason="tecmag test data not available (download nmrglue test data)")
def test_tecmag_load_time_domain():
    """read() returns data matching the LiCl_ref1 reference text file."""
    ref1, data = ng.tecmag.read(_LICL_TNT)
    real, imag, usec = np.loadtxt(_LICL_TXT, skiprows=3, unpack=True)
    assert_array_almost_equal(data.real.squeeze(), real, decimal=3)
    assert_array_almost_equal(data.imag.squeeze(), imag, decimal=3)
    assert_array_almost_equal(
        np.arange(ref1['npts'][0]) * ref1['dwell'][0] * 1e6,
        usec, decimal=3)


# =============================================================================
# Unit tests — read()
# =============================================================================

def test_read_returns_dic_and_data(tmp_path):
    """read() returns a (dict, ndarray) pair with the correct shape and dtype."""
    path = _write_tnt(tmp_path, n_points=16, n_records=4)
    dic, data = ng.tecmag.read(path)
    assert isinstance(dic, dict)
    assert data.shape == (16, 4, 1, 1)
    assert data.dtype == np.dtype('<c8')
    assert dic['actual_npts'][0] == 16
    assert dic['actual_npts'][1] == 4
    assert_allclose(dic['sw'][0],      20000.0)
    assert_allclose(dic['ob_freq'][0], 21.31)


def test_read_rejects_non_tnt_magic(tmp_path):
    """read() raises ValueError on a file with an invalid magic header."""
    path = tmp_path / "bad.tnt"
    path.write_bytes(b'NOT_TNT_' + b'\x00' * 128)
    with pytest.raises(ValueError, match="Invalid magic number"):
        ng.tecmag.read(str(path))


def test_guess_udic_basic(tmp_path):
    """guess_udic() returns axis metadata matching the file's TMAG values."""
    path = _write_tnt(tmp_path, n_points=8, n_records=4)
    dic, data = ng.tecmag.read(path)
    udic = ng.tecmag.guess_udic(dic, data)
    assert udic[0]['size'] == 8
    assert udic[0]['sw']   == 20000.0
    assert_allclose(udic[0]['obs'], 21.31 * 1e6)


# =============================================================================
# Unit tests — guess_tables()
# =============================================================================

def test_guess_tables_basic_numeric(tmp_path):
    """A multi-value numeric table is found and parsed correctly."""
    path = _write_tnt(tmp_path, [('rfamp', '1.0 4.222 7.444 10.667 13.889 17.111')])
    tables = guess_tables(path)
    assert 'rfamp' in tables
    assert isinstance(tables['rfamp'], np.ndarray)
    assert len(tables['rfamp']) == 6
    assert_allclose(tables['rfamp'], [1.0, 4.222, 7.444, 10.667, 13.889, 17.111])


def test_guess_tables_si_microseconds(tmp_path):
    """SI-prefix 'u' tokens are converted to seconds (×1e-6)."""
    path = _write_tnt(tmp_path, [('tau', '5u 10u 20u 50u 100u')])
    assert_allclose(guess_tables(path)['tau'], [5e-6, 10e-6, 20e-6, 50e-6, 100e-6])


def test_guess_tables_si_nanoseconds(tmp_path):
    """SI-prefix 'n' tokens are converted to seconds (×1e-9)."""
    path = _write_tnt(tmp_path, [('tp', '100n 200n 500n 750n 999n')])
    assert_allclose(guess_tables(path)['tp'], [100e-9, 200e-9, 500e-9, 750e-9, 999e-9])


def test_guess_tables_si_milliseconds(tmp_path):
    """SI-prefix 'm' tokens are converted to seconds (×1e-3)."""
    path = _write_tnt(tmp_path, [('d1', '1m 2m 5m 10m 50m 100m')])
    assert_allclose(guess_tables(path)['d1'], [1e-3, 2e-3, 5e-3, 10e-3, 50e-3, 100e-3])


def test_guess_tables_si_seconds(tmp_path):
    """SI-prefix 's' tokens are converted to base seconds (×1.0)."""
    path = _write_tnt(tmp_path, [('last_delay', '1s 2s 5s 10s')])
    assert_allclose(guess_tables(path)['last_delay'], [1.0, 2.0, 5.0, 10.0])


def test_guess_tables_negative_values(tmp_path):
    """Negative float values (negative-only shape table) are parsed correctly."""
    path = _write_tnt(tmp_path, [('shape_neg', '-0.1 -0.5 -1.0 -0.5 -0.1')])
    got = guess_tables(path)['shape_neg']
    assert np.all(got <= 0)
    assert_allclose(got, [-0.1, -0.5, -1.0, -0.5, -0.1])


def test_guess_tables_bipolar_values(tmp_path):
    """Bipolar (positive + negative) shape tables are parsed correctly."""
    path = _write_tnt(tmp_path, [('shape_bip', '-0.5 -0.1 0.0 0.5 1.0 0.5 0.0 -0.1 -0.5')])
    got = guess_tables(path)['shape_bip']
    assert np.any(got > 0) and np.any(got < 0)
    assert_allclose(got, [-0.5, -0.1, 0.0, 0.5, 1.0, 0.5, 0.0, -0.1, -0.5])


def test_guess_tables_high_precision_floats(tmp_path):
    """High-precision floats (9 decimal places) round-trip within float64 tolerance."""
    values = [0.0, 0.123456789, 0.344343467, 0.577350269,
              0.707106781, 0.866025404, 0.951056516, 1.0]
    content = " ".join(f"{v:.9f}" for v in values)
    path = _write_tnt(tmp_path, [('rfamp_hp', content)])
    assert_allclose(guess_tables(path)['rfamp_hp'], values, atol=1e-7)


def test_guess_tables_signed_hz_offsets(tmp_path):
    """Plain signed float Hz values (no SI suffix) are parsed correctly."""
    values = [0.0, 1234.567890123, -987.654321098, 500.123456789, -333.444555666]
    content = " ".join(f"{v:.9f}" for v in values)
    path = _write_tnt(tmp_path, [('freq_offset', content)])
    assert_allclose(guess_tables(path)['freq_offset'], values, atol=1e-7)


def test_guess_tables_skips_scalars(tmp_path):
    """Single-token (scalar) entries are excluded from the result."""
    path = _write_tnt(tmp_path, [
        ('rfamp',    '1.0 5.0 10.0'),
        ('tp',       '50u'),
        ('deadtime', '20u'),
    ])
    tables = guess_tables(path)
    assert 'rfamp'    in tables
    assert 'tp'       not in tables
    assert 'deadtime' not in tables


def test_guess_tables_name_filter(tmp_path):
    """The `names` argument restricts results to the requested subset."""
    path = _write_tnt(tmp_path, [
        ('rfamp', '1.0 2.0 3.0'),
        ('d1',    '1m 2m 3m'),
        ('d2',    '1s 2s 3s'),
    ])
    tables = guess_tables(path, names=['rfamp', 'd1'])
    assert set(tables.keys()) == {'rfamp', 'd1'}


def test_guess_tables_convert_si_disabled(tmp_path):
    """When convert_si=False, values are returned as raw string token lists."""
    path = _write_tnt(tmp_path, [('d1', '1m 2m 3m')])
    assert guess_tables(path, convert_si=False)['d1'] == ['1m', '2m', '3m']


def test_guess_tables_no_pseq_section(tmp_path):
    """A file with no PSEQ section returns an empty dict without raising."""
    path = tmp_path / "no_pseq.tnt"
    path.write_bytes(b'TNT1.008' + b'\x00' * 256)
    assert guess_tables(str(path)) == {}


def test_guess_tables_rejects_invalid_length_prefix(tmp_path):
    """A name with a mismatched length prefix is silently skipped."""
    bogus = struct.pack('<i', 99) + b'rfamp' + b'extra junk'
    pseq  = _tlv_section('PSEQ', b'\x00' * 32 + bogus + b'\x00' * 32)
    path  = tmp_path / "bad_len.tnt"
    path.write_bytes(b'TNT1.008' + b'\x00' * 64 + pseq + b'\x00' * 64)
    assert 'rfamp' not in guess_tables(str(path))


def test_guess_tables_skips_non_ascii_names(tmp_path):
    """Names containing control bytes are rejected."""
    bad  = _pascal_string(b'\x01\x02bad') + _pascal_string('1.0 2.0 3.0')
    pseq = _tlv_section('PSEQ', b'\x00' * 32 + bad + b'\x00' * 32)
    path = tmp_path / "bad_name.tnt"
    path.write_bytes(b'TNT1.008' + b'\x00' * 64 + pseq + b'\x00' * 64)
    assert len(guess_tables(str(path))) == 0


def test_guess_tables_mixed_separators(tmp_path):
    """Content with mixed CRLF and space separators is parsed correctly."""
    path = _write_tnt(tmp_path, [('rfamp', '1.0\r\n2.0\r\n3.0 4.0 5.0')])
    assert_allclose(guess_tables(path)['rfamp'], [1.0, 2.0, 3.0, 4.0, 5.0])


def test_guess_tables_alongside_read(tmp_path):
    """read() and guess_tables() both work correctly on the same file."""
    path = _write_tnt(tmp_path, [('rfamp', '1.0 5.0 10.0 15.0')],
                      n_points=8, n_records=4)
    dic, data = ng.tecmag.read(path)
    tables    = guess_tables(path)
    assert data.shape == (8, 4, 1, 1)
    assert_allclose(tables['rfamp'], [1.0, 5.0, 10.0, 15.0])


# =============================================================================
# Unit tests — internal helpers
# =============================================================================

def test_convert_si_plain_floats():
    """_convert_si() passes plain floats through unchanged."""
    assert_allclose(_convert_si(['1.0', '2.5', '-3.14', '1e-6']),
                    [1.0, 2.5, -3.14, 1e-6])


def test_convert_si_with_prefixes():
    """_convert_si() applies SI prefix multipliers correctly."""
    assert_allclose(_convert_si(['100u', '5m', '21.31M', '2k']),
                    [100e-6, 5e-3, 21.31e6, 2e3])


def test_convert_si_all_time_prefixes():
    """_convert_si() handles n, u, m, s prefixes across many decades."""
    assert_allclose(_convert_si(['100n', '1u', '5m', '2s']),
                    [100e-9, 1e-6, 5e-3, 2.0])


def test_convert_si_invalid_token_raises():
    """_convert_si() raises ValueError on an unparseable token."""
    with pytest.raises(ValueError):
        _convert_si(['not_a_number'])


def test_read_pascal_pair_valid():
    """_read_pascal_pair() correctly parses a well-formed Pascal-string pair."""
    buf    = struct.pack('<i', 5) + b'rfamp' + struct.pack('<i', 3) + b'1 2'
    result = _read_pascal_pair(buf, 0)
    assert result is not None
    name, content, next_pos = result
    assert name    == 'rfamp'
    assert content == '1 2'
    assert next_pos == len(buf)


def test_read_pascal_pair_invalid_length():
    """_read_pascal_pair() returns None when the length prefix is out of range."""
    buf = struct.pack('<i', 99) + b'rfamp' + struct.pack('<i', 3) + b'1 2'
    assert _read_pascal_pair(buf, 0) is None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
