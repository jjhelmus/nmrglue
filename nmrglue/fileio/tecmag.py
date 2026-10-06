"""
Functions for reading Tecmag .tnt data files.
"""
__developer_doc__ = """
Tecmag .tnt file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Tecmag .tnt file format is documented with C pseudo-code in the file
"A1 - TNMR File Format.doc" distributed with the TNMR software.

This file is based on the
`pytnt module <https://www.github.com/chatcannon/pytnt>`_.
Please inform upstream if you find a bug or to request additional features.

"""

import io
import os
import struct
import re

import numpy as np

from . import fileiobase


TNTMAGIC_RE = re.compile(br"^TNT1\.\d\d\d$")

TNTMAGIC = np.dtype('S8')
TNTTLV = np.dtype([('tag', 'S4'), ('bool', '<u4'), ('length', '<u4')])

TNTTMAG = np.dtype([
    ('npts', '<i4', 4),
    ('actual_npts', '<i4', 4),
    ('acq_points', '<i4'),
    ('npts_start', '<i4', 4),
    ('scans', '<i4'),
    ('actual_scans', '<i4'),
    ('dummy_scans', '<i4'),
    ('repeat_times', '<i4'),
    ('sadimension', '<i4'),
    ('samode', '<i4'),
    # ('space1', 'S0'),

    ('magnet_field', '<f8'),
    ('ob_freq', '<f8', 4),
    ('base_freq', '<f8', 4),
    ('offset_freq', '<f8', 4),
    ('ref_freq', '<f8'),
    ('NMR_frequency', '<f8'),
    ('obs_channel', '<i2'),
    ('space2', 'S42'),

    ('sw', '<f8', 4),
    ('dwell', '<f8', 4),
    ('filter', '<f8'),
    ('experiment_time', '<f8'),
    ('acq_time', '<f8'),
    ('last_delay', '<f8'),
    ('spectrum_direction', '<i2'),
    ('hardware_sideband', '<i2'),
    ('Taps', '<i2'),
    ('Type', '<i2'),
    ('bDigRec', '<u4'),
    ('nDigitalCenter', '<i4'),
    ('space3', 'S16'),

    ('transmitter_gain', '<i2'),
    ('receiver_gain', '<i2'),
    ('NumberOfReceivers', '<i2'),
    ('RG2', '<i2'),
    ('receiver_phase', '<f8'),
    ('space4', 'S4'),

    ('set_spin_rate', '<u2'),
    ('actual_spin_rate', '<u2'),

    ('lock_field', '<i2'),
    ('lock_power', '<i2'),
    ('lock_gain', '<i2'),
    ('lock_phase', '<i2'),
    ('lock_freq_mhz', '<f8'),
    ('lock_ppm', '<f8'),
    ('H2O_freq_ref', '<f8'),
    ('space5', 'S16'),

    ('set_temperature', '<f8'),
    ('actual_temperature', '<f8'),

    ('shim_units', '<f8'),
    ('shims', '<i2', 36),
    ('shim_FWHM', '<f8'),

    ('HH_dcpl_attn', '<i2'),
    ('DF_DN', '<i2'),
    ('F1_tran_mode', '<i2', 7),
    ('dec_BW', '<i2'),
    ('grd_orientation', 'S4'),
    ('LatchLP', '<i4'),
    ('grd_Theta', '<f8'),
    ('grd_Phi', '<f8'),
    ('space6', 'S264'),

    ('start_time', '<u4'),
    ('finish_time', '<u4'),
    ('elapsed_time', '<i4'),

    ('date', 'S32'),
    ('nuclei', 'S16', 4),
    ('sequence', 'S32'),
    ('lock_solvent', 'S16'),
    ('lock_nucleus', 'S16')
])


TNTGRIDANDAXIS = np.dtype([
    ('majorTickInc', '<f8', 12),
    ('minorIntNum', '<i2', 12),
    ('labelPrecision', '<i2', 12),
    ('gaussPerCentimeter', '<f8'),
    ('gridLines', '<i2'),
    ('axisUnits', '<i2'),
    ('showGrid', '<u4'),
    ('showGridLabels', '<u4'),
    ('adjustOnZoom', '<u4'),
    ('showDistanceUnits', '<u4'),
    ('axisName', 'S32'),
    ('space', 'S52'),
])


TNTTMG2 = np.dtype([
    ('real_flag', '<u4'),
    ('imag_flag', '<u4'),
    ('magn_flag', '<u4'),
    ('axis_visible', '<u4'),
    ('auto_scale', '<u4'),
    ('line_display', '<u4'),
    ('show_shim_units', '<u4'),

    ('integral_display', '<u4'),
    ('fit_display', '<u4'),
    ('show_pivot', '<u4'),
    ('label_peaks', '<u4'),
    ('keep_manual_peaks', '<u4'),
    ('label_peaks_in_units', '<u4'),
    ('integral_dc_average', '<u4'),
    ('integral_show_multiplier', '<u4'),
    ('Boolean_space', '<u4', 9),

    ('all_ffts_done', '<u4', 4),
    ('all_phase_done', '<u4', 4),

    ('amp', '<f8'),
    ('ampbits', '<f8'),
    ('ampCtl', '<f8'),
    ('offset', '<i4'),

    ('axis_set', TNTGRIDANDAXIS),

    ('display_units', '<i2', 4),
    ('ref_point', '<i4', 4),
    ('ref_value', '<f8', 4),
    ('z_start', '<i4'),
    ('z_end', '<i4'),
    ('z_select_start', '<i4'),
    ('z_select_end', '<i4'),
    ('last_zoom_start', '<i4'),
    ('last_zoom_end', '<i4'),
    ('index_2D', '<i4'),
    ('index_3D', '<i4'),
    ('index_4D', '<i4'),

    ('apodization_done', '<i4', 4),
    ('linebrd', '<f8', 4),
    ('gaussbrd', '<f8', 4),
    ('dmbrd', '<f8', 4),
    ('sine_bell_shift', '<f8', 4),
    ('sine_bell_width', '<f8', 4),
    ('sine_bell_skew', '<f8', 4),
    ('Trapz_point_1', '<i4', 4),
    ('Trapz_point_2', '<i4', 4),
    ('Trapz_point_3', '<i4', 4),
    ('Trapz_point_4', '<i4', 4),
    ('trafbrd', '<f8', 4),
    ('echo_center', '<i4', 4),

    ('data_shift_points', '<i4'),
    ('fft_flag', '<i2', 4),
    ('unused', '<f8', 8),
    ('pivot_point', '<i4', 4),
    ('cumm_0_phase', '<f8', 4),
    ('cumm_1_phase', '<f8', 4),
    ('manual_0_phase', '<f8'),
    ('manual_1_phase', '<f8'),
    ('phase_0_value', '<f8'),
    ('phase_1_value', '<f8'),
    ('session_phase_0', '<f8'),
    ('session_phase_1', '<f8'),

    ('max_index', '<i4'),
    ('min_index', '<i4'),
    ('peak_threshold', '<f4'),
    ('peak_noise', '<f4'),
    ('integral_dc_points', '<i2'),
    ('integral_label_type', '<i2'),
    ('integral_scale_factor', '<f4'),
    ('auto_integrate_shoulder', '<i4'),
    ('auto_integrate_noise', '<f8'),
    ('auto_integrate_threshold', '<f8'),
    ('s_n_peak', '<i4'),
    ('s_n_noise_start', '<i4'),
    ('s_n_noise_end', '<i4'),
    ('s_n_calculated', '<f4'),

    ('Spline_point', '<i4', 14),
    ('Spline_point_avr', '<i2'),
    ('Poly_point', '<i4', 8),
    ('Poly_point_avr', '<i2'),
    ('Poly_order', '<i2'),

    ('space', 'S610'),

    ('line_simulation_name', 'S32'),
    ('integral_template_name', 'S32'),
    ('baseline_template_name', 'S32'),
    ('layout_name', 'S32'),
    ('relax_information_name', 'S32'),
    ('username', 'S32'),
    ('user_string_1', 'S16'),
    ('user_string_2', 'S16'),
    ('user_string_3', 'S16'),
    ('user_string_4', 'S16')
])


def read(filename):
    """
    Read a Tecmag .tnt data file.

    Parameters
    ----------
    filename : str
        Name of file to read from

    Returns
    -------
    dic : dict
        Dictionary of Tecmag parameters.
    data : ndarray
        Array of NMR data.

    """
    tnt_sections = dict()

    with open(filename, 'rb') as tntfile:

        tntmagic = np.frombuffer(tntfile.read(TNTMAGIC.itemsize),
                                 TNTMAGIC, count=1)[0]

        if not TNTMAGIC_RE.match(tntmagic):
            err = ("Invalid magic number (is '%s' really TNMR file?): %s" %
                   (filename, tntmagic))
            raise ValueError(err)

        # Read in the section headers
        tnthdrbytes = tntfile.read(TNTTLV.itemsize)
        while TNTTLV.itemsize == len(tnthdrbytes):
            tlv = np.frombuffer(tnthdrbytes, TNTTLV)[0]
            data_length = tlv['length']
            hdrdict = {'offset': tntfile.tell(),
                       'length': data_length,
                       'bool': bool(tlv['bool'])}
            if data_length <= 4096:
                hdrdict['data'] = tntfile.read(data_length)
                assert len(hdrdict['data']) == data_length
            else:
                tntfile.seek(data_length, os.SEEK_CUR)
            tnt_sections[tlv['tag'].decode()] = hdrdict
            tnthdrbytes = tntfile.read(TNTTLV.itemsize)

    assert tnt_sections['TMAG']['length'] == TNTTMAG.itemsize
    tmag = np.frombuffer(tnt_sections['TMAG']['data'], TNTTMAG, count=1)[0]

    assert (tnt_sections['DATA']['length'] ==
            tmag['actual_npts'].prod() * 8)
    #  For some reason we can't set offset and shape together
    # DATA = np.memmap(tntfilename,np.dtype('<c8'), mode='r',
    #                  offset=self.tnt_sections['DATA']['offset'],
    #                  shape=self.TMAG['actual_npts'].tolist(),order='F')
    data = np.memmap(filename, np.dtype('<c8'), mode='c',
                     offset=tnt_sections['DATA']['offset'],
                     shape=tmag['actual_npts'].prod())
    data = np.reshape(data, tmag['actual_npts'], order='F')

    assert tnt_sections['TMG2']['length'] == TNTTMG2.itemsize
    tmg2 = np.frombuffer(tnt_sections['TMG2']['data'], TNTTMG2, count=1)[0]

    dic = dict()
    for name in TNTTMAG.names:
        if not name.startswith('space'):
            dic[name] = tmag[name]
    for name in TNTTMG2.names:
        if name not in ['Boolean_space', 'unused', 'space', 'axis_set']:
            dic[name] = tmg2[name]
    for name in TNTGRIDANDAXIS.names:
        dic[name] = tmg2['axis_set'][name]

    return dic, data


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of Tecmag parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.

    """

    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(4)

    # update default values
    for i in range(4):
        udic[i]["size"] = dic['actual_npts'][i]
        udic[i]["sw"] = dic['sw'][i]
        udic[i]["complex"] = True
        udic[i]["obs"] = dic['ob_freq'][i] * 1e6
        # Not sure what the difference is here
        # N.B. base_freq is some bogus value like 1.4e-13
        udic[i]["car"] = dic['ob_freq'][i] * 1e6
        udic[i]["time"] = not bool(dic['fft_flag'][i])
        udic[i]["freq"] = bool(dic['fft_flag'][i])

    return udic


# -----------------------------------------------------------------------------
# Sequence-table reading from PSEQ section (added for issue #123)
# -----------------------------------------------------------------------------

# SI prefix table used to convert tokens like '100u' -> 1e-4.
# Same set of prefixes used by pytnt (chatcannon/pytnt, GPL-3.0).
_SI_PREFIX = {
    'y': 1e-24, 'z': 1e-21, 'a': 1e-18, 'f': 1e-15, 'p': 1e-12,
    'n': 1e-9,  'u': 1e-6,  'm': 1e-3,  'c': 1e-2,  'd': 1e-1,
    's': 1.0,   'k': 1e3,   'M': 1e6,   'G': 1e9,   'T': 1e12,
    'P': 1e15,  'E': 1e18,  'Z': 1e21,  'Y': 1e24,
}

# Maximum allowed length of a single table name (characters)
_NAME_MAX_LEN = 64
# Maximum allowed length of a single table's content (bytes)
_CONTENT_MAX_LEN = 1_000_000


def _convert_si(tokens):
    """Convert a list of strings to a numpy float64 array, handling SI suffixes.

    Each token is parsed as a plain float; if that fails and the last character
    is an SI prefix (y, z, ..., Y), the rest is parsed as float and multiplied
    by the prefix value. Raises ValueError on unparseable tokens.
    """
    out = np.empty(len(tokens), dtype=np.float64)
    for i, t in enumerate(tokens):
        try:
            out[i] = float(t)
        except ValueError:
            if t and t[-1] in _SI_PREFIX:
                out[i] = _SI_PREFIX[t[-1]] * float(t[:-1])
            else:
                raise ValueError(
                    "could not convert {0!r} to a number".format(t))
    return out


def _read_pascal_pair(buf, pos):
    """Try to read two consecutive Pascal strings at byte offset ``pos``.

    A Pascal string is a 4-byte little-endian length followed by that many
    bytes of ASCII content. Returns ``(name, content, next_pos)`` if both
    strings validate, else ``None``.

    Validation rules:
      - Name length must be 1..64 inclusive.
      - Name bytes must be printable ASCII (0x20..0x7e).
      - Name must not begin with a digit (table names start with a letter
        or underscore in TNMR).
      - Content length must be 1..1_000_000 inclusive.
      - Content bytes must be printable ASCII or whitespace.
    """
    if pos < 0 or pos + 4 > len(buf):
        return None
    name_len = struct.unpack_from('<i', buf, pos)[0]
    if not 1 <= name_len <= _NAME_MAX_LEN:
        return None
    name_bytes = buf[pos + 4:pos + 4 + name_len]
    if len(name_bytes) != name_len:
        return None
    if not all(0x20 <= b <= 0x7e for b in name_bytes):
        return None
    if name_bytes[0:1].isdigit():
        return None
    try:
        name = name_bytes.decode('ascii')
    except UnicodeDecodeError:
        return None

    content_pos = pos + 4 + name_len
    if content_pos + 4 > len(buf):
        return None
    content_len = struct.unpack_from('<i', buf, content_pos)[0]
    if not 1 <= content_len <= _CONTENT_MAX_LEN:
        return None
    content_bytes = buf[content_pos + 4:content_pos + 4 + content_len]
    if len(content_bytes) != content_len:
        return None
    if not all(b == 0 or 0x09 <= b <= 0x7e for b in content_bytes):
        return None
    content = content_bytes.decode('latin1').rstrip('\x00')

    return name, content, content_pos + 4 + content_len


def guess_tables(filename, names=None, convert_si=True):
    """Discover sequence tables (rfamp, tp, delay tables, ...) in a .tnt file.

    Tables are stored in the PSEQ section as pairs of length-prefixed Pascal
    strings: ``[name_length: int32][name: bytes][content_length: int32]
    [content: bytes]``. The content is either whitespace- or CRLF-separated
    numeric tokens, optionally with SI-prefix suffixes (e.g. ``'100u'`` for
    100 microseconds).

    Parameters
    ----------
    filename : str
        Name of the .tnt file to read.
    names : list of str, optional
        If given, restrict results to tables with these names. If None,
        every detectable table is returned.
    convert_si : bool, optional
        If True (default), convert tokens with SI prefixes (e.g. ``'5m'``) to
        base units (``5e-3``). If False, table values are returned as the
        raw whitespace-split tokens (list of str).

    Returns
    -------
    tables : dict
        Mapping from table name (str) to either:

        - ``numpy.ndarray`` of float values when ``convert_si=True``, or
        - ``list`` of str tokens when ``convert_si=False``.

        Only multi-valued tables are returned; single-value (scalar) entries
        are excluded since they typically represent sequence variable
        defaults rather than swept tables.

    Notes
    -----
    The Pascal-string layout was first identified by the pytnt project
    (chatcannon/pytnt, GPL-3.0). pytnt's reader only matches tables whose
    name follows the default TNMR delay-table pattern ``de[0-9]+:[0-9]``;
    this implementation extends to arbitrary user-given names and adds
    strict length-prefix validation to reject coincidental substring
    matches.

    The function locates the PSEQ section by finding the ``b'PSEQ'`` tag
    and skipping the 8-byte TLV header (bool + length). Inside the section,
    it byte-walks looking for valid Pascal-string pairs.

    Examples
    --------
    >>> import nmrglue as ng
    >>> tables = ng.tecmag.guess_tables('data.tnt')
    >>> sorted(tables.keys())
    ['d1', 'rfamp']
    >>> tables['rfamp']
    array([ 1.  ,  4.22,  7.44, 10.67, 13.89, 17.11, 20.33, 23.56, 26.78, 30.  ])
    >>> ng.tecmag.guess_tables('data.tnt', names=['rfamp'])
    {'rfamp': array([1.  , 4.22, 7.44, ...])}
    """
    with open(filename, 'rb') as f:
        data = f.read()

    # Locate PSEQ section payload (skip 4-byte tag + 4-byte bool + 4-byte length)
    pseq_pos = data.find(b'PSEQ')
    if pseq_pos < 0:
        return {}
    region_start = pseq_pos + 4 + 4 + 4

    name_filter = set(names) if names is not None else None
    tables = {}

    pos = region_start
    while pos < len(data) - 8:
        pair = _read_pascal_pair(data, pos)
        if pair is None:
            pos += 1
            continue
        name, content, next_pos = pair

        # Skip if user requested a name filter and this isn't on it
        if name_filter is not None and name not in name_filter:
            pos = next_pos
            continue
        # Skip duplicate names; keep first occurrence
        if name in tables:
            pos = next_pos
            continue
        # Skip scalars: whitespace-split content must have >1 token
        tokens = content.split()
        if len(tokens) <= 1:
            pos = next_pos
            continue

        if convert_si:
            try:
                tables[name] = _convert_si(tokens)
            except ValueError:
                # Non-numeric content; keep as token list
                tables[name] = tokens
        else:
            tables[name] = tokens
        pos = next_pos

    return tables
