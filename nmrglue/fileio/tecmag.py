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
import re

import numpy as np

from . import fileiobase


TNTMAGIC_RE = re.compile(br"^TNT1\.\d\d\d$")

TNTMAGIC = np.dtype('a8')
TNTTLV = np.dtype([('tag', 'a4'), ('bool', '<u4'), ('length', '<u4')])

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
    # ('space1', 'a0'),

    ('magnet_field', '<f8'),
    ('ob_freq', '<f8', 4),
    ('base_freq', '<f8', 4),
    ('offset_freq', '<f8', 4),
    ('ref_freq', '<f8'),
    ('NMR_frequency', '<f8'),
    ('obs_channel', '<i2'),
    ('space2', 'a42'),

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
    ('space3', 'a16'),

    ('transmitter_gain', '<i2'),
    ('receiver_gain', '<i2'),
    ('NumberOfReceivers', '<i2'),
    ('RG2', '<i2'),
    ('receiver_phase', '<f8'),
    ('space4', 'a4'),

    ('set_spin_rate', '<u2'),
    ('actual_spin_rate', '<u2'),

    ('lock_field', '<i2'),
    ('lock_power', '<i2'),
    ('lock_gain', '<i2'),
    ('lock_phase', '<i2'),
    ('lock_freq_mhz', '<f8'),
    ('lock_ppm', '<f8'),
    ('H2O_freq_ref', '<f8'),
    ('space5', 'a16'),

    ('set_temperature', '<f8'),
    ('actual_temperature', '<f8'),

    ('shim_units', '<f8'),
    ('shims', '<i2', 36),
    ('shim_FWHM', '<f8'),

    ('HH_dcpl_attn', '<i2'),
    ('DF_DN', '<i2'),
    ('F1_tran_mode', '<i2', 7),
    ('dec_BW', '<i2'),
    ('grd_orientation', 'a4'),
    ('LatchLP', '<i4'),
    ('grd_Theta', '<f8'),
    ('grd_Phi', '<f8'),
    ('space6', 'a264'),

    ('start_time', '<u4'),
    ('finish_time', '<u4'),
    ('elapsed_time', '<i4'),

    ('date', 'a32'),
    ('nuclei', 'a16', 4),
    ('sequence', 'a32'),
    ('lock_solvent', 'a16'),
    ('lock_nucleus', 'a16')
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
    ('axisName', 'a32'),
    ('space', 'a52'),
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

    ('space', 'a610'),

    ('line_simulation_name', 'a32'),
    ('integral_template_name', 'a32'),
    ('baseline_template_name', 'a32'),
    ('layout_name', 'a32'),
    ('relax_information_name', 'a32'),
    ('username', 'a32'),
    ('user_string_1', 'a16'),
    ('user_string_2', 'a16'),
    ('user_string_3', 'a16'),
    ('user_string_4', 'a16')
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

        tntmagic = np.fromstring(tntfile.read(TNTMAGIC.itemsize),
                                 TNTMAGIC, count=1)[0]

        if not TNTMAGIC_RE.match(tntmagic):
            err = ("Invalid magic number (is '%s' really TNMR file?): %s" %
                   (filename, tntmagic))
            raise ValueError(err)

        # Read in the section headers
        tnthdrbytes = tntfile.read(TNTTLV.itemsize)
        while(TNTTLV.itemsize == len(tnthdrbytes)):
            tlv = np.fromstring(tnthdrbytes, TNTTLV)[0]
            data_length = tlv['length']
            hdrdict = {'offset': tntfile.tell(),
                       'length': data_length,
                       'bool': bool(tlv['bool'])}
            if data_length <= 4096:
                hdrdict['data'] = tntfile.read(data_length)
                assert(len(hdrdict['data']) == data_length)
            else:
                tntfile.seek(data_length, os.SEEK_CUR)
            tnt_sections[tlv['tag'].decode()] = hdrdict
            tnthdrbytes = tntfile.read(TNTTLV.itemsize)

    assert(tnt_sections['TMAG']['length'] == TNTTMAG.itemsize)
    tmag = np.fromstring(tnt_sections['TMAG']['data'], TNTTMAG, count=1)[0]

    assert(tnt_sections['DATA']['length'] ==
           tmag['actual_npts'].prod() * 8)
    #  For some reason we can't set offset and shape together
    # DATA = np.memmap(tntfilename,np.dtype('<c8'), mode='r',
    #                  offset=self.tnt_sections['DATA']['offset'],
    #                  shape=self.TMAG['actual_npts'].tolist(),order='F')
    data = np.memmap(filename, np.dtype('<c8'), mode='c',
                     offset=tnt_sections['DATA']['offset'],
                     shape=tmag['actual_npts'].prod())
    data = np.reshape(data, tmag['actual_npts'], order='F')

    assert(tnt_sections['TMG2']['length'] == TNTTMG2.itemsize)
    tmg2 = np.fromstring(tnt_sections['TMG2']['data'], TNTTMG2, count=1)[0]

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
