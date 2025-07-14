"""
Functions for reading `nmrML <http://nmrml.org/>`_ files.
"""

import base64
import zlib
from warnings import warn

import numpy as np


def read(filename, data_dtype=None, read_processed=False):
    """
    Read a nmrML file.

    Parameters
    ----------
    filename : str
        Name of nmrML file to read.
    data_dtype : str, optional
        NumPy data type of the data. None, the default, will determine this
        data type from the information in the file. Occasionally this
        information is incorrect and this argument can be used to explicitly
        supply this information.
    read_processed : bool
        Option to read processed data instead of FID (dafault=False)

    Returns
    -------
    dic : dict
        Dictionary of spectra parameters.
    data : ndarray, 1D
        Array of spectral data.

    """
    import xmltodict  # delay import so that xmltodict is optional
    with open(filename, 'rb') as nmrml_file:
        doc = xmltodict.parse(nmrml_file.read())

    if read_processed:
        data_dict = doc['nmrML']['spectrumList']['spectrum1D']['spectrumDataArray']
    else:
        data_dict = doc['nmrML']['acquisition']['acquisition1D']['fidData']
    data = _get_nmrml_data(data_dict, data_dtype)

    dic = doc['nmrML']
    return dic, data


BYTE_FORMAT_TO_DTYPE = {
    'float32': '<f4',
    'float64': '<f8',
    'complex64': '<c8',
    'complex128': '<c16',
}


def _get_nmrml_data(fid_dict, data_dtype):
    """ Return a NumPy array for spectral data in a nmrML data file. """
    # read required parameters from dictionary
    byte_format = fid_dict['@byteFormat'].lower()
    compressed = fid_dict['@compressed'].lower() == 'true'
    encoded_length = int(fid_dict['@encodedLength'])
    base64_text = fid_dict['#text']
    if len(base64_text) != encoded_length:
        warn("Size of encoded text (%i) does not match stated length (%i)"
             % (len(base64_text), encoded_length))

    # decode and uncompress
    compressed_bytes = base64.b64decode(base64_text)
    if compressed:
        uncompressed_bytes = zlib.decompress(compressed_bytes)
    else:
        uncompressed_bytes = compressed_bytes

    # convert to ndarray
    if data_dtype is None:
        data_dtype = BYTE_FORMAT_TO_DTYPE[byte_format]
    data = np.frombuffer(uncompressed_bytes, dtype=data_dtype)

    return data
