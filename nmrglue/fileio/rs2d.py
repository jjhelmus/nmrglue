"""
Functions for reading RS2D data format.
"""

import os
from warnings import warn

import numpy as np

from . import fileiobase

__developer_info__ = """
RS2D file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
NMR spectra in RS2D format are stored in two files:
- "data.dat"    Binary data in complex64 or float32 (big-endian)
- "header.xml"  Acquisition/processing parameters in XML format

The complete format specification with test files is distributed
in nmrglue test data package.
"""


def read(path='.', datafile="data.dat", headerfile="header.xml"):
    """
    Reads RS2D files from a directory.


    Parameters
    ----------
    path : str
        Directory to read from
    datafile : str, optional
        Filename to import spectral data from. Default name is used if not given.
    headerfile : str, optional
        Filename for acquisition/processing parameters. Default name is used if not given.


    Returns
    -------
    dic : dict
        All parameters read from the header xml file
    data : ndarray
        Array of NMR data, either FID or processed. The data will be in shape of
        [RECEIVER COUNT, DIM4, DIM3, DIM2, DIM1], but if any of these dimensions
        is 1, it is left out from shape.

    """

    datafilepath = os.path.join(path, datafile)
    if os.path.isfile(datafilepath) is not True:
        raise OSError(f"data file {datafilepath} does not exist")

    headerfilepath = os.path.join(path, headerfile)
    if os.path.isfile(headerfilepath) is not True:
        raise OSError(f"header file {headerfilepath} does not exist")

    import xmltodict  # delay import so that xmltodict is optional
    with open(headerfilepath, 'rb') as nmrml_file:
        rawdict = xmltodict.parse(nmrml_file.read())

    # clean parameter dictionary

    # only read the following relevant parameters from each entry:
    params_to_read = ("name", "description", "group",
                      "category", "numberEnum", "value",
                      "uuid", "uuidBaseFrequency")

    dic = {}
    for entry in rawdict["header"]["params"]["entry"]:
        key = entry["key"]
        dic[key] = {}
        for param in entry["value"]:
            if param not in params_to_read:
                continue
            if param in dic[key]:
                # encountering list (multiple same param names)
                if not isinstance(dic[key][param], list):
                    # turn this param to list
                    dic[key][param] = [dic[key][param]]
                dic[key][param].append(entry["value"][param])
            # first or singular value
            dic[key][param] = entry["value"][param]

    # find dimensions
    try:
        receiver_count = int(dic['RECEIVER_COUNT']["value"])
        dim1d = int(dic['MATRIX_DIMENSION_1D']["value"])
        dim2d = int(dic['MATRIX_DIMENSION_2D']["value"])
        dim3d = int(dic['MATRIX_DIMENSION_3D']["value"])
        dim4d = int(dic['MATRIX_DIMENSION_4D']["value"])
    except (KeyError, ValueError) as e:
        raise ValueError("Cannot find data dimensions from RS2D header") from e

    # make shape
    shape = [receiver_count, dim4d, dim3d, dim2d, dim1d]
    shape = [i for i in shape if i > 1]

    # check data type
    is_real = False  # default = complex64
    try:
        datareps = dic["DATA_REPRESENTATION"]["value"]
        if isinstance(datareps, list):
            is_real = datareps[0] == "REAL"
    except KeyError:
        pass

    # read and reshape data
    if is_real:
        data = np.fromfile(datafilepath, dtype=">f4")
    else:
        data = np.fromfile(datafilepath, dtype=">c8")
    data = np.reshape(data, shape)

    return dic, data


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Clean dictionary of RS2D parameters from rs2d.read().
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.
    """

    uuid_dict = {}
    for key, entry in dic.items():
        uuid = entry.get("uuid", None)
        if uuid:
            uuid_dict[uuid] = entry

    # create an empty universal dictionary
    dims = len(data.shape)
    # shape may have extra dimension from possible multiple receivers.
    # in that case decrease one dimension from udic.
    receiver_count = int(dic['RECEIVER_COUNT']["value"])
    if receiver_count > 1:
        dims -= 1
    udic = fileiobase.create_blank_udic(dims)

    # set udic parameters:

    # "size"
    # note: shape is already in correct order (direct dimension last)
    i = 0
    for i in range(dims):
        # skip receiver count dimension if present
        dim = data.shape[i+1] if receiver_count > 1 else data.shape[i]
        udic[i]["size"] = dim

    # "label"
    for i, key in enumerate(["NUCLEUS_1", "NUCLEUS_2", "NUCLEUS_3", "NUCLEUS_4"]):

        label = None
        try:
            label = dic[key]["value"]
            udic[dims-i-1]["label"] = label
        except KeyError:
            pass

    # "obs"
    for i, key in enumerate(["BASE_FREQ_1", "BASE_FREQ_2", "BASE_FREQ_3", "BASE_FREQ_4"]):
        obs_freq = None
        try:
            obs_freq = float(dic[key]["value"])
            udic[dims-i-1]["obs"] = obs_freq/1e6
        except ValueError:
            warn(f"Cannot parse {key}")
        except KeyError:
            pass

    # "sw"
    for i, key in enumerate(["SPECTRAL_WIDTH", "SPECTRAL_WIDTH_2D",
                             "SPECTRAL_WIDTH_3D", "SPECTRAL_WIDTH_4D"]):
        sw = None
        try:
            sw = float(dic[key]["value"])
            unitname = dic[key].get("numberEnum", None)
            # sw needs conversion to Hz if stored in ppm:
            if unitname == "FrequencyShift":
                # conversion frequency may be user defined. find
                # correct one by given uuid:
                try:
                    uuid_basefreq = dic[key]["uuidBaseFrequency"]
                    freq_entry = uuid_dict[uuid_basefreq]
                    basefreq = float(freq_entry["value"])
                    sw = sw * (basefreq / 1e6)
                except (KeyError, ValueError):
                    warn("Cannot set sw: ppm to Hz conversion failed")
                    sw = None
            if sw:
                udic[dims-i-1]["sw"] = sw
        except ValueError:
            warn(f"Cannot parse {key}")
        except KeyError:
            pass

    # "car"
    # most conveniently found by calculating
    # TRANSMIT_FREQ - BASE_FREQ, which are already in Hz
    for i in range(dims):
        car = None
        tag_basefreq = "BASE_FREQ_"+str(i+1)
        tag_transmitfreq = "TRANSMIT_FREQ_"+str(i+1)
        try:
            basefreqhz = float(dic[tag_basefreq]["value"])
            transmitfreqhz = float(dic[tag_transmitfreq]["value"])
            car = transmitfreqhz - basefreqhz
            udic[dims-i-1]["car"] = car
        except ValueError:
            warn(f"Cannot parse {tag_basefreq} or {tag_transmitfreq}")
        except KeyError:
            pass

    # "time" and "freq"

    # default is FID i.e. time domain:
    for i in range(dims):
        udic[i]["freq"] = False
        udic[i]["time"] = True

    # if STATE is populated, dimensions with value 1 are in frequency domain:
    try:
        states = dic["STATE"]["value"]
        if isinstance(states, list):
            for i, state in enumerate(states):
                is_processed = int(state) == 1
                udic[dims-i-1]["freq"] = is_processed
                udic[dims-i-1]["time"] = not is_processed
    except KeyError:
        pass

    # "complex"

    # default representation is complex values.
    # if DATA_REPRESENTATION is populated, dimensions with value "REAL"
    # are represented in real values:
    try:
        datareps = dic["DATA_REPRESENTATION"]["value"]
        if isinstance(datareps, list):
            for i, datarep in enumerate(datareps):
                is_real = datarep == "REAL"
                udic[dims-i-1]["complex"] = not is_real
    except KeyError:
        pass

    return udic
