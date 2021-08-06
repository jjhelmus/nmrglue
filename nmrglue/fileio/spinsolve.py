"""
Functions for reading Magritek Spinsolve binary (dx/1d) files and 
parameter (acqu.par/proc.par) files.
"""

from warnings import warn
from pathlib import Path

import numpy as np

from . import fileiobase
from . import jcampdx

__developer_info__ = """
Spinsolve is the software used on the Magritek benchtop NMR devices. 

A spectrum is saved in a folder with several files. The spectral data is
stored in these files: 'data.1d' (FID), 'spectrum.1d' (Fourier transformed)
and 'spectrum_processed.1d' (FT + processed by spinsolve)
Optional spectral data (System->Prefs->Setup->Global data storage):
'nmr_fid.dx' (FID stored in `JCAMP-DX standard <http://www.jcamp-dx.org/>`), 
'spectrum.csv' and 'spectrum_processed.csv' (FT + processed by Spinsolve with ppm for each 
point and intensity delimited by ';')
Other files:
'acqu.par' - all parameters that are used for acquisition
'protocol.par' - text file used to reload data back into the Spinsolve software
'processing.script' - text file to transfer Spinsolve software protocol settings 
into MNOVA

The Spinsolve Expert software has a slightly different output:
[Needs to be double checked as I do not have access to this software -LCageman]
- Output into JCAMP-DX is not possible
- 'spectrum_processed.1d' is not generated
- (new) 'fid.1d' - seems to be the same as 'data.1d'
- (new) 'proc.par' - contains processing parameters in the same style as 'acqu.par'
- (new) .pt1 files - seem to be plot files specific for the expert software, cannot
be read by NMRglue
"""


def parse_spinsolve_acqu_line(line):
    """
    Parse lines in acqu.par and return a tuple (paramter name, parameter value)
    """
    line = line.strip()  # Drop newline
    name, value = line.split("=")  # Split at equal sign

    # remove spaces
    name = name.strip()
    value = value.strip()

    # Detect value type
    if value[0] == value[-1] == '"':  # String
        return name, str(value[1:-1])  # Drop quote marks
    if "." in value:  # Float
        return name, float(value)
    else:
        return name, int(value)


def read(dir='.', specfile=None, acqupar="acqu.par", procpar="proc.par"):
    """
    Reads spinsolve files from a directory
    When no spectrum filename is given (specfile), the following list is tried, in 
    that specific order
    ["nmr_fid.dx", "data.1d", "fid.1d", "spectrum.1d", "spectrum_processed.1d"]
    To use the resolution enhanced spectrum use the './Enhanced' folder as input.
    Note that spectrum.1d and spectrum_processed.1d contain only data in the 
    frequency domain, so no Fourier transformation is needed. Also, use 
    dic["spectrum"]["xaxis"] to plot the x-axis

    Parameters
    ----------
    dir : str, Path
        Directory to read from
    specfile : str, optional
        Filename to import spectral data from. None uses standard filename from: 
        ["nmr_fid.dx", "data.1d", "fid.1d", "spectrum.1d", "spectrum_processed.1d"]
    acqupar : str, optional
        Filename for acquisition parameters. None uses standard name.
    procpar : str, optional
        Filename for processing parameters. None uses standard name.

    Returns
    -------
    dic : dict
        All parameters that can be present in the data folder:
        dic["spectrum"] - First bytes of spectrum(_processed).1d
        dic["acqu"] - Parameters present in acqu.par
        dic["proc"] - Parameters present in proc.par
        dic["dx"] - - Parameters present in the header of nmr_fid.dx
    data : ndarray
        Array of NMR data

    """

    dir = Path(dir)
    if not dir.is_dir():
        raise IOError(f"Directory {dir} does not exist or is not a directory!")

    # Create empty dic
    dic = {
        "spectrum": {},
        "acqu": {},
        "proc": {},
        "dx": {}
    }

    # Read in acqu.par and write to dic
    acqupar = dir / acqupar
    if acqupar.is_file():
        with acqupar.open() as f:
            info = f.readlines()
        for line in info:
            par_name, par_value = parse_spinsolve_acqu_line(line)

            if par_name is not None:
                dic["acqu"][par_name] = par_value

    # Read in proc.par and write to dic
    procpar = dir / procpar
    if procpar.is_file():
        with acqupar.open() as f:
            info = f.readlines()
        for line in info:
            line = line.replace("\n", "")
            k, v = line.split("=")
            dic["proc"][k.strip()] = v.strip()

    # Define which spectrumfile to take, using 'specfile' when defined, otherwise 
    # the files in 'priority_list' are tried, in that particular order
    priority_list = ["nmr_fid.dx", "data.1d", "fid.1d", "spectrum.1d", "spectrum_processed.1d"]
    inputfile = None
    if specfile:
        inputfile = dir / specfile
        if not inputfile.is_file():
            raise FileNotFoundError(f"File {inputfile} does not exist")
    else:
        for priority in priority_list:
            inputfile = dir / priority
            if inputfile.is_file():
                break
    if inputfile is None:
        raise IOError(f"Directory {dir} does not contain spectral data")

    # Detect which file we are dealing with from the extension and read in the spectral data

    # Reading .dx file using existing nmrglue.fileio.jcampdx module
    if inputfile.suffix == ".dx":
        dic["dx"], raw_data = jcampdx.read(str(inputfile.absolute()))
        data = raw_data[0][:] + 1j * raw_data[1][:]

    # Reading .1d files
    elif inputfile.suffix == ".1d":
        with inputfile.open("rb") as f:
            raw_data = f.read()

        # Write out parameters from the first 32 bytes into dic["spectrum"]
        keys = ["owner", "format", "version", "dataType", "xDim", "yDim", "zDim", "qDim"]
        for i, k in enumerate(keys):
            start = i * 4
            end = start + 4
            value = int.from_bytes(raw_data[start:end], "little")
            dic["spectrum"][k] = value
        data = np.frombuffer(raw_data[32:], "<f")

        # The first 1/3 of the file is xaxis data (s or ppm)
        split = data.shape[-1] // 3
        xscale = data[0: split]
        dic["spectrum"]["xaxis"] = xscale

        # The rest is real and imaginary data points interleaved
        data = data[split:: 2] + 1j * data[split + 1:: 2]

    else:
        raise IOError(f"File {inputfile} cannot be interpreted, use .dx or .1d instead")

    return dic, data


def get_udic_from_acqu_dict(param: dict):
    """ Returns an udic from the parameters in acqu dictionary provided """
    return_dict = dict()

    # Spectral Width (sw, in Hz)
    return_dict['sw'] = float(param.get('bandwidth')) * 1000

    # Observation Frequency (obs, in MHz)
    return_dict['obs'] = float(param.get('b1Freq'))

    # Carrier frequency (car, in Hz)
    return_dict['car'] = float(param.get('lowestFrequency')) + (return_dict['sw'] / 2)

    # Label (str describing the axis name)
    return_dict['label'] = param['rxChannel']

    return return_dict


def get_udic_from_jcamp_dict(param):
    """ Returns an udic from the parameters in the JCAMP-DX dict provided """
    return_dict = dict()

    # Observation Frequency (obs, in MHz)
    try:
        return_dict['obs'] = float(param['$BF1'][0])
    except KeyError:
        warn("Cannot set observation frequency - set manually using: 'udic[0]['obs'] = x' "
             "where x is magnetic field in MHz")

    # Spectral Width (sw, in Hz)
    try:
        return_dict['sw'] = float(param['$SW'][0]) * return_dict['obs']
    except KeyError:
        warn("Cannot set spectral width - set manually using: 'udic[0]['sw'] = x' "
             "where x is the spectral width in Hz")

    # Carrier frequency (car, in Hz)
    try:
        return_dict['car'] = -float(param['$REFERENCEPOINT'][0]) + (return_dict['sw'] / 2)
    except KeyError:
        pass
    try:
        return_dict['car'] = (return_dict['obs'] - float(param['$SF'][0])) * 1e6
    except KeyError:
        warn("Cannot set carrier - try: 'udic[0]['car'] = x * udic[0]['obs']' "
             "where x is the center of the spectrum in ppm")

    # Label (str describing the axis name)
    try:
        return_dict['label'] = param[".OBSERVENUCLEUS"][0].replace("^", "")
    except KeyError:
        warn("Cannot set observed nucleus label")

    return return_dict


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of JCAMP-DX, acqu, proc and spectrum parameters.
    data : ndarray
        Array of NMR data.
    
    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.
    """

    # Create an empty universal dictionary
    udic = fileiobase.create_blank_udic(1)

    # Size can be calculated from data len
    if data is not None:
        udic[0]["size"] = len(data)
    else:
        warn('No data, cannot set udic size')

    # For the other parameters, first try to use info from acqu.par (dic['acqu']), then JCAMP-DX header parameters
    if dic["acqu"]:
        udic[0].update(get_udic_from_acqu_dict(dic["acqu"]))
    else:
        udic[0].update(get_udic_from_jcamp_dict(dic["dx"]))

    # keys left to default
    # udic[0]['complex']
    # udic[0]['encoding']
    # udic[0]['time'] = True
    # udic[0]['freq'] = False
    return udic
