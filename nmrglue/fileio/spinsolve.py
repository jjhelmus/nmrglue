"""
Functions for reading Magritek Spinsolve binary (dx/1d) files and 
parameter (acqu.par/proc.par) files.
"""

import os
from warnings import warn

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
'spectrum.csv' and 'spectrum_processed.csv' (FT + processed by Spinsovle with ppm for each 
point and intensity delimited by ';')
Other files:
'acqu.par' - all parameters that are used for acquisition
'Protocol.par' - text file used to reload data back into the Spinsolve software
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
    dir : str
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
    
    if os.path.isdir(dir) is not True:
        raise IOError("directory %s does not exist" % (dir))

    # Create empty dic
    dic = {"spectrum": {}, "acqu": {}, "proc":{}, "dx":{}} 
    
    # Read in acqu.par and write to dic
    acqupar = os.path.join(dir, acqupar)
    if os.path.isfile(acqupar):
        with open(acqupar, "r") as f:
            info = f.readlines()
        for line in info:
            line = line.replace("\n", "")
            k, v = line.split("=")
            dic["acqu"][k.strip()] = v.strip()

    # Read in proc.par and write to dic
    procpar = os.path.join(dir,procpar)
    if os.path.isfile(procpar):
        with open(procpar, "r") as f:
            info = f.readlines()
        for line in info:
            line = line.replace("\n", "")
            k, v = line.split("=")
            dic["proc"][k.strip()] = v.strip()

    # Define which spectrumfile to take, using 'specfile' when defined, otherwise 
    # the files in 'priority_list' are tried, in that particular order
    priority_list = ["nmr_fid.dx", "data.1d", "fid.1d", "spectrum.1d", "spectrum_processed.1d", None]
    if specfile: 
        inputfile = os.path.join(dir, specfile)
        if not os.path.isfile(inputfile):
            raise IOError("File %s does not exist" % (inputfile))
    else:
        for priority in priority_list:
            if priority == None:
                raise IOError("directory %s does not contain spectral data" % (dir))
            inputfile = os.path.join(dir, priority)
            if os.path.isfile(inputfile):
                break


    # Detect which file we are dealing with from the extension and read in the spectral data
    # Reading .dx file using existing nmrglue.fileio.jcampdx module
    if inputfile.split('.')[-1] == "dx":
        dic["dx"], raw_data = jcampdx.read(inputfile)
        data = np.empty((int(dic["dx"]["$TD"][0]), ), dtype='complex128')
        data = raw_data[0][:] + 1j * raw_data[1][:]
        
    # Reading .1d files
    elif inputfile.split('.')[-1] == "1d":
        with open(inputfile, "rb") as f:
            raw_data = f.read()  
        
        # Write out parameters from the first 32 bytes into dic["spectrum"]
        keys = ["owner", "format", "version", "dataType", "xDim", "yDim", "zDim", "qDim"]
        for i, k in enumerate(keys):
            start = i * 4
            end = start + 4
            value = int.from_bytes( raw_data[start:end], "little")
            dic["spectrum"][k] = value
        data = np.frombuffer(raw_data[end:], "<f")

        # The first 1/3 of the file is xaxis data (s or ppm)
        split = data.shape[-1] // 3
        xscale = data[0 : split]
        dic["spectrum"]["xaxis"] = xscale

        # The rest is real and imaginary data points interleaved
        data = data[split : : 2] + 1j * data[split + 1 : : 2]

    else:
        raise IOError("File %s cannot be interpreted, use .dx or .1d instead" % (inputfile))

    return dic,data



def guess_udic(dic,data):
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

    # Update defalt parameters, first acqu.par parameters in dic are tried, then JCAMP-DX header parameters
    # size
    if data is not None:
        udic[0]["size"] = len(data)
    else:
        warn('No data, cannot set udic size')
    
    # sw
    try:
        udic[0]['sw'] = float(dic['acqu']['bandwidth']) * 1000
    except KeyError:
        try:
            udic[0]['sw'] = float(dic['dx']['$SW'][0]) * float(dic['dx']['$BF1'][0])
        except KeyError:
            try:
                if dic["spectrum"]["freqdata"]:
                    udic[0]['sw'] = dic["spectrum"]["xaxis"][-1] - dic["spectrum"]["xaxis"][0]
                elif data is not None:
                    udic[0]['sw'] = len(data) / dic["spectrum"]["xaxis"][-1]
                else:
                    warn("Cannot set spectral width - set manually using: 'udic[0]['sw'] = x' where x is the spectral width in Hz")
            except KeyError:
                warn("Cannot set spectral width - set manually using: 'udic[0]['sw'] = x' where x is the spectral width in Hz")
    
    # obs
    try:
        udic[0]['obs'] = float(dic['acqu']['b1Freq'])
    except KeyError:
        try:
            udic[0]['obs'] = float(dic['dx']['$BF1'][0])
        except KeyError:
            warn("Cannot set observe frequency - set manually using: 'udic[0]['obs'] = x' where x is magnetic field in MHz")
    
    # car
    try:
        udic[0]['car'] = float(dic['acqu']['lowestFrequency']) + (float(dic['acqu']['bandwidth']) * 1000 / 2)
    except KeyError:
        try:
            udic[0]['car'] = (float(dic['dx']['$REFERENCEPOINT'][0]) * -1 ) + (float(dic['dx']['$SW'][0]) * udic[0]['obs'] / 2)
        except KeyError:
            try:
                udic[0]['car'] = (float(dic['dx']['$BF1'][0]) - float(dic['dx']['$SF'][0])) * 1000000
            except KeyError:
                warn("Cannot set carrier - try: 'udic[0]['car'] = x * udic[0]['obs']' where x is the center of the spectrum in ppm")

    # label
    try:
        udic[0]['label'] = dic['acqu']['rxChannel']
    except KeyError:
        try:
            label_value = dic['dx'][".OBSERVENUCLEUS"][0].replace("^", "")
            udic[0]["label"] = label_value
        except KeyError:
            warn("Cannot set observed nucleus label")
    
    #keys left to default
        # udic[0]['complex']
        # udic[0]['encoding']
        # udic[0]['time'] = True
        # udic[0]['freq'] = False
    return udic
