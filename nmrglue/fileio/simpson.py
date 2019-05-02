"""
Functions for reading and writing files created by the `SIMPSON simulations
program <http://www.bionmr.chem.au.dk/bionmr/software/simpson.php>`_.
"""

import math
from warnings import warn

import numpy as np


def read(filename, ftype=None, ndim=None, NP=None, NI=None, spe=None):
    """
    Read a SIMPSON file.

    Read a NMR data file saved using in a number of formats produced using
    the SIMPSON simulation program.

    Parameters
    ----------
    filename : str
        Name of SIMPSON file to read data from.
    ftype : {None, 'TEXT', 'BINARY', 'XREIM', 'XYREIM', 'RAWBIN'}, optional
        A string indicating the type of SIMPSON file:

        ``TEXT`` : SIMPSON text format, no fsave arguments.

        ``BINARY`` : SIMPSON binary format, -binary fsave argument.

        ``XREIM`` : Indexed 1D format, rows of frequency/time, real and
        imaginary parts of the data.  Saved with -xreim argument.

        ``XYREIM`` : Indexed 2D format, rows of frequency/time (indirect
        dimension then direct dimension), real and imaginary parts of the
        data. Saved with -xyreim argument.

        ``RAWBIN`` : Raw binary format.  Saved with -raw_bin argument.  ndim
        and spe must also be provided.  In addition if ndim is 2, NP and NI
        must be defined.

        ``None`` : Automatically determine file type.  If this fails the file
        type should be implicitly provided.

        Other formats of files may be created by SIMPSON, but are not currently
        supported.

    Other Parameters
    ----------------
    ndim : {None, 1, 2}, optional
        Dimensionality of the data in the file, only used when ftype is
        "RAWBIN".
    NP : int, optional
        Number of points {R|I} in the direct dimension.  Only used when
        ftype is "RAWBIN" and ndim is 2.
    NI : int, optional
        Number of points in the indirect dimension.  Only used when ftype is
        "RAWBIN" and ndim is 2.
    spe : bool, optional
        True when the data is in the frequency domain, False for time domain
        data.  Only used when ftype is "RAWBIN"

    Returns
    -------
    dic : dict
        Dictionary of spectra parameters.  For some file formats this may be
        empty.
    data : ndarray
        Complex array of spectral data.

    """
    if ftype is None:
        ftype = guess_ftype(filename)

    if ftype == "TEXT":
        return read_text(filename)
    elif ftype == "BINARY":
        return read_binary(filename)
    elif ftype == "XREIM":
        return read_xreim(filename)
    elif ftype == "XYREIM":
        return read_xyreim(filename)
    elif ftype == "RAWBIN":
        if spe is None:
            raise ValueError('spe must be True or False for raw_bin data')
        if ndim == 1:
            return read_raw_bin_1d(filename, spe)
        elif ndim == 2:
            if NP is None or NI is None:
                raise ValueError("NP and NI must be given for raw_bin data")
            return read_raw_bin_2d(filename, NP, NI, spe)
        else:
            raise ValueError('ndim must be 1 or 2 for raw_bin data')
    else:
        raise ValueError("unknown ftype")


def guess_ftype(filename):
    """ Determine a SIMPSON file type from the first few lines of the file. """
    f = open(filename, 'rb')
    line = f.readline()

    if line[:4] == b"SIMP":   # either text or binary
        f.seek(0)
        for line in f:  # look for format line
            line = line.decode('ascii', 'ignore').strip('\n')
            if line[:6] == "FORMAT":
                f.close()
                return "BINARY"
            elif line == "DATA":    # early end if we screw up
                break
            else:
                continue
        f.close()
        return "TEXT"

    else:   # xyreim, xreim, or binary
        # determine the number of columns in the first 4 lines of the file
        ll1 = len(line.split())
        ll2 = len(f.readline().split())
        ll3 = len(f.readline().split())
        ll4 = len(f.readline().split())
        f.close()
        if ll1 == ll2 == ll3 == ll4 == 3:   # all have 3 columns
            return "XREIM"
        elif ll1 == ll2 == ll3 == ll4 == 4:  # all have 4 columns
            return "XYREIM"
        else:
            return "RAWBIN"

###############################
# Read SIMPSON text file data #
###############################


def read_text(filename):
    """Read a SIMPSON text file. See :py:func:`read`."""
    f = open(filename, 'r')  # open the file

    dic = {}                # initalize dictionary of parameters

    # parse the header of the file, storing parameters in dictionary, stop
    # when we hit the data line
    for line in f:
        line = line.strip('\n')
        if line == "SIMP":
            continue
        elif line == "DATA":
            break
        elif "=" in line:
            key, value = line.split("=")
            dic[key] = value
        else:
            warn("Warning, skipping line: %s" % (line))

    # convert float and int keys
    for key in ['SW1', 'SW']:    # convert keys to floats
        if key in dic:
            dic[key] = float(dic[key])
    for key in ['NP', 'NI', 'NELEM']:    # convert keys to ints
        if key in dic:
            dic[key] = int(dic[key])

    # set NELEM to 1 if its a single 1D/2D dataset
    if 'NELEM' not in dic:
        dic['NELEM'] = 1

    # DEBUGGING
    # return dic

    if "NI" in dic:     # 2D data
        data = np.empty((dic['NI']*dic['NELEM'], dic['NP']), dtype='complex64')

        # loop over remaining lines, extracting data values
        for iline, line in enumerate(f):
            if line[0:3] == "END":   # END marks end of data block
                break
            ni_idx, np_idx = divmod(iline, dic['NP'])
            r_val, i_val = [float(i) for i in line.split()]
            data.real[ni_idx, np_idx] = r_val
            data.imag[ni_idx, np_idx] = i_val

    else:   # 1D data
        data = np.empty((dic['NP']*dic['NELEM']), dtype='complex64')

        # loop over remaining lines, extracting data values
        for iline, line in enumerate(f):
            if line[0:3] == "END":
                break
            r_val, i_val = [float(i) for i in line.split()]
            data.real[iline] = r_val
            data.imag[iline] = i_val
        data = data.reshape(dic['NELEM'], -1)

    f.close()
    return dic, data


#############################
# Read xreim or xyreim data #
#############################


def read_xreim(filename):
    """Read a 1D indexed SIMPSON file.  See :py:func:`read`."""
    f = open(filename, 'r')
    # first pass, determine NP
    for iline, line in enumerate(f):
        pass
    NP = int(iline + 1)

    # second pass, place text-data into data and unit arrays
    f.seek(0)
    data = np.empty((NP, ), dtype='complex64')
    units = np.empty((NP, ), dtype='float32')

    for l_idx, line in enumerate(f):
        np_unit, r_val, i_val = [float(i) for i in line.split()]
        units[l_idx] = np_unit
        data.real[l_idx] = r_val
        data.imag[l_idx] = i_val

    f.close()
    return {'units': units}, data


def read_xyreim(filename):
    """Read a 2D indexed SIMPSON file. See :py:func:`read`."""
    f = open(filename, 'r')

    # first pass, determine NP, NI
    blines = 0  # number of blank lines, gives us num. of indirect points
    for nline, line in enumerate(f):
        if line == '\n':
            if blines == 0:     # location of first blank line gives us NP
                NP = int(nline)
            blines += 1
    NI = blines

    # second pass, place text-data into data and unit arrays
    f.seek(0)
    data = np.empty((NI, NP), dtype='complex64')
    units = np.recarray((NI, NP), dtype=[('ni_unit', 'f8'), ('np_unit', 'f8')])

    for l_idx, line in enumerate(f):
        ni_idx, np_idx = divmod(l_idx, NP + 1)  # determine indicies
        if np_idx == NP:    # skip blank line between blocks
            continue
        # unpack the line and store
        ni_unit, np_unit, r_val, i_val = [float(i) for i in line.split()]
        data.real[ni_idx, np_idx] = r_val
        data.imag[ni_idx, np_idx] = i_val
        units[ni_idx, np_idx].ni_unit = ni_unit
        units[ni_idx, np_idx].np_unit = np_unit

    f.close()
    return {'units': units}, data


######################
# Read raw_bin files #
######################


def read_raw_bin_1d(filename, spe=False):
    """Read a 1D raw binary SIMPSON file. See :py:func:`read`. """
    data = unappend_data(np.fromfile(filename, dtype='float32'))
    if spe:
        return {}, data[::-1]   # freq domain data is flipped
    else:
        return {}, data


def read_raw_bin_2d(filename, NP, NI, spe=False):
    """Read a 2D raw binary SIMPSON file. See :py:func:`read`. """
    data = np.fromfile(filename, dtype='float32').reshape(NI, NP * 2)
    if spe:
        # freq domain is saved reversed and shifted by 1 pt in each dimension,
        # so the first points is 0 in each vector, then N-1, N-2, ... 1
        # A bug? in SIMPSON causes the direct dimension shift to be left off
        # the imaginary data.  In addition the first real point in each vector
        # is the first real point in the spectrum.  Therefore only data[:,1:]
        # will match with other file formats.
        cdata = np.empty((NI, NP), dtype='complex64')
        cdata.real = np.roll(np.roll(data[::-1, NP - 1::-1], 1, 1), 1, 0)
        cdata.imag = np.roll(data[::-1, 2 * NP:NP - 1:-1], 1, 0)
        return {}, cdata
    else:
        # time domain is written as vectors of NP, with the imag appended
        return {}, unappend_data(data)


def unappend_data(data):
    """
    Return complex data with last axis unappended.

    Data should have imaginary data vector appended to real data vector.
    """
    h = int(data.shape[-1] / 2.0)
    return np.array(data[..., :h] + data[..., h:] * 1.j, dtype='complex64')


#############################
# Read SIMPSON binary files #
#############################


def read_binary(filename):
    """Read a binary SIMPSON file. See :py:func:`read`."""
    f = open(filename, 'r')  # open the file

    dic = {}                # initalize dictionary of parameters

    # parse the header of the file, storing parameters in dictionary, stop
    # when we hit the data line
    for line in f:
        line = line.strip('\n')
        if line == "SIMP":
            continue
        elif line == "DATA":
            break
        elif "=" in line:
            key, value = line.split("=")
            dic[key] = value
        else:
            warn("Warning, skipping line: %s" % (line))

    # convert float and int keys
    for key in ['SW1', 'SW']:    # convert keys to floats
        if key in dic:
            dic[key] = float(dic[key])
    for key in ['NP', 'NI', 'NELEM']:    # convert keys to ints
        if key in dic:
            dic[key] = int(dic[key])

    if not 'NELEM' in dic.keys():
        dic['NELEM'] = 1

    # DEBUGGING
    # return dic, f

    # extract characters from data block
    chardata = "".join([line.strip('\n') for line in f])
    chardata = chardata[:-3]    # remove END
    f.close()

    # convert every 4 character to 3 'bytes'
    nquads, mod = divmod(len(chardata), 4)
    assert mod == 0     # character should be in blocks of 4
    bytes = []
    for i in range(nquads):
        bytes += chars2bytes(chardata[i * 4:(i + 1) * 4])

    # DEBUGGING
    # return dic, bytes
    # convert every 4 'bytes' to a float, then complex data
    num_points, _num_pad = divmod(len(bytes), 4)
    data = np.empty((num_points, ), dtype='float32')

    for i in range(num_points):
        data[i] = bytes2float(bytes[i * 4: (i + 1) * 4])
    data = data.view('complex64')

    # reorder data according to dimensionality and domain
    if 'NI' in dic:  # 2D data
        return dic, data.reshape(dic['NI']*dic['NELEM'], dic['NP'])

    else:   # 1D data
        return dic, data.reshape(dic['NELEM'], dic['NP'])

BASE = 33
FIRST = lambda f, x: ((x) & ~(~0 << f))
LAST = lambda f, x: ((x) & (~0 << (8 - f)))


def chars2bytes(chars):
    """Convert four characters from a data block into 3 'bytes'."""
    c0, c1, c2, c3 = [ord(c) - BASE for c in chars]

    return [FIRST(6, c0) | LAST(2, c1 << 2),
            FIRST(4, c1) | LAST(4, c2 << 2),
            FIRST(2, c2) | LAST(6, c3 << 2)]


def bytes2float(bytes):
    """ Convert four bytes to a string. """
    # the first 23 bits (0-22) define the mantissa
    # the next 8 bits (23-30) the exponent,
    # the last bit (31) defines the sign
    b0, b1, b2, b3 = bytes
    mantissa = ((b2 % 128) << 16) + (b1 << 8) + b0
    exponent = (b3 % 128) * 2 + (b2 >= 128) * 1
    negative = b3 >= 128

    e = exponent - 0x7f
    m = np.abs(mantissa) / np.float64(1 << 23)

    if negative:
        return -math.ldexp(m, e)
    return math.ldexp(m, e)
