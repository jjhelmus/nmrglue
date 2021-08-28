"""
Functions for reading and writing NMRPipe files and table (.tab) files
"""

from __future__ import print_function, division

__developer_info__ = """
NMRPipe file structure is described in the NMRPipe man pages and fdatap.h
"""
import io
import struct
import datetime
import os
from typing import Tuple, Union

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from warnings import warn

import numpy as np

from . import fileiobase
from .table import pipe2glue, glue2pipe, guess_pformat

#########################
# table reading/writing #
#########################


def read_table(filename):
    """
    Read a NMRPipe database table (.tab) file.

    Parameters
    ----------
    filename : str
        Filename of NMRPipe table file to read.

    Returns
    -------
    pcomments : list
        List of NMRPipe comment lines
    pformat: list
        List of NMRPipe table column format strings.
    rec : recarray
        Records array with named fields.

    See Also
    --------
    write_table : Write a NMRPipe table file.

    """
    # divide up into comment lines and data lines
    specials = ["VARS", "FORMAT", "NULLSTRING", "NULLVALUE", "REMARK", "DATA"]
    f = open(filename, 'r')
    cl = []
    dl = []
    for line in f:
        for k in specials:
            if line[:len(k)] == k:
                cl.append(line)
                break
        else:
            dl.append(line)
    f.close()

    # pull out and parse the VARS line
    vl = [i for i, l in enumerate(cl) if l[:4] == "VARS"]
    if len(vl) != 1:
        raise IOError("%s has no/more than one VARS line" % (filename))
    dtd = {'names': cl.pop(vl[0]).split()[1:]}

    # pull out and parse the FORMAT line
    fl = [i for i, l in enumerate(cl) if l[:6] == "FORMAT"]
    if len(fl) != 1:
        raise IOError("%s has no/more than one FORMAT line" % (filename))
    pformat = cl.pop(fl[0]).split()[1:]
    p2f = {'d': 'i4', 'f': 'f8', 'e': 'f8', 's': 'S256'}  # pipe -> format
    dtd['formats'] = [p2f[i[-1]] for i in pformat]

    # DEBUG
    # print(dtd['names'],dtd['formats'])
    s = [l.encode('utf-8') for l in dl]

    rec = np.recfromtxt(s, dtype=dtd, comments='XXXXXXXXXXX')
    return cl, pformat, np.atleast_1d(rec)


def write_table(filename, pcomments, pformats, rec, overwrite=False):
    """
    Write a NMRPipe database table (.tab) file.

    Parameters
    ----------
    filename : str
        Filename of file to write to.
    pcomments: list
        List of NMRPipe comment lines.
    pformats :
        List of NMRPipe table column formats strings.
    rec : recarray
        Records array of table.
    overwrite: bool, optional
        True to overwrite file if it exists, False will raise a Warning if the
        file exists.

    See Also
    --------
    read_table : Read a NMRPipe table file.

    """
    if len(rec[0]) != len(pformats):
        s = "number of rec columns %i and pformat elements %i do not match"
        raise ValueError(s % (len(rec[0]), len(pformats)))

    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite)

    # write out the VARS line
    names = rec.dtype.names
    s = "VARS   " + " ".join(names) + "\n"
    f.write(s.encode('utf-8'))

    # write out the FORMAT line
    s = "FORMAT " + " ".join(pformats) + "\n"
    f.write(s.encode('utf-8'))

    # write out any comment lines
    for c in pcomments:
        f.write(c.encode('utf-8'))

    # write out each line of the records array
    s = " ".join(pformats) + "\n"
    for row in rec:
        drow = [i.decode('utf-8') if i.dtype.kind == 'S' else i for i in row]
        f.write((s % tuple(drow)).encode('utf-8'))
    f.close()
    return

###################
# unit conversion #
###################


def make_uc(dic, data, dim=-1):
    """
    Create a unit conversion object

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    dim : int, optional
        Dimension number to create unit conversion object for. Default is for
        last (direct) dimension.

    Returns
    -------
    uc : unit conversion object
        Unit conversion object for given dimension.

    """
    if dim == -1:
        dim = data.ndim - 1     # last dimention

    fn = "FDF" + str(int(dic["FDDIMORDER"][data.ndim - 1 - dim]))
    size = float(data.shape[dim])

    # check for quadrature in indirect dimentions
    if (dic[fn + "QUADFLAG"] != 1) and (dim != data.ndim - 1):
        size = size / 2.
        cplx = True
    else:
        cplx = False

    sw = dic[fn + "SW"]
    if sw == 0.0:
        sw = 1.0
    obs = dic[fn + "OBS"]
    if obs == 0.0:
        obs = 1.0

    # calculate the carrier from the origin, the left most point which has a
    # frequency of CAR*OBS - SW * (N/2 - 1) / 2,
    # see Fig 3.1 on p.36 of Hoch and Stern
    # The carried should have units on MHz so solve the above for CAR*OBS
    orig = dic[fn + "ORIG"]
    car = orig + sw / 2. - sw / size
    return fileiobase.unit_conversion(size, cplx, sw, obs, car)

############################
# dictionary/data creation #
############################

fd2dphase_dic = {"magnitude": 0, "tppi": 1, "states": 2, "image": 3}


def create_data(data):
    """
    Create a NMRPipe data array (recast into float32 or complex64)
    """
    if np.iscomplexobj(data):   # check quadrature
        return np.array(data, dtype="complex64")
    else:
        return np.array(data, dtype="float32")

########################
# universal dictionary #
########################


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.

    """
    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in range(data.ndim):
        udic[i]["size"] = data.shape[i]     # size from data shape

        # determind NMRPipe axis name
        fn = "FDF" + str(int(dic["FDDIMORDER"][data.ndim - 1 - i]))

        # directly corresponding
        udic[i]["sw"] = dic[fn + "SW"]
        udic[i]["obs"] = dic[fn + "OBS"]
        udic[i]["car"] = dic[fn + "CAR"] * dic[fn + "OBS"]  # ppm->hz
        udic[i]["label"] = dic[fn + "LABEL"]

        if dic[fn + "QUADFLAG"] == 1:   # real data
            udic[i]["complex"] = False
        else:
            udic[i]["complex"] = True

        if dic[fn + "FTFLAG"] == 0:     # time domain
            udic[i]["time"] = True
            udic[i]["freq"] = False
        else:
            udic[i]["time"] = False
            udic[i]["freq"] = True

        if i != 0:
            if dic["FD2DPHASE"] == 0:
                udic[i]["encoding"] = "magnitude"
            elif dic["FD2DPHASE"] == 1:
                udic[i]["encoding"] = "tppi"
            elif dic["FD2DPHASE"] == 2:
                udic[i]["encoding"] = "states"
            elif dic["FD2DPHASE"] == 3:
                udic[i]["encoding"] = "image"
            elif dic["FD2DPHASE"] == 4:
                udic[i]["encoding"] = "array"
            else:
                udic[i]["encoding"] = "unknown"
    return udic


def create_dic(udic, datetimeobj=datetime.datetime.now()):
    """
    Crate a NMRPipe parameter dictionary from universal dictionary

    This function does not update the dictionary keys that are unknown such as
    MIN/MAX, apodization and processing parameters, and sizes in none-current
    domain. Also rounding of parameter is different than NMRPipe.

    Parameters
    ----------
    udic : dict
        Universal dictionary of spectral parameters.
    datetimeobj : datetime object, optional
        Datetime to record in NMRPipe dictionary

    Returns
    -------
    dic : dict
        Dictionary NMRPipe parameters.

    """
    # create the blank dictionary
    dic = create_empty_dic()    # create the empty dictionary
    dic = datetime2dic(datetimeobj, dic)  # add the datetime to the dictionary

    # fill global dictionary parameters
    dic["FDDIMCOUNT"] = float(udic["ndim"])

    # FD2DPHASE
    if udic[0]["encoding"] == "tppi":
        dic["FD2DPHASE"] = 1.0
    elif (udic[0]["encoding"] == "complex" or
          udic[0]["encoding"] == "states" or
          udic[0]["encoding"] == "states-tppi"):
        dic["FD2DPHASE"] = 2.0
    else:
        dic["FD2DPHASE"] = 0.0

    # fill in parameters for each dimension
    for i, adic in enumerate([udic[k] for k in range(udic["ndim"])]):
        n = int((dic["FDDIMCOUNT"] - 1) - i)
        dic = add_axis_to_dic(dic, adic, n)

    if dic["FDDIMCOUNT"] >= 3:  # at least 3D
        dic["FDFILECOUNT"] = dic["FDF3SIZE"] * dic["FDF4SIZE"]

    if ((dic["FDF1QUADFLAG"] == dic["FDF2QUADFLAG"] == dic["FDF3QUADFLAG"]) and
            (dic["FDF1QUADFLAG"] == dic["FDF4QUADFLAG"] == 1)):
                dic["FDQUADFLAG"] = 1.0

    return dic


def add_axis_to_dic(dic, adic, n):
    """
    Add an axis dictionary (adic) to a NMRPipe dictionary (dic) as axis n.
    """
    # determind F1,F2,F3,...
    fn = ["FDF2", "FDF1", "FDF3", "FDF4"][n]

    # parameter directly in dictionary
    dic[fn + "SW"] = float(adic["sw"])
    dic[fn + "OBS"] = float(adic["obs"])
    dic[fn + "CAR"] = float(adic["car"] / adic["obs"])
    dic[fn + "LABEL"] = adic["label"]

    if adic["complex"]:
        dic[fn + "QUADFLAG"] = 0.0
    else:
        dic[fn + "QUADFLAG"] = 1.0

    # determine R|I size
    if adic["complex"] and n != 0:
        psize = adic["size"] / 2.
    else:
        psize = adic["size"] / 1.

    # origin calculation size
    osize = psize

    # set FT/TD SIZE and FTFLAG depending on domain
    if adic["time"]:
        dic[fn + "TDSIZE"] = psize
        dic[fn + "FTFLAG"] = 0.0
    else:
        dic[fn + "FTSIZE"] = psize
        dic[fn + "FTFLAG"] = 1.0

    # apodization and center
    dic[fn + "APOD"] = dic[fn + "TDSIZE"]

    if n == 0 or dic["FD2DPHASE"] != 1:
        dic[fn + "CENTER"] = int(psize / 2.) + 1.
    else:   # TPPI requires division by 4
        dic[fn + "CENTER"] = int(psize / 4.) + 1
        osize = psize / 2.

    # origin (last point) is CAR*OBS-SW*(N/2-1)/N
    # see Fig 3.1 on p.36 of Hoch and Stern
    # print("fn:",n)
    # print("CAR:",dic[fn+"CAR"])
    # print("OBS:",dic[fn+"OBS"])
    # print("SW:",dic[fn+"SW"])
    # print("osize:",osize)
    # print("CENTER:",dic[fn+"CENTER"])
    dic[fn + "ORIG"] = (dic[fn + "CAR"] * dic[fn + "OBS"] - dic[fn + "SW"] *
                        (osize - dic[fn + "CENTER"]) / osize)

    if n == 0:  # direct dim
        dic["FDSIZE"] = psize
        dic["FDREALSIZE"] = psize

    if n == 1:  # first indirect
        dic["FDSPECNUM"] = float(adic["size"])  # R+I
        if adic["encoding"] == 'complex':
            dic["FDF1AQSIGN"] = 0
        if adic["encoding"] == 'states':
            dic["FDF1AQSIGN"] = 0  # should this be 2?
        elif adic["encoding"] == 'states-tppi':
            dic["FDF1AQSIGN"] = 16

    if n == 2:  # second indirect
        if adic["complex"]:
            dic["FDF3SIZE"] = psize * 2
        else:
            dic["FDF3SIZE"] = psize

    if n == 3:  # third indirect
        if adic["complex"]:
            dic["FDF4SIZE"] = psize * 2
        else:
            dic["FDF3SIZE"] = psize
    return dic


def create_empty_dic():
    """
    Creates a NMRPipe dictionary with default values
    """
    dic = fdata2dic(np.zeros((512), dtype="float32"))

    # parameters which are 1
    dic["FDF1CENTER"] = 1.
    dic["FDF2CENTER"] = 1.
    dic["FDF3CENTER"] = 1.
    dic["FDF4CENTER"] = 1.

    dic["FDF3SIZE"] = 1.
    dic["FDF4SIZE"] = 1.

    dic["FDF1QUADFLAG"] = 1.
    dic["FDF2QUADFLAG"] = 1.
    dic["FDF3QUADFLAG"] = 1.
    dic["FDF4QUADFLAG"] = 1.

    dic["FDSPECNUM"] = 1.
    dic["FDFILECOUNT"] = 1.
    dic["FD2DVIRGIN"] = 1.
    # dimention ordering

    dic["FDDIMORDER1"] = 2.0
    dic["FDDIMORDER2"] = 1.0
    dic["FDDIMORDER3"] = 3.0
    dic["FDDIMORDER4"] = 4.0
    dic["FDDIMORDER"] = [2.0, 1.0, 3.0, 4.0]

    # string and such
    dic["FDF1LABEL"] = "Y"
    dic["FDF2LABEL"] = "X"
    dic["FDF3LABEL"] = "Z"
    dic["FDF4LABEL"] = "A"

    # misc values
    dic["FDFLTFORMAT"] = struct.unpack('f', b'\xef\xeenO')[0]
    dic["FDFLTORDER"] = float(2.3450000286102295)

    return dic


def datetime2dic(dt, dic):
    """
    Add datatime object to a NMRPipe dictionary
    """
    dic["FDYEAR"] = float(dt.year)
    dic["FDMONTH"] = float(dt.month)
    dic["FDDAY"] = float(dt.day)
    dic["FDHOURS"] = float(dt.hour)
    dic["FDMINS"] = float(dt.minute)
    dic["FDSECS"] = float(dt.second)
    return dic


def dic2datetime(dic):
    """
    Create a datetime object from a NMRPipe dictionary
    """
    year = int(dic["FDYEAR"])
    month = int(dic["FDMONTH"])
    day = int(dic["FDDAY"])
    hour = int(dic["FDHOURS"])
    minute = int(dic["FDMINS"])
    second = int(dic["FDSECS"])
    return datetime.datetime(year, month, day, hour, minute, second)

################
# file reading #
################


def read(filename: Union[str, bytes, io.BytesIO]) -> Tuple[dict, np.array]:
    """
    Read a NMRPipe file.

    For standard multi-file 3D/4D NMRPipe data sets, filename should be a
    filemask (for example "/ft/test%03d.ft3") with a "%" formatter.  If only
    one file of a 3D/4D data set is provided only that 2D slice of the data is
    read (for example "/ft/test001.ft3" results in a 2D data set being read).

    NMRPipe data streams stored as files (one file 3D/4D data sets made using
    xyz2pipe) can be read by providing the file name of the stream.  The entire
    data set is read into memory.

    Parameters
    ----------
    filename : str | bytes | io.BytesIO
        Filename or filemask of NMRPipe file(s) to read. Binary io.BytesIO stream 
        (e.g. from open(filename, "rb")) or bytes buffer can also be provided

    Returns
    --------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.

    See Also
    --------
    read_lowmem : NMRPipe file reading with minimal memory usage.
    write : Write a NMRPipe data to file(s).

    """
    if (type(filename) is bytes):
        filemask = None
    elif hasattr(filename, "read"):
        filename = filename.read()
        filemask = None
    elif filename.count("%") == 1:
        filemask = filename
        filename = filename % 1
    elif filename.count("%") == 2:
        filemask = filename
        filename = filename % (1, 1)
    else:
        filemask = None

    fdata = get_fdata(filename)
    dic = fdata2dic(fdata)
    order = dic["FDDIMCOUNT"]

    if order == 1:
        return read_1D(filename)
    if order == 2:
        return read_2D(filename)
    if dic["FDPIPEFLAG"] != 0:  # open streams
        return read_stream(filename)
    if filemask is None:     # if no filemask open as 2D
        return read_2D(filename)
    if order == 3:
        return read_3D(filemask)
    if order == 4:
        return read_4D(filemask)
    raise ValueError('unknown dimensionality: %s' % order)


def read_lowmem(filename):
    """
    Read a NMRPipe file with minimal memory usage.

    See :py:func:`read` for Parameters and information.

    Returns
    -------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    See Also
    --------
    read : Read NMRPipe files.
    write_lowmem : Write NMRPipe files using minimal amounts of memory.

    """
    if filename.count("%") == 1:
        filemask = filename
        filename = filename % 1
    elif filename.count("%") == 2:
        filemask = filename
        filename = filename % (1, 1)
    else:
        filemask = None

    fdata = get_fdata(filename)
    dic = fdata2dic(fdata)
    order = dic["FDDIMCOUNT"]

    if order == 1:
        return read_1D(filename)    # there is no 1D low memory option
    if order == 2:
        return read_lowmem_2D(filename)
    if dic["FDPIPEFLAG"] != 0:  # open streams
        return read_lowmem_stream(filename)
    if filemask is None:    # if no filemask open as 2D
        return read_lowmem_2D(filename)
    if order == 3:
        return read_lowmem_3D(filemask)
    if order == 4:
        return read_lowmem_4D(filemask)

    raise ValueError('unknown dimentionality: %s' % order)


# dimension specific reading
def read_1D(filename):
    """
    Read a 1D NMRPipe file.

    See :py:func:`read` for documentation.

    """
    fdata, data = get_fdata_data(filename)   # get the fdata and data arrays
    dic = fdata2dic(fdata)  # convert the fdata block to a python dictionary
    data = reshape_data(data, find_shape(dic))    # reshape data

    # unappend imaginary data if needed
    if dic["FDF2QUADFLAG"] != 1:
        data = unappend_data(data)

    return (dic, data)


def read_2D(filename):
    """
    Read a 2D NMRPipe file or NMRPipe data stream.

    See :py:func:`read` for documentation.

    """
    fdata, data = get_fdata_data(filename)   # get the fdata and data arrays
    dic = fdata2dic(fdata)  # convert the fdata block to a python dictionary
    data = reshape_data(data, find_shape(dic))    # reshape data

    # unappend imaginary data if needed
    if dic["FDTRANSPOSED"] == 1 and dic["FDF1QUADFLAG"] != 1:
        data = unappend_data(data)
    elif dic["FDTRANSPOSED"] == 0 and dic["FDF2QUADFLAG"] != 1:
        data = unappend_data(data)

    return (dic, data)


def read_lowmem_2D(filename):
    """
    Read a 2D NMRPipe file or NMRPipe data stream using minimal memory.

    See :py:func:`read_lowmem` for documentation

    """
    dic = fdata2dic(get_fdata(filename))
    order = dic["FDDIMCOUNT"]
    if order == 2:
        data = pipe_2d(filename)
    if order == 3:
        data = pipestream_3d(filename)
    if order == 4:
        data = pipestream_4d(filename)
    return dic, data


def read_stream(filename):
    """
    Read a NMRPipe data stream (one file 3D or 4D files).

    See :py:func:`read` for documentation.

    """
    return read_2D(filename)


def read_lowmem_stream(filename):
    """
    Read a NMRPipe data stream using minimal memory.

    See :py:func:`read_lowmem` for documentation.
    """
    return read_lowmem_2D(filename)


def read_3D(filemask):
    """
    Read a 3D NMRPipe file.

    See :py:func:`read` for documentation.

    """
    dic, data = read_lowmem_3D(filemask)
    data = data[:, :, :]  # read all the data
    return dic, data


def read_lowmem_3D(filemask):
    """
    Read a 3D NMRPipe file using minimal memory.

    See :py:func:`read_lowmem` for documentation

    """
    if '%' not in filemask:  # data streams should be read with read_stream
        return read_lowmem_stream(filemask)
    data = pipe_3d(filemask)    # create a new pipe_3d object
    dic = fdata2dic(get_fdata(filemask % (1)))
    return dic, data


def read_4D(filemask):
    """
    Read a 3D NMRPipe file.

    See :py:func:`read` for documentation.

    Notes
    -----
    This function should not be used to read NMRPipe data streams stored in a
    single file (one file 3D/4D data sets made using xyz2pipe),
    :py:func:`read_2D` should be used.

    """
    dic, data = read_lowmem_4D(filemask)
    data = data[:, :, :, :]  # read all the data
    return dic, data


def read_lowmem_4D(filemask):
    """
    Read a NMRPipe file using minimal memory.

    See :py:func:`read_lowmem` for documentation

    Notes
    -----
    This function should not be used to read NMRPipe data streams stored in a
    single file (one file 3D/4D data sets made using xyz2pipe),
    :py:func:`read_lowmem_2D` should be used.

    """
    if '%' not in filemask:  # data streams should be read with read_stream
        return read_lowmem_stream(filemask)

    data = pipe_4d(filemask)    # create a new pipe_3d object
    if data.singleindex:
        dic = fdata2dic(get_fdata(filemask % (1)))
    else:
        dic = fdata2dic(get_fdata(filemask % (1, 1)))
    return (dic, data)

#####################
# writing functions #
#####################


def write(filename, dic, data, overwrite=False):
    """
    Write a NMRPipe file to disk.

    Parameters
    ----------
    filename : str
        Filename of NMRPipe to write to.  See Notes.
    dic : dict
        Dictionary of NMRPipe parameters.
    data : array_like
        Array of NMR data.
    overwrite : bool, optional.
        Set True to overwrite files, False will raise a Warning if file
        exists.

    Notes
    -----
    For 3D data if filename has no '%' formatter then the data is written as a
    3D NMRPipe data stream.  When the '%' formatter is provided the data is
    written out as a standard NMRPipe 3D multi-file 3D.

    For 4D data, filename can have one, two or no '%' formatters resulting in
    a single index file (test%03d.ft), two index file(test%02d%03d.ft), or
    one file data stream (test.ft4).

    dic["FDPIPEFLAG"] is not changed or checked when writing, please check
    that this value is 0.0 for standard non-data stream files, and 1.0 for data
    stream files or an file may be written with an incorrect header.

    Set overwrite to True to overwrite files that exist.

    See Also
    --------
    write_lowmem : Write NMRPipe files using minimal amounts of memory.
    read : Read NMRPipe files.

    """
    # load all data if the data is not a numpy ndarray
    if not isinstance(data, np.ndarray):
        data = data[:]

    if filename.count("%") == 0:
        return write_single(filename, dic, data, overwrite)
    elif data.ndim == 3:
        return write_3D(filename, dic, data, overwrite)
    elif data.ndim == 4:
        return write_4D(filename, dic, data, overwrite)

    raise ValueError('unknown filename/dimension')


def write_single(filename, dic, data, overwrite=False):
    """
    Write data to a single NMRPipe file from memory.

    Write 1D and 2D files completely as well as NMRPipe data streams.
    2D planes of 3D and 4D files should be written with this function.

    See :py:func:`write` for documentation.

    """
    # append imaginary and flatten
    if data.dtype == "complex64":
        data = append_data(data)
    data = unshape_data(data)

    # create the fdata array
    fdata = dic2fdata(dic)

    # write the file
    put_data(filename, fdata, data, overwrite)
    return


def write_3D(filemask, dic, data, overwrite=False):
    """
    Write a standard multi-file 3D NMRPipe file

    See :py:func:`write` for documentation.

    """
    lenZ, lenY, lenX = data.shape
    for zi in range(lenZ):
        fn = filemask % (zi + 1)
        plane = data[zi]
        write_single(fn, dic, plane, overwrite)
    return


def write_4D(filemask, dic, data, overwrite=False):
    """
    Write a one or two index 4D NMRPipe file.

    See :py:func:`write` for documentation.

    """
    lenA, lenZ, lenY, lenX = data.shape
    for ai in range(lenA):
        for zi in range(lenZ):
            if filemask.count("%") == 2:
                fn = filemask % (ai + 1, zi + 1)
            else:
                fn = filemask % (ai * lenZ + zi + 1)

            plane = data[ai, zi]

            # update dictionary if needed
            if dic["FDSCALEFLAG"] == 1:
                dic["FDMAX"] = plane.max()
                dic["FDDISPMAX"] = dic["FDMAX"]
                dic["FDMIN"] = plane.min()
                dic["FDDISPMIN"] = dic["FDMIN"]
            write_single(fn, dic, plane, overwrite)
    return


def write_lowmem(filename, dic, data, overwrite=False):
    """
    Write a NMRPipe file to disk using minimal memory (trace by trace).

    Parameters
    ----------
    filename : str
        Filename of NMRPipe to write to.  See :py:func:`write` for details.
    dic : dict
        Dictionary of NMRPipe parameters.
    data : array_like
        Array of NMR data.
    overwrite : bool, optional.
        Set True to overwrite files, False will raise a Warning if file
        exists.

    See Also
    --------
    write : Write a NMRPipe file to disk.
    read_lowmem : Read a NMRPipe file using minimal memory.

    """
    if data.ndim == 1:
        return write_single(filename, dic, data, overwrite)
    if data.ndim == 2:
        return write_lowmem_2D(filename, dic, data, overwrite)
    if data.ndim == 3:
        if "%" in filename:
            return write_lowmem_3D(filename, dic, data, overwrite)
        else:
            return write_lowmem_3Ds(filename, dic, data, overwrite)
    if data.ndim == 4:
        if "%" in filename:
            return write_lowmem_4D(filename, dic, data, overwrite)
        else:
            return write_lowmem_4Ds(filename, dic, data, overwrite)

    raise ValueError('unknown dimensionality: %s' % data.ndim)


def write_lowmem_2D(filename, dic, data, overwrite=False):
    """
    Write a 2D NMRPipe file using minimal memory (trace by trace)

    See :py:func:`write_lowmem` for documentation.

    """
    fh = fileiobase.open_towrite(filename, overwrite=overwrite)

    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh, fdata)

    # put data trace by trace
    lenY, lenX = data.shape
    for y in range(lenY):
        put_trace(fh, data[y])
    fh.close()
    return


def write_lowmem_3D(filename, dic, data, overwrite=False):
    """
    Write a standard multi-file 3D NMRPipe file using minimal memory.

    See :py:func:`write_lowmem` for documentation.

    Notes
    -----
    MIN/MAX parameters are not updated in the NMRPipe headers.

    """
    # create the fdata array
    fdata = dic2fdata(dic)

    # put data trace by trace
    lenZ, lenY, lenX = data.shape
    for z in range(lenZ):
        # open the file to store the 2D plane
        fh = fileiobase.open_towrite(filename % (z + 1), overwrite=overwrite)
        put_fdata(fh, fdata)
        for y in range(lenY):
            put_trace(fh, data[z, y])
        fh.close()
    return


def write_lowmem_3Ds(filename, dic, data, overwrite=False):
    """
    Write 3D NMRPipe data stream file using minimal memory (trace by trace)

    See :py:func:`write_lowmem` for documentation.

    """
    fh = fileiobase.open_towrite(filename, overwrite=overwrite)

    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh, fdata)

    # put data trace by trace
    lenZ, lenY, lenX = data.shape
    for z in range(lenZ):
        for y in range(lenY):
            put_trace(fh, data[z, y])
    fh.close()
    return


def write_lowmem_4D(filename, dic, data, overwrite=False):
    """
    Write a multi-file (single or double index) 4D NMRPipe file using
    minimal memory.

    See :py:func:`write_lowmem` for documentation.

    Notes
    -----
    MIN/MAX parameters are not updated in the NMRPipe headers.

    """
    # create the fdata array
    fdata = dic2fdata(dic)

    # put data trace by trace
    lenA, lenZ, lenY, lenX = data.shape
    for a in range(lenA):
        for z in range(lenZ):
            # open the file to store the 2D plane
            if filename.count("%") == 1:
                fname = filename % (a * lenZ + z + 1)
            else:
                fname = filename % (a + 1, z + 1)
            fh = fileiobase.open_towrite(fname, overwrite=overwrite)
            put_fdata(fh, fdata)
            for y in range(lenY):
                put_trace(fh, data[a, z, y])
            fh.close()
    return


def write_lowmem_4Ds(filename, dic, data, overwrite=False):
    """
    Write 4D NMRPipe data stream file using minimal memory (trace by trace)

    See :py:func:`write_lowmem` for documentation.

    """
    fh = fileiobase.open_towrite(filename, overwrite=overwrite)

    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh, fdata)

    # put data trace by trace
    lenA, lenZ, lenY, lenX = data.shape
    for a in range(lenA):
        for z in range(lenZ):
            for y in range(lenY):
                put_trace(fh, data[a, z, y])
    fh.close()
    return


###############
# put to disk #
###############


def put_fdata(fh, fdata):
    """
    Put NMR data, fdata, to a NMRPipe file described by file object fh.
    """
    if fdata.dtype != 'float32':
        raise TypeError('fdata.dtype is not float32')
    fh.write(fdata.tobytes())
    return


def put_trace(fh, trace):
    """
    Put a trace (real or complex) to NMRPipe file described by file object fh.
    """
    if trace.dtype == 'complex64':
        trace = append_data(trace)
    if trace.dtype != 'float32':
        raise TypeError('trace.dtype is not float32')
    fh.write(trace.tobytes())
    return


def put_data(filename, fdata, data, overwrite=False):
    """
    Put fdata and data to 2D NMRPipe.
    """
    if data.dtype != 'float32':
        # print(data.dtype)
        raise TypeError('data.dtype is not float32')
    if fdata.dtype != 'float32':
        raise TypeError('fdata.dtype is not float32')

    # write the file
    f = fileiobase.open_towrite(filename, overwrite=overwrite)
    f.write(fdata.tobytes())
    f.write(data.tobytes())
    f.close()
    return


def write_slice_3D(filemask, dic, data, shape, slices):
    """
    Write a slice of a 3D data array to file.

    Opens (or if necessary creates) a 2D NMRPipe file(s) to write
    data, where the total 3D file size is given by shape.

    Parameters
    ----------
    filemask : str
        String of NMRPipe file with single formatting operator (%).
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        3D array of NMR data.
    shape : tuple
        Tuple of 3 integers indicating the overall matrix shape.
    (sz, sy, sx) : slices
        Slice objects which specify the location of the to be written data.

    Notes
    -----
    This function memmaps 2D NMRPipe files for speed. It only writes
    dictionaries to file when created, leaving them unmodified if the file
    exists. Only error checking is that data is 3D.

    See Also
    --------
    iter3D : Users should use this object, not this function.

    """
    sz, sy, sx = slices
    if data.ndim != 3:
        raise ValueError("passed array must be 3D")

    # unpack the shape
    dz, dy, dx = shape

    # create list of file names
    fnames = [filemask % i for i in range(1, dz + 1)]

    # loop over the requested z-slice
    for i, f in enumerate(fnames[sz]):

        # print("i:",i,"f:",f)
        if os.path.isfile(f) is False:
            # file doesn't exist, create a empty one
            ndata = np.zeros((dy, dx), dtype=data.dtype)
            write_single(f, dic, data, False)
            del(ndata)

        # mmap the [new] file
        mdata = np.memmap(f, dtype='float32', offset=512 * 4, mode='r+')
        # reshape
        mdata = mdata.reshape((dy, dx))

        # unpack into rdata,[idata] depending on quadrature
        if data.dtype == 'complex64':
            h = mdata.shape[-1] // 2.0
            rdata = mdata[..., :h]
            idata = mdata[..., h:]
        else:
            rdata = mdata

        # write the data out, flush and close
        rdata[sy, sx] = data.real[i]
        rdata.flush()
        if data.dtype == 'complex64':
            idata[sy, sx] = data.imag[i]
            idata.flush()
            del(idata)

        # clean up
        del(rdata)
        del(mdata)

# iter3D tools (xyz2pipe and pipe2xyz replacements)
# Notes for iter3D implementation
#
# 'x'/'y' in_lead
# ==============
# Reading
# -------
# - When passed x must transposed 1,2 if dic["FDTRANSPOSED"] == 1
#  (might need to call pipe_proc.tp)
# - if 'y' passed then cann pipe_proc.tp unless dic["FDTRANSPOED"]
# - save 'good' dictionary and return each loop
#
# Looping
# -------
# - will loop until data.shape[0] reached
# - returns dic, XY or YX plane
#
# Writing
# -------
# - if 'y' out then need final pipe_proc.tp of data, if 'x' do nothing
# - reshape data to 1,plane.shape[0],plane.shape[1]
# - size becomes data.shape[0],plane.shape[0],plane.shape[1]
# - sz = slice(i,i+1,1) sy=sx=slice(None)
#
# 'z' in_lead
# ===========
# Reading
# -------
# - Untranspose if dic["TRANSPOSED"] == 1 (call pipe_proc.tp)
# - transpose (1,2,0)
# - ORDER 1,2,3 = 3,1,2 and array
# - update "FDSLICECOUNT" and "FDSIZE" taking into accound complex packing
# - also update "FDSPECNUM"
# - call write_slice3D
# - store shape as self.max_iter
#
# Looping
# -------
# - grab the slice and pack_complex if needed
# - returns dic,ZX-plane
#
# Writing
# -------
# - if out_lead = 'x' needs final pipe_proc.tp of data, if 'z' do nothing
# - reshape data to 1,plane.shape[0],plane.shape[1]
# - transposed data to 2,0,1 (or combine with above step
# - update "FDSIZE" and "FDSPECNUM"
# - remove min/max
# - update FDDIMORDER and ORDER1,2,3
# - size plane.shape[0],self.max_iter,plane.shape[2]
# - sz = slice(None)=sx
# - sy = slice(i,i+1,1)


def pack_complex(data):
    """
    Pack inteleaved real,imag array into complex array.
    """
    return np.array(data[..., ::2] + data[..., 1::2] * 1.j, dtype="complex64")


def transpose_3D(dic, data, axes=(2, 1, 0)):
    """
    Transpose pipe_3d object and dictionary
    """
    a1, a2, a3 = axes
    rdic = dict(dic)    # create a copy of the dictionary

    # transpose the data
    data = data.transpose((a1, a2, a3))

    # transpose the dictionary
    s3 = "FDDIMORDER" + str(int(3 - a1))    # 3rd axis is 0th axis in data_nd
    s2 = "FDDIMORDER" + str(int(3 - a2))    # 2nd axis is 1st axis in data_nd
    s1 = "FDDIMORDER" + str(int(3 - a3))    # 1st axis is 3nd axis in data_nd

    rdic["FDDIMORDER1"] = dic[s1]
    rdic["FDDIMORDER2"] = dic[s2]
    rdic["FDDIMORDER3"] = dic[s3]

    rdic['FDDIMORDER'] = [rdic["FDDIMORDER1"], rdic["FDDIMORDER2"],
                          rdic["FDDIMORDER3"], rdic["FDDIMORDER4"]]

    # set the shape dictionary parameters
    fn = "FDF" + str(int(rdic["FDDIMORDER1"]))
    if rdic[fn + "QUADFLAG"] != 1.0:   # last axis is complex
        rdic["FDSIZE"] = data.shape[2] / 2.
    else:   # last axis is singular
        rdic["FDSIZE"] = data.shape[2]

    rdic["FDSLICECOUNT"] = data.shape[1]
    rdic["FDSPECNUM"] = rdic["FDSLICECOUNT"]
    return rdic, data


class iter3D(object):
    """
    Object which allows for graceful iteration over 3D NMRPipe files.

    iter3D.iter() returns a (dic,plane) tuple which can be written using
    the x.writeplane function.

    When processing 3D files with iter3D object(s) the following dictionary
    parameters may not have the same values as NMRPipe processing scripts
    return:

    * FDSLICECOUNT
    * FDMAX,FDDISMAX,FDMIN,FDDISPMIN when FDSCALEFLAG == 0

    Example::

        #3D data processing
        xiter = iter3D("data/test%03d.fid","x","x")
        for dic,YXplane in xiter:
            # process X and Y axis
            xiter.write("ft/test%03d.ft2",YXplane,dic)
        ziter = iter3D("ft/test%03d.ft2","z","z")
        for dic,XZplane in ziter:
            # process Z axis
            ziter.write("ft/test%03d.ft3",XZplane,dic)

    """
    def __init__(self, filemask, in_lead="x", out_lead="DEFAULT"):
        """
        Create a iter3D object

        Parameters
        ----------
        filemask : str
            String file with single formatter (%) which represents which
            indicates which NMRPipe file to read.
        in_lead : ('x', 'y', 'z'), optional
            Axis name of last (1st) axis in outputted 2D
        out_lead : ('x', 'y', 'z', 'DEFAULT'), optional
            Axis name of axis to be written, typically this is the same as
            in_load, which is the used if 'DEFAULT' is given.

        Notes
        -----
        =======     ===============
        In-lead     Iterated Planes
        =======     ===============
        "x"         ('y','x')
        "y"         ('x','y')
        "z"         ('x','z')
        =======     ===============

        """
        # check for invalid in_lead, out_lead
        if in_lead not in ["x", "y", "z"]:
            raise ValueError("in_lead must be 'x','y' or 'z'")

        if out_lead not in ["x", "y", "z", "DEFAULT"]:
            raise ValueError("out_lead must be 'x','y','z' or 'DEFAULT'")

        if out_lead == "DEFAULT":
            out_lead = in_lead

        if in_lead in ["x", "y"] and out_lead not in ["x", "y"]:
            raise ValueError("Invalid in_lead, out_lead pair")

        if in_lead == "z" and out_lead not in ["x", "z"]:
            raise ValueError("Invalid in_lead, out_lead pair")

        self.in_lead = in_lead
        self.out_lead = out_lead

        self.dic, self.pipe_3d = read_3D(filemask)

        # uptranspose data if needed
        if self.dic["FDTRANSPOSED"] == 1.0:
            # need to switch X and Y (0,2,1)
            self.dic, self.pipe_3d = transpose_3D(self.dic, self.pipe_3d,
                                                  (0, 2, 1))

        # self.pipe_3d and self.dic are now REALLY ZYX order
        # now prep pipe_3d for slicing and make idic the iterator dictionary
        self.i = -1  # counter

        if self.in_lead == "x":
            # leave as is Z(YX)
            self.needs_pack_complex = False
            self.idic = dict(self.dic)
            self.i_max = int(self.pipe_3d.shape[0])
        elif self.in_lead == "y":
            # transpose to Z(XY)
            self.idic, self.pipe_3d = transpose_3D(self.dic, self.pipe_3d,
                                                   (0, 2, 1))
            self.needs_pack_complex = False
            self.i_max = int(self.pipe_3d.shape[0])
        elif self.in_lead == "z":
            # transpose to Y(XZ)
            self.idic, self.pipe_3d = transpose_3D(self.dic, self.pipe_3d,
                                                   (1, 2, 0))
            fn = "FDF" + str(int(self.idic["FDDIMORDER1"]))
            if self.idic[fn + "QUADFLAG"] != 1.0:   # z axis is complex
                self.needs_pack_complex = True
            else:
                self.needs_pack_complex = False
            self.i_max = int(self.pipe_3d.shape[0])
        else:
            raise ValueError("Invalid in_lead")  # this should never be raised.

    def __iter__(self):
        """
        x.__iter__() <==> iter(x)
        """
        return self

    def __next__(self):
        """ next iterator. """
        return self.next()

    def next(self):
        """
        Return the next dic, plane or raise StopIteration
        """
        self.i = self.i + 1
        if self.i >= self.i_max:
            raise StopIteration
        else:
            plane = self.pipe_3d[self.i]
            if self.needs_pack_complex:
                plane = pack_complex(plane)
            return (dict(self.idic), plane)

    def reinitialize(self):
        """
        Restart iterator at first dic,plane.
        """
        self.i = -1

    def write(self, filemask, plane, dic):
        """
        Write out current plane.
        """
        # make the plane a 3D array
        plane = plane.reshape(1, plane.shape[0], plane.shape[1])

        if self.in_lead != self.out_lead:
            # transpose the last two axes
            dic, plane = transpose_3D(dic, plane, (0, 2, 1))

        if self.in_lead == "x" or self.in_lead == "y":
            shape = (self.i_max, plane.shape[1], plane.shape[2])
            sz = slice(self.i, self.i + 1, 1)
            sx = slice(None)
            sy = slice(None)
        elif self.in_lead == "z":
            # reorder from YXZ -> ZYX
            dic, plane = transpose_3D(dic, plane, (2, 0, 1))

            # turn scale flag off
            dic["FDSCALEFLAG"] = 0.0
            # the Y size is incorrect
            dic["FDSPECNUM"] = self.i_max

            # update the file count
            # XXX these could be done better
            dic["FDFILECOUNT"] = plane.shape[0]
            dic["FDF3SIZE"] = plane.shape[0]

            shape = (plane.shape[0], self.i_max, plane.shape[2])
            sx = slice(None)
            sy = slice(self.i, self.i + 1, 1)
            sz = slice(None)
        else:
            raise ValueError("invalid in_lead")  # this should never be raised

        # DEBUGGING
        # print("Writing out slice :",self.i)
        # print("shape:",shape)
        # print("plane.shape",plane.shape)
        # print("sx,sy,sz",sx,sy,sz)
        # print(dic["FDFILECOUNT"])
        write_slice_3D(filemask, dic, plane, shape, (sz, sy, sx))

#####################
# Shaping functions #
#####################


def find_shape(dic):
    """
    Find the shape (tuple) of data in a NMRPipe file from parameters.

    1-tuple is returned for 1D data, 2-tuple for 2D and non-stream 3D/4D data,
    3-tuple or 4-tuple for stream 3D/4D data.

    The last dimension of the tuple is length of the data in the file, the
    actual length of the data matrix may be half of this if the data is
    complex.

    """
    if dic["FDDIMCOUNT"] == 1:  # 1D Data
        if dic["FDF2QUADFLAG"] == 1:
            multi = 1.0
        else:
            multi = 2.0

        dim1 = int(dic["FDSIZE"] * multi)
        return (dim1)
    else:  # 2D+ Data
        if dic["FDF1QUADFLAG"] == 1 and dic["FDTRANSPOSED"] == 1:
            multi = 1.0
        elif dic["FDF2QUADFLAG"] == 1 and dic["FDTRANSPOSED"] == 0:
            multi = 1.0
        else:
            multi = 2.0

        dim1 = int(dic["FDSIZE"] * multi)
        dim2 = int(dic["FDSPECNUM"])

        # when the direct dim is singular and the indirect
        # dim is complex FDSPECNUM is half of the correct value
        if dic["FDQUADFLAG"] == 0 and multi == 1.0:
            dim2 = dim2 * 2

        # check for 3D/4D data stream format files (made using xyz2pipe)
        if dic["FDDIMCOUNT"] == 3 and dic["FDPIPEFLAG"] != 0:
            dim3 = int(dic["FDF3SIZE"])
            return (dim3, dim2, dim1)
        if dic["FDDIMCOUNT"] == 4 and dic["FDPIPEFLAG"] != 0:
            dim3 = int(dic["FDF3SIZE"])
            dim4 = int(dic["FDF4SIZE"])
            return (dim4, dim3, dim2, dim1)

        return (dim2, dim1)


def reshape_data(data, shape):
    """
    Reshape data or return 1D data after warning.
    """
    try:
        return data.reshape(shape)
    except ValueError:
        warn(str(data.shape) + "cannot be shaped into" + str(shape))
        return data


def unshape_data(data):
    """
    Return 1D version of data.
    """
    return data.flatten()


def unappend_data(data):
    """
    Return complex data with last axis (-1) unappended.

    Data should have imaginary data vector appended to real data vector

    """
    h = int(data.shape[-1] / 2)
    return np.array(data[..., :h] + data[..., h:] * 1.j, dtype="complex64")


def append_data(data):
    """
    Return data with last axis (-1) appeneded.

    Data should be complex

    """
    return np.concatenate((data.real, data.imag), axis=-1)

###################
# fdata functions #
###################


def fdata2dic(fdata):
    """
    Convert a fdata array to fdata dictionary.

    Converts the raw 512x4-byte NMRPipe header into a python dictionary
    with keys as given in fdatap.h

    """
    dic = dict()

    # Populate the dictionary with FDATA which contains numbers
    for key in fdata_dic.keys():
        dic[key] = float(fdata[int(fdata_dic[key])])

    # make the FDDIMORDER
    dic["FDDIMORDER"] = [dic["FDDIMORDER1"], dic["FDDIMORDER2"],
                         dic["FDDIMORDER3"], dic["FDDIMORDER4"]]

    def _unpack_str(fmt, d):
        return struct.unpack(fmt, d)[0].decode().strip('\x00')
        #return struct.unpack(fmt, d)[0].decode('utf-8', 'replace')
        #return struct.unpack(fmt, d)[0]

    # Populate the dictionary with FDATA which contains strings
    dic["FDF2LABEL"] = _unpack_str('8s', fdata[16:18])
    dic["FDF1LABEL"] = _unpack_str('8s', fdata[18:20])
    dic["FDF3LABEL"] = _unpack_str('8s', fdata[20:22])
    dic["FDF4LABEL"] = _unpack_str('8s', fdata[22:24])
    dic["FDSRCNAME"] = _unpack_str('16s', fdata[286:290])
    dic["FDUSERNAME"] = _unpack_str('16s', fdata[290:294])
    dic["FDTITLE"] = _unpack_str('60s', fdata[297:312])
    dic["FDCOMMENT"] = _unpack_str('160s', fdata[312:352])
    dic["FDOPERNAME"] = _unpack_str('32s', fdata[464:472])
    return dic


def dic2fdata(dic):
    """
    Converts a NMRPipe dictionary into an array.
    """
    # A 512 4-byte array to hold the nmrPipe header data
    fdata = np.zeros(512, 'float32')

    # Populate the array with the simple numbers
    for key in fdata_nums.keys():
        fdata[int(fdata_dic[key])] = float(dic[key])

    # Check that FDDIMORDER didn't overwrite FDDIMORDER1
    fdata[int(fdata_dic["FDDIMORDER1"])] = dic["FDDIMORDER1"]

    # Pack the various strings into terminated strings of the correct length
    # then into floats in the fdata array
    fdata[16:18] = struct.unpack(
        '2f', struct.pack('8s', dic["FDF2LABEL"].encode()))
    fdata[18:20] = struct.unpack(
        '2f', struct.pack('8s', dic["FDF1LABEL"].encode()))
    fdata[20:22] = struct.unpack(
        '2f', struct.pack('8s', dic["FDF3LABEL"].encode()))
    fdata[22:24] = struct.unpack(
        '2f', struct.pack('8s', dic["FDF4LABEL"].encode()))

    # and the longer strings (typically blank)
    fdata[286:290] = struct.unpack(
        '4f', struct.pack('16s', dic["FDSRCNAME"].encode()))
    fdata[290:294] = struct.unpack(
        '4f', struct.pack('16s', dic["FDUSERNAME"].encode()))
    fdata[297:312] = struct.unpack(
        '15f', struct.pack('60s', dic["FDTITLE"].encode()))
    fdata[312:352] = struct.unpack(
        '40f', struct.pack('160s', dic["FDCOMMENT"].encode()))
    fdata[464:472] = struct.unpack(
        '8f', struct.pack('32s', dic["FDOPERNAME"].encode()))

    return fdata

#################################
# raw reading of data from file #
#################################


def get_fdata(filename: Union[str, bytes]):
    """
    Get an array of length 512-bytes holding NMRPipe header.
    """
    if type(filename) is bytes:
        return get_fdata_bytes(filename)
    else:
        fdata = np.fromfile(filename, 'float32', 512)
        if fdata[2] - 2.345 > 1e-6:    # fdata[2] should be 2.345
            fdata = fdata.byteswap()
        return fdata


def get_data(filename: Union[str, bytes]):
    """
    Get array of data
    """
    if type(filename) is bytes:
        return get_data_bytes(filename)
    else:
        data = np.fromfile(filename, 'float32')
        if data[2] - 2.345 > 1e-6:  # check for byteswap
            data = data.byteswap()
        return data[512:]


def get_fdata_data(filename: Union[str, bytes]):
    """
    Get fdata and data array, return (fdata, data)
    """
    if type(filename) is bytes:
        return get_fdata_data_bytes(filename)
    else:
        data = np.fromfile(filename, 'float32')
        if data[2] - 2.345 > 1e-6:  # check for byteswap
            data = data.byteswap()
        return data[:512], data[512:]


def get_fdata_bytes(data: bytes) -> np.array:
    """
    Get an array of length 512-bytes holding NMRPipe header.
    """
    fdata = np.frombuffer(data, dtype=np.float32, count=512)
    if fdata[2] - 2.345 > 1e-6:    # fdata[2] should be 2.345
        fdata = fdata.byteswap()
    return fdata


def get_data_bytes(data: bytes) -> np.array:
    """
    Get array of data from binary stream
    """
    data = np.frombuffer(data, dtype=np.float32)
    if data[2] - 2.345 > 1e-6:  # check for byteswap
        data = data.byteswap()
    return data[512:]


def get_fdata_data_bytes(data: bytes) -> Tuple[np.array, np.array]:
    """
    Get fdata and data array from binary stream, return (fdata, data)
    """
    data = np.frombuffer(data, dtype=np.float32)
    if data[2] - 2.345 > 1e-6:  # check for byteswap
        data = data.byteswap()
    return data[:512], data[512:]


##############################################
# low memory numpy.ndarray emulating objects #
##############################################


def get_trace(fhandle, ntrace, pts, bswap, cplex):
    """
    Get a single trace from a NMRPipe file

    Parameters
    ----------
    fhandle : file object
        File object of open NMRPipe file.
    ntrace : int
        Trace numbers (starting from 0).
    pts : int
        Number of points in trace, R|I.
    bswap : bool
        True to perform byteswap on trace.
    cplex : bool
        True to unappend imaginary data.

    """
    if cplex:
        tpts = pts * 2  # read twice as many points if data is complex
    else:
        tpts = pts

    fhandle.seek(4 * (512 + ntrace * tpts))  # seek to the start of the trace
    trace = np.fromfile(fhandle, 'float32', tpts)

    if bswap:
        trace = trace.byteswap()
    if cplex:
        return unappend_data(trace)
    else:
        return trace


class pipe_2d(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of 2D NMRPipe files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of 2D NMRPipe file.
    order : tuple
        Ordering of axes against file.

    """

    def __init__(self, filename, order=(0, 1)):
        """
        Create and set up object
        """
        # read and parse the NMRPipe header
        fdata = get_fdata(filename)  # get the header data
        if fdata[2] - 2.345 > 1e-6:  # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # set object attributes
        self.filename = filename
        self.order = order

        # check last axis quadrature
        fn = "FDF" + str(int(dic["FDDIMORDER1"]))
        if dic[fn + "QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype = np.dtype('float32')
        else:
            self.cplex = True
            self.dtype = np.dtype('complex64')
            fshape[1] = fshape[1] // 2

        # finalize
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = pipe_2d(self.filename, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values.

        (sY, sX) is a well formated tuple of slices
        """
        sY, sX = slices
        f = open(self.filename, 'rb')  # open the file for reading

        # determine which objects should be selected
        lenY, lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]

        # create an empty array to store the selected slice
        out = np.empty((len(ych), len(xch)), dtype=self.dtype)

        # read in the data trace by trace
        for yi, y in enumerate(ych):
            ntrace = y
            trace = get_trace(f, ntrace, lenX, self.bswap, self.cplex)
            out[yi] = trace[sX]
        f.close()
        return out

# There are two types of NMRPipe 3D files:
# 1) streams which are single file data sets made with xyz2pipe.
# 2) multiple file data test, names test%03d.ft3, etc.
# Low memory objects exist for both, choose the correct one, or let read
# do it for you.


class pipe_3d(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of 3D NMRPipe files (multiple file data sets).

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filemask : str
        Filename of 3D NMRPipe file. Should contain one formatter '%'
        operator.
    order : tuple
        Ordering of axes against file.
    fcheck : bool, optional.
        True to perform a basic check to see if all files expected for the data
        set exist.  Raises a IOError if files are missing. Default is False.

    """

    def __init__(self, filemask, order=(0, 1, 2), fcheck=False):
        """
        Create and set up object, check that files exist if fcheck is True
        """
        filename = filemask % 1

        # read and parse the NMRPipe header in the first file of the 3D
        fdata = get_fdata(filename)  # get the header data
        if fdata[2] - 2.345 > 1e-6:  # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        # find the shape of the first two dimensions
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))[-2:]

        # find the length of the third dimension
        f3 = "FDF" + str(int(dic["FDDIMORDER3"]))
        quadrature_factor = [2, 1][int(dic[f3 + 'QUADFLAG'])]

        #Checking whether "nmrPipe -fn EXT ..." has been applied to z-dim or not.
        #If EXT has been applied, FDF*XN is not zero.
        #If z-dim is in time-domain, data-size given by FDF*X1 and FDF*XN has to be doubled.
        if dic[f3 + 'FTFLAG']:

            if int(dic[f3 + 'XN']) == 0:
                lenZ = int(dic[f3 + 'FTSIZE'] * quadrature_factor)
            else:
                lenZ = int(dic[f3 + 'XN']) - int(dic[f3 + 'X1']) + 1

        else:
            if int(dic[f3 + 'XN']) == 0:
                lenZ = int(dic[f3 + 'TDSIZE'] * quadrature_factor)
            else:
                lenZ = 2*(int(dic[f3 + 'XN']) - int(dic[f3 + 'X1']) + 1)

        fshape.insert(0, lenZ)   # insert as leading size of fshape

        # check that all files exist if fcheck is set
        if fcheck:
            for i in range(1, lenZ + 1):
                if os.path.exists(filemask % i) is False:
                    raise IOError("File not found: " + str(filemask % i))

        # check last axis quadrature
        fn = "FDF" + str(int(dic["FDDIMORDER1"]))
        if dic[fn + "FTFLAG"] == 1.0:
            self.cplex = False
            self.dtype = np.dtype('float32')
        else:
            self.cplex = True
            self.dtype = np.dtype('complex64')
            fshape[2] = fshape[2] // 2

        # finalize
        self.filemask = filemask
        self.order = order
        self.fshape = fshape
        self.__setdimandshape__()  # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = pipe_3d(self.filemask, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values

        (sZ, sY, sX) is a well formated tuple of slices
        """
        sZ, sY, sX = slices
        # determine which objects should be selected
        lenZ, lenY, lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]

        # create an empty array to store the selected slice
        out = np.empty((len(zch), len(ych), len(xch)), dtype=self.dtype)

        # read in the data file by file and trace by trace
        for zi, z in enumerate(zch):
            # open the Z axis file
            f = open(self.filemask % (z + 1), 'rb')
            for yi, y in enumerate(ych):
                ntrace = y
                trace = get_trace(f, ntrace, lenX, self.bswap, self.cplex)
                out[zi, yi] = trace[sX]
            f.close()
        return out


class pipestream_3d(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of 3D NMRPipe data stream files (one file data sets).

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of 3D NMRPipe stream file.
    order : tuple
        Ordering of axes against file.

    """
    def __init__(self, filename, order=(0, 1, 2)):
        """
        Create and set up object
        """
        # read and parse the NMRPipe header
        fdata = get_fdata(filename)  # get the header data
        if fdata[2] - 2.345 > 1e-6:  # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # check last axis quadrature
        fn = "FDF" + str(int(dic["FDDIMORDER1"]))
        if dic[fn + "QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype = np.dtype('float32')
        else:
            self.cplex = True
            self.dtype = np.dtype('complex64')
            fshape[2] = fshape[2] // 2

        # finalize
        self.filename = filename
        self.order = order
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = pipestream_3d(self.filename, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values

        (sZ, sY, sX) is a well formated tuple of slices
        """
        sZ, sY, sX = slices
        f = open(self.filename, 'rb')  # open the file for reading

        # determine which objects should be selected
        lenZ, lenY, lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]

        # create an empty array to store the selected slice
        out = np.empty((len(zch), len(ych), len(xch)), dtype=self.dtype)

        # read in the data trace by trace
        for zi, z in enumerate(zch):
            for yi, y in enumerate(ych):
                ntrace = y + z * lenY
                trace = get_trace(f, ntrace, lenX, self.bswap, self.cplex)
                out[zi, yi] = trace[sX]
        f.close()
        return out

# There are three types of NMRPipe 4D files:
# 1) streams which are single file data sets made with xyz2pipe.
# 2) single index multiple file data sets, named test%03d.ft4, etc.
# 3) two index muttuple file data sets, named test%02d%03d.ft2, made with
# pipe2xyz and conversion binary.
# Low memory objects exist for all three, choose the correct one, or let read
# do it for you.


class pipe_4d(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of single/two index 4D NMRPipe data files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filemask : str
        Filename of 4D NMRPipe file with one or two formatter (%) operators.
    order : tuple
        Ordering of axes against file.
    fcheck : bool, optional.
        True to perform a basic check to see if all files expected for the data
        set exist.  Raises a IOError if files are missing. Default is False.

    """
    def __init__(self, filemask, order=(0, 1, 2, 3), fcheck=False):
        """
        Create and set up object, check that files exist if fcheck is True
        """
        if filemask.count("%") == 1:
            self.singleindex = True
            filename = filemask % (1)
        elif filemask.count("%") == 2:
            self.singleindex = False
            filename = filemask % (1, 1)
        else:
            raise ValueError("bad filemask")

        # read and parse the NMRPipe header in the first file of the 3D
        fdata = get_fdata(filename)  # get the header data
        if fdata[2] - 2.345 > 1e-6:  # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        # find the shape of the first two dimensions
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))[-2:]

        # find the length of the third dimension
        f3 = "FDF" + str(int(dic["FDDIMORDER3"]))
        quadrature_factor = [2, 1][int(dic[f3 + 'QUADFLAG'])]
        if dic[f3 + 'QUADFLAG']:
            lenZ = int(dic[f3 + 'FTSIZE'] * quadrature_factor)
        else:
            lenZ = int(dic[f3 + 'TDSIZE'] * quadrature_factor)
        fshape.insert(0, lenZ)   # insert as leading size of fshape

        # find the length of the fourth dimension
        f4 = "FDF" + str(int(dic["FDDIMORDER4"]))
        quadrature_factor = [2, 1][int(dic[f4 + 'QUADFLAG'])]
        if dic[f4 + 'QUADFLAG']:
            lenA = int(dic[f4 + 'FTSIZE'] * quadrature_factor)
        else:
            lenA = int(dic[f4 + 'TDSIZE'] * quadrature_factor)
        fshape.insert(0, lenA)   # insert as leading size of fshape

        # check that all files exist if fcheck is set
        if fcheck:
            for ai in range(1, lenA + 1):
                for zi in range(1, lenZ + 1):
                    if self.singleindex:
                        fname = filemask % (ai * lenZ + zi + 1)
                    else:
                        fname = filemask % (ai + 1, zi + 1)
                    if os.path.exists(fname) is False:
                        raise IOError("File not found: " + str(fname))

        # check last axis quadrature
        fn = "FDF" + str(int(dic["FDDIMORDER1"]))
        if dic[fn + "QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype = np.dtype('float32')
        else:
            self.cplex = True
            self.dtype = np.dtype('complex64')
            fshape[3] = fshape[3] // 2

        # finalize
        self.filemask = filemask
        self.order = order
        self.fshape = fshape
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = pipe_4d(self.filemask, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values

        (sZ, sY, sX) is a well formated tuple of slices

        """
        sA, sZ, sY, sX = slices
        # determine which objects should be selected
        lenA, lenZ, lenY, lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]
        ach = range(lenA)[sA]

        # create an empty array to store the selected slice
        out = np.empty((len(ach), len(zch), len(ych), len(xch)),
                       dtype=self.dtype)

        # read in the data file by file, trace by trace
        for ai, a in enumerate(ach):
            for zi, z in enumerate(zch):
                if self.singleindex:   # single index
                    f = open(self.filemask % (a * lenZ + z + 1), 'rb')
                else:   # two index
                    f = open(self.filemask % (a + 1, z + 1), 'rb')
                for yi, y in enumerate(ych):
                    ntrace = y
                    trace = get_trace(f, ntrace, lenX, self.bswap, self.cplex)
                    out[ai, zi, yi] = trace[sX]
                f.close()
        return out


class pipestream_4d(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of 4D NMRPipe data steams (one file 4D data sets).

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of 4D NMRPipe stream file.
    order : tuple
        Ordering of axes against file.

    """

    def __init__(self, filename, order=(0, 1, 2, 3)):
        """
        Create and set up object
        """
        # read and parse the NMRPipe header
        fdata = get_fdata(filename)  # get the header data
        if fdata[2] - 2.345 > 1e-6:  # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # set object attributes
        self.filename = filename
        self.order = order

        # check last axis quadrature
        fn = "FDF" + str(int(dic["FDDIMORDER1"]))
        if dic[fn + "QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype = np.dtype('float32')
        else:
            self.cplex = True
            self.dtype = np.dtype('complex64')
            fshape[3] = fshape[3] // 2

        # finalize
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = pipestream_4d(self.filename, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values

        (sA, sZ, sY, sX) is a well formated tuple of slices

        """
        sA, sZ, sY, sX = slices
        f = open(self.filename, 'rb')  # open the file for reading

        # determine which objects should be selected
        lenA, lenZ, lenY, lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]
        ach = range(lenA)[sA]

        # create an empty array to store the selected slice
        out = np.empty((len(ach), len(zch), len(ych), len(xch)),
                       dtype=self.dtype)

        # read in the data trace by trace
        for ai, a in enumerate(ach):
            for zi, z in enumerate(zch):
                for yi, y in enumerate(ych):
                    ntrace = y + z * lenY + a * lenY * lenZ
                    trace = get_trace(f, ntrace, lenX, self.bswap, self.cplex)
                    out[ai, zi, yi] = trace[sX]
        f.close()
        return out

# data, see fdata.h
fdata_nums = {
    'FDMAGIC': '0',
    'FDFLTFORMAT': '1',
    'FDFLTORDER': '2',

    'FDSIZE': '99',
    'FDREALSIZE': '97',
    'FDSPECNUM': '219',
    'FDQUADFLAG': '106',
    'FD2DPHASE': '256',

    'FDTRANSPOSED': '221',
    'FDDIMCOUNT': '9',

    'FDDIMORDER1': '24',
    'FDDIMORDER2': '25',
    'FDDIMORDER3': '26',
    'FDDIMORDER4': '27',

    'FDNUSDIM': '45',

    'FDPIPEFLAG': '57',
    'FDCUBEFLAG': '447',
    'FDPIPECOUNT': '75',
    'FDSLICECOUNT': '443',  # Also FDSLICECOUNT0
    'FDSLICECOUNT1': '446',
    'FDFILECOUNT': '442',

    'FDTHREADCOUNT': '444',
    'FDTHREADID': '445',

    'FDFIRSTPLANE': '77',
    'FDLASTPLANE': '78',
    'FDPARTITION': '65',

    'FDPLANELOC': '14',

    'FDMAX': '247',
    'FDMIN': '248',
    'FDSCALEFLAG': '250',
    'FDDISPMAX': '251',
    'FDDISPMIN': '252',
    'FDPTHRESH': '253',
    'FDNTHRESH': '254',

    'FDUSER1': '70',
    'FDUSER2': '71',
    'FDUSER3': '72',
    'FDUSER4': '73',
    'FDUSER5': '74',
    'FDUSER6': '76',

    'FDLASTBLOCK': '359',
    'FDCONTBLOCK': '360',
    'FDBASEBLOCK': '361',
    'FDPEAKBLOCK': '362',
    'FDBMAPBLOCK': '363',
    'FDHISTBLOCK': '364',
    'FD1DBLOCK': '365',

    'FDMONTH': '294',
    'FDDAY': '295',
    'FDYEAR': '296',
    'FDHOURS': '283',
    'FDMINS': '284',
    'FDSECS': '285',

    'FDMCFLAG': '135',
    'FDNOISE': '153',
    'FDRANK': '180',
    'FDTEMPERATURE': '157',
    'FDPRESSURE': '158',
    'FD2DVIRGIN': '399',
    'FDTAU': '199',
    'FDDOMINFO': '266',
    'FDMETHINFO': '267',

    'FDSCORE': '370',
    'FDSCANS': '371',

    'FDF2APOD': '95',
    'FDF2SW': '100',
    'FDF2OBS': '119',
    'FDF2OBSMID': '378',
    'FDF2ORIG': '101',
    'FDF2UNITS': '152',
    'FDF2QUADFLAG': '56',
    'FDF2FTFLAG': '220',
    'FDF2AQSIGN': '64',
    'FDF2CAR': '66',
    'FDF2CENTER': '79',
    'FDF2OFFPPM': '480',
    'FDF2P0': '109',
    'FDF2P1': '110',
    'FDF2APODCODE': '413',
    'FDF2APODQ1': '415',
    'FDF2APODQ2': '416',
    'FDF2APODQ3': '417',
    'FDF2LB': '111',
    'FDF2GB': '374',
    'FDF2GOFF': '382',
    'FDF2C1': '418',
    'FDF2APODDF': '419',
    'FDF2ZF': '108',
    'FDF2X1': '257',
    'FDF2XN': '258',
    'FDF2FTSIZE': '96',
    'FDF2TDSIZE': '386',

    'FDDMXVAL': '40',
    'FDDMXFLAG': '41',
    'FDDELTATR': '42',

    'FDF1APOD': '428',
    'FDF1SW': '229',
    'FDF1OBS': '218',
    'FDF1OBSMID': '379',
    'FDF1ORIG': '249',
    'FDF1UNITS': '234',
    'FDF1FTFLAG': '222',
    'FDF1AQSIGN': '475',
    'FDF1QUADFLAG': '55',
    'FDF1CAR': '67',
    'FDF1CENTER': '80',
    'FDF1OFFPPM': '481',
    'FDF1P0': '245',
    'FDF1P1': '246',
    'FDF1APODCODE': '414',
    'FDF1APODQ1': '420',
    'FDF1APODQ2': '421',
    'FDF1APODQ3': '422',
    'FDF1LB': '243',
    'FDF1GB': '375',
    'FDF1GOFF': '383',
    'FDF1C1': '423',
    'FDF1ZF': '437',
    'FDF1X1': '259',
    'FDF1XN': '260',
    'FDF1FTSIZE': '98',
    'FDF1TDSIZE': '387',

    'FDF3APOD': '50',
    'FDF3OBS': '10',
    'FDF3OBSMID': '380',
    'FDF3SW': '11',
    'FDF3ORIG': '12',
    'FDF3FTFLAG': '13',
    'FDF3AQSIGN': '476',
    'FDF3SIZE': '15',
    'FDF3QUADFLAG': '51',
    'FDF3UNITS': '58',
    'FDF3P0': '60',
    'FDF3P1': '61',
    'FDF3CAR': '68',
    'FDF3CENTER': '81',
    'FDF3OFFPPM': '482',
    'FDF3APODCODE': '400',
    'FDF3APODQ1': '401',
    'FDF3APODQ2': '402',
    'FDF3APODQ3': '403',
    'FDF3LB': '372',
    'FDF3GB': '376',
    'FDF3GOFF': '384',
    'FDF3C1': '404',
    'FDF3ZF': '438',
    'FDF3X1': '261',
    'FDF3XN': '262',
    'FDF3FTSIZE': '200',
    'FDF3TDSIZE': '388',

    'FDF4APOD': '53',
    'FDF4OBS': '28',
    'FDF4OBSMID': '381',
    'FDF4SW': '29',
    'FDF4ORIG': '30',
    'FDF4FTFLAG': '31',
    'FDF4AQSIGN': '477',
    'FDF4SIZE': '32',
    'FDF4QUADFLAG': '54',
    'FDF4UNITS': '59',
    'FDF4P0': '62',
    'FDF4P1': '63',
    'FDF4CAR': '69',
    'FDF4CENTER': '82',
    'FDF4OFFPPM': '483',
    'FDF4APODCODE': '405',
    'FDF4APODQ1': '406',
    'FDF4APODQ2': '407',
    'FDF4APODQ3': '408',
    'FDF4LB': '373',
    'FDF4GB': '377',
    'FDF4GOFF': '385',
    'FDF4C1': '409',
    'FDF4ZF': '439',
    'FDF4X1': '263',
    'FDF4XN': '264',
    'FDF4FTSIZE': '201',
    'FDF4TDSIZE': '389',
}

fdata_dic = {
    'FDMAGIC': '0',
    'FDFLTFORMAT': '1',
    'FDFLTORDER': '2',

    'FDSIZE': '99',
    'FDREALSIZE': '97',
    'FDSPECNUM': '219',
    'FDQUADFLAG': '106',
    'FD2DPHASE': '256',

    'FDTRANSPOSED': '221',
    'FDDIMCOUNT': '9',
    'FDDIMORDER': '24',

    'FDDIMORDER1': '24',
    'FDDIMORDER2': '25',
    'FDDIMORDER3': '26',
    'FDDIMORDER4': '27',

    'FDNUSDIM': '45',

    'FDPIPEFLAG': '57',
    'FDCUBEFLAG': '447',
    'FDPIPECOUNT': '75',
    'FDSLICECOUNT': '443',  # Also FDSLICECOUNT0
    'FDSLICECOUNT1': '446',
    'FDFILECOUNT': '442',

    'FDTHREADCOUNT': '444',
    'FDTHREADID': '445',

    'FDFIRSTPLANE': '77',
    'FDLASTPLANE': '78',
    'FDPARTITION': '65',

    'FDPLANELOC': '14',

    'FDMAX': '247',
    'FDMIN': '248',
    'FDSCALEFLAG': '250',
    'FDDISPMAX': '251',
    'FDDISPMIN': '252',
    'FDPTHRESH': '253',
    'FDNTHRESH': '254',

    'FDUSER1': '70',
    'FDUSER2': '71',
    'FDUSER3': '72',
    'FDUSER4': '73',
    'FDUSER5': '74',
    'FDUSER6': '76',

    'FDLASTBLOCK': '359',
    'FDCONTBLOCK': '360',
    'FDBASEBLOCK': '361',
    'FDPEAKBLOCK': '362',
    'FDBMAPBLOCK': '363',
    'FDHISTBLOCK': '364',
    'FD1DBLOCK': '365',

    'FDMONTH': '294',
    'FDDAY': '295',
    'FDYEAR': '296',
    'FDHOURS': '283',
    'FDMINS': '284',
    'FDSECS': '285',

    'FDMCFLAG': '135',
    'FDNOISE': '153',
    'FDRANK': '180',
    'FDTEMPERATURE': '157',
    'FDPRESSURE': '158',
    'FD2DVIRGIN': '399',
    'FDTAU': '199',
    'FDDOMINFO': '266',
    'FDMETHINFO': '267',

    'FDSCORE': '370',
    'FDSCANS': '371',

    'FDSRCNAME': '286',
    'FDUSERNAME': '290',
    'FDOPERNAME': '464',
    'FDTITLE': '297',
    'FDCOMMENT': '312',

    'FDF2LABEL': '16',
    'FDF2APOD': '95',
    'FDF2SW': '100',
    'FDF2OBS': '119',
    'FDF2OBSMID': '378',
    'FDF2ORIG': '101',
    'FDF2UNITS': '152',
    'FDF2QUADFLAG': '56',
    'FDF2FTFLAG': '220',
    'FDF2AQSIGN': '64',
    'FDF2CAR': '66',
    'FDF2CENTER': '79',
    'FDF2OFFPPM': '480',
    'FDF2P0': '109',
    'FDF2P1': '110',
    'FDF2APODCODE': '413',
    'FDF2APODQ1': '415',
    'FDF2APODQ2': '416',
    'FDF2APODQ3': '417',
    'FDF2LB': '111',
    'FDF2GB': '374',
    'FDF2GOFF': '382',
    'FDF2C1': '418',
    'FDF2APODDF': '419',
    'FDF2ZF': '108',
    'FDF2X1': '257',
    'FDF2XN': '258',
    'FDF2FTSIZE': '96',
    'FDF2TDSIZE': '386',

    'FDDMXVAL': '40',
    'FDDMXFLAG': '41',
    'FDDELTATR': '42',

    'FDF1LABEL': '18',
    'FDF1APOD': '428',
    'FDF1SW': '229',
    'FDF1OBS': '218',
    'FDF1OBSMID': '379',
    'FDF1ORIG': '249',
    'FDF1UNITS': '234',
    'FDF1FTFLAG': '222',
    'FDF1AQSIGN': '475',
    'FDF1QUADFLAG': '55',
    'FDF1CAR': '67',
    'FDF1CENTER': '80',
    'FDF1OFFPPM': '481',
    'FDF1P0': '245',
    'FDF1P1': '246',
    'FDF1APODCODE': '414',
    'FDF1APODQ1': '420',
    'FDF1APODQ2': '421',
    'FDF1APODQ3': '422',
    'FDF1LB': '243',
    'FDF1GB': '375',
    'FDF1GOFF': '383',
    'FDF1C1': '423',
    'FDF1ZF': '437',
    'FDF1X1': '259',
    'FDF1XN': '260',
    'FDF1FTSIZE': '98',
    'FDF1TDSIZE': '387',

    'FDF3LABEL': '20',
    'FDF3APOD': '50',
    'FDF3OBS': '10',
    'FDF3OBSMID': '380',
    'FDF3SW': '11',
    'FDF3ORIG': '12',
    'FDF3FTFLAG': '13',
    'FDF3AQSIGN': '476',
    'FDF3SIZE': '15',
    'FDF3QUADFLAG': '51',
    'FDF3UNITS': '58',
    'FDF3P0': '60',
    'FDF3P1': '61',
    'FDF3CAR': '68',
    'FDF3CENTER': '81',
    'FDF3OFFPPM': '482',
    'FDF3APODCODE': '400',
    'FDF3APODQ1': '401',
    'FDF3APODQ2': '402',
    'FDF3APODQ3': '403',
    'FDF3LB': '372',
    'FDF3GB': '376',
    'FDF3GOFF': '384',
    'FDF3C1': '404',
    'FDF3ZF': '438',
    'FDF3X1': '261',
    'FDF3XN': '262',
    'FDF3FTSIZE': '200',
    'FDF3TDSIZE': '388',

    'FDF4LABEL': '22',
    'FDF4APOD': '53',
    'FDF4OBS': '28',
    'FDF4OBSMID': '381',
    'FDF4SW': '29',
    'FDF4ORIG': '30',
    'FDF4FTFLAG': '31',
    'FDF4AQSIGN': '477',
    'FDF4SIZE': '32',
    'FDF4QUADFLAG': '54',
    'FDF4UNITS': '59',
    'FDF4P0': '62',
    'FDF4P1': '63',
    'FDF4CAR': '69',
    'FDF4CENTER': '82',
    'FDF4OFFPPM': '483',
    'FDF4APODCODE': '405',
    'FDF4APODQ1': '406',
    'FDF4APODQ2': '407',
    'FDF4APODQ3': '408',
    'FDF4LB': '373',
    'FDF4GB': '377',
    'FDF4GOFF': '385',
    'FDF4C1': '409',
    'FDF4ZF': '439',
    'FDF4X1': '263',
    'FDF4XN': '264',
    'FDF4FTSIZE': '201',
    'FDF4TDSIZE': '389',
}
