"""
Functions for reading and writing Agilent/Varian binary (fid) files and
parameter (procpar) files.
"""


__developer_doc__ = """
Agilent/Varian file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Both the Agilent/Varian binary and parameter formats are documented in:

* Varian MR News 2005-04-18 Importing Binary VnmrJ / VNMR FIDs into Third
    Party Software and VnmrJ / VNMR FID Data Format
* VnmrJ User Programming - Chapter 5: Parameters and Data

These are available (as of 04/2011) online from
`Agilent <http://agilent.com>`_.

"""

import os
import struct
import inspect
from warnings import warn
from functools import reduce
import operator

import numpy as np

from . import fileiobase

############################
# dictionary/data creation #
############################


def create_data(data):
    """
    Create a Agilent/Varian data array (recast into complex64 array)
    """
    return np.array(data, dtype="complex64")

########################
# universal dictionary #
########################


def guess_udic(dic, data):
    """
    Guess parameter of a universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of Agilent/Varian parameters.
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
        udic[i]["size"] = data.shape[i]

    return udic


def create_dic(udic):
    """
    Create a Agilent/Varian parameter dictionary from a universal dictionary.

    Parameters
    ----------
    udic : dict
        Universal dictionary of spectral parameters.

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian parameters
    """
    ddic = udic[udic["ndim"] - 1]   # direct dimension dictionary

    # number of points in direct dimension is R+I
    if ddic["complex"]:
        num_points = ddic["size"] * 2
    else:
        num_points = ddic["size"]

    # number of blocks is product of R+I of all other dimensions
    num_blocks = 1
    for f in [udic[k]["size"] for k in range(udic["ndim"] - 1)]:
        num_blocks = num_blocks * f

    # create the varian dictionary
    dic = dict()

    # default flag bits
    dic["S_DATA"] = 1
    dic["S_SPEC"] = 0
    dic["S_32"] = 0
    dic["S_FLOAT"] = 1
    dic["S_COMPLEX"] = 0
    dic["S_HYPERCOMPLEX"] = 0
    dic["S_ACQPAR"] = 1
    dic["S_SECND"] = 0
    dic["S_TRANSF"] = 0
    dic["S_NP"] = 0
    dic["S_NF"] = 0
    dic["S_NI"] = 0
    dic["S_NI2"] = 0

    # parameters
    dic["nblocks"] = num_blocks
    dic["ntraces"] = 1
    dic["np"] = num_points
    dic["ebytes"] = 4
    dic["tbytes"] = 4 * num_points
    dic["bbytes"] = dic["tbytes"] + 28
    dic["vers_id"] = 0
    dic["status"] = 201     # 6th bit set
    dic["nbheaders"] = 1

    # fake procpar dictionary
    pdic = dict()

    # direct dimension
    ndim = udic["ndim"]
    pdic["np"] = create_pdic_param('np', [str(udic[ndim - 1]['size'] * 2)])

    if ndim >= 2:
        if udic[ndim - 2]['complex']:
            pdic["ni"] = create_pdic_param('ni',
                                           [str(udic[ndim - 2]['size'] // 2)])
            pdic["phase"] = create_pdic_param('phase', ['0', '1'])
        else:
            pdic["ni"] = create_pdic_param('ni', [str(udic[ndim - 2]['size'])])
            pdic["phase"] = create_pdic_param('phase', ['0'])

    if ndim >= 3:
        pdic['array'] = create_pdic_param('array', ['phase,phase2'])
        if udic[ndim - 3]['complex']:
            pdic["ni2"] = create_pdic_param(
                'ni2', [str(udic[ndim - 3]['size'] // 2)])
            pdic["phase2"] = create_pdic_param('phase', ['0', '1'])
        else:
            pdic["ni2"] = create_pdic_param('ni2',
                                            [str(udic[ndim - 3]['size'])])
            pdic["phase2"] = create_pdic_param('phase2', ['0'])

    if ndim >= 4:
        pdic['array'] = create_pdic_param('array', ['phase,phase2,phase3'])
        if udic[ndim - 4]['complex']:
            pdic["ni3"] = create_pdic_param(
                'ni3', [str(udic[ndim - 4]['size'] // 2)])
            pdic["phase3"] = create_pdic_param('phase3', ['0', '1'])
        else:
            pdic["ni3"] = create_pdic_param(
                'ni3', [str(udic[ndim - 4]['size'])])
            pdic["phase3"] = create_pdic_param('phase3', ['0'])

    dic["procpar"] = pdic
    return dic


def create_pdic_param(name, values):
    """
    Create a fake procpar dictionary element of with given name and values.
    """
    dic = dict()
    dic["Dgroup"] = '1'
    dic["Ggroup"] = '2'
    dic["active"] = '1'
    dic["basictype"] = '1'
    dic["enumerable"] = '0'
    dic["intptr"] = '64'
    dic["maxvalue"] = '32767'   # 2^15-1
    dic["minvalue"] = '0'
    dic["name"] = name
    dic["protection"] = '0'
    dic["stepsize"] = '0'
    dic["subtype"] = '7'
    dic["values"] = values
    return dic

########################
# file reading/writing #
########################


def read(dir=".", fid_file="fid", procpar_file="procpar", read_blockhead=False,
         shape=None, torder=None, as_2d=False):
    """
    Read Agilent/Varian files in a directory.

    Parameters
    ----------
    dir : str, optional
        Directory holding Agilent/Varian data. Default is the current working
        directory.
    fid_file : str, optional
        Filename of binary (fid) file in directory.
    procpar_file : str, optional
        Filename of procpar file in directory.
    read_blockhead : bool, optional
        True to read blockheader(s) and return then in the Agilent/Varian
        parameter dictionary. False (default) does not perform this reading.
    shape : tuple of ints, optional
        Shape of data in binary file. None (default)  will attempt to finds
        this automatically.
    torder : {None, 'r', 'o', 'f' or a Python function} , optional
        Description of the mapping of traces in the file to the NMR data
        matrix. None (the default) will attempt to find this automatically
        which is typically fine for most NMR experiments. See below for
        additional details.
    as_2d : bool, optional
        True to return data as a 2D array ignoring the shape and torder
        parameters.

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian parameters.
    data : ndarray
        Array of NMR data

    Notes
    -----
    The torder parameter describes how the traces on disk should be
    re-organized to form the NMR data matrix.  In most cases this can be
    determined automatically by examining the order of phase parameters in the
    procpar file.  This is done if torder is set to None.  In some cases
    it must be explicitly provided.  Three common cases are:

    ======== ======================  ================= =======================
    Name     Ordering of arrays      torder            Notes
    ======== ======================  ================= =======================
    regular  d3, d2, phase2, phase   'regular' or 'r'
    opposite d3, d2, phase, phase2   'opposite' or 'o'
    flat     As data exists in file  'flat' or 'f'     Valid for 1D or 2D data
    ======== ======================  ================= =======================

    In addition a function which maps indirect dimension index tuples to/from
    trace numbers as stored on disk can be provided.  For reading this
    function should take 2 arguments: shape, index_tuple and return an integer
    trace number.  For writing this function should again take 2 arguments:
    shape,trace_number and return the indirect dimension index tuple for the
    given trace.

    See Also
    --------
    read_lowmem : Read Agilent/Varian files using minimal amounts of memory.
    write : Write Agilent/Varian files.

    """
    if os.path.isdir(dir) is False:
        raise OSError("directory %s does not exist" % (dir))

    # read in the procpar file
    pdic = read_procpar(os.path.join(dir, procpar_file))

    # determine shape and torder is needed
    if as_2d is False:
        if shape is None:
            shape = find_shape(pdic)
        if torder is None and shape is not None and len(shape) >= 3:
            torder = find_torder(pdic, shape)

    # read in the fid file
    fname = os.path.join(dir, fid_file)
    dic, data = read_fid(fname, shape, torder, as_2d, read_blockhead)

    # add the procpar dictionary to the main dictionary
    dic["procpar"] = pdic

    return dic, data


def read_lowmem(dir=".", fid_file="fid", procpar_file="procpar",
                read_blockhead=False, shape=None, torder=None):
    """
    Read Agilent/Varian files in a directory using minimal amount of memory.

    Parameters
    ----------
    dir : str, optional
        Directory holding Agilent/Varian data. Default is the current working
        directory.
    fid_file : str, optional
        Filename of binary (fid) file in directory.
    procpar_file : str, optional
        Filename of procpar file in directory.
    read_blockhead : bool, optional
        True to read blockheader(s) and return then in the Agilent/Varian
        parameter dictionary. False (default) does not perform this reading.
    shape : tuple of ints, optional
        Shape of data in binary file. None (default)  will attempt to finds
        this automatically.
    torder : {None, 'r', 'o', 'f' or a Python function} , optional
        Description of the mapping of traces in the file to the NMR data
        matrix. None (the default) will attempt to find this automatically
        which is typically fine for most NMR experiments. See :py:func:`read`
        for additional details.

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    See Also
    --------
    read : Read Agilent/Varian files.
    write_lowmem : Write Agilent/Varian files using minimal memory

    """
    if os.path.isdir(dir) is False:
        raise OSError("directory %s does not exist" % (dir))

    # read in the procpar file
    pdic = read_procpar(os.path.join(dir, procpar_file))

    if shape is None:
        shape = find_shape(pdic)

    if torder is None:    # always try to find trace order
        torder = find_torder(pdic, shape)

    # read in the fid file
    fname = os.path.join(dir, fid_file)
    dic, data = read_fid_lowmem(fname, shape, torder, read_blockhead)

    # add the procpar dictionary to the main dictionary
    dic["procpar"] = pdic

    return dic, data


def write(dir, dic, data, fid_file="fid", procpar_file="procpar",
          torder=None, repack=False, overwrite=False):
    """
    Write Agilent/Varian files to a directory.

    Parameters
    ----------
    dir : str
        Name of the directory to write to.
    dic : dict
        Dictionary of Agilent/Varian parameters.
    data : array_like
        Array of NMR data to write.
    fid_file : str, optional
        Filename of binary (fid) file in directory to write to.
    procpar_file : str, optional
        Filename of procpar file in directory to write to.
    torder : {None, 'r', 'o', 'f' or a Python function} , optional
        Description of the mapping of traces in the file to the NMR data
        matrix. None (the default) will attempt to find this automatically
        which is typically fine for most NMR experiments. See :py:func:`read`
        for additional details.
    repack : bool, optional
        True to repack file and block headers. False (default) leaves these as
        is.
    overwrite : bool, optional
        Set to True to overwrite existing files, False will raise a Warning if
        files already exist.

    See Also
    --------
    write_lowmem : Write Agilent/Varian files using minimal memory
    read : Read Agilent/Varian files.

    """
    if torder is None and data.ndim >= 3:
        torder = find_torder(dic["procpar"], data.shape)

    # write out the fid file
    fname = os.path.join(dir, fid_file)
    write_fid(fname, dic, data, torder, repack=repack, overwrite=overwrite)

    # write out procpar file
    write_procpar(os.path.join(dir, procpar_file), dic["procpar"], overwrite)


def write_lowmem(dir, dic, data, fid_file="fid", procpar_file="procpar",
                 torder=None, repack=False, overwrite=False):
    """
    Write Agilent/Varian files to a directory using minimal amounts of memory.

    Parameters
    ----------
    dir : str
        Name of the directory to write to.
    dic : dict
        Dictionary of Agilent/Varian parameters.
    data : array_like
        Array of NMR data to write.
    fid_file : str, optional
        Filename of binary (fid) file in directory to write to.
    procpar_file : str, optional
        Filename of procpar file in directory to write to.
    torder : {None, 'r', 'o', 'f' or a Python function} , optional
        Description of the mapping of traces in the file to the NMR data
        matrix. None (the default) will attempt to find this automatically
        which is typically fine for most NMR experiments. See :py:func:`read`
        for additional details.
    repack : bool, optional
        True to repack file and block headers. False (default) leaves these as
        is.
    overwrite : bool, optional
        Set to True to overwrite existing files, False will raise a Warning if
        files already exist.

    See Also
    --------
    write : Write Agilent/Varian files.
    read_lowmem : Read Agilent/Varian using minimal amounts of memory.

    """
    # always find trace ordering
    if torder is None:
        torder = find_torder(dic["procpar"], data.shape)

    # write out the fid file
    fname = os.path.join(dir, fid_file)
    write_fid_lowmem(fname, dic, data, torder, repack=repack,
                     overwrite=overwrite)

    # write out procpar file
    write_procpar(os.path.join(dir, procpar_file), dic["procpar"], overwrite)

############
# ordering #
############


def find_torder(dic, shape):
    """
    Find the torder from the procpar dictionary.

    If procpar dictionary is incomplete a UserWarning is issued and 'r' is
    returned.

    Parameters
    ----------
    dic : dict
        Dictionary of parameters in the procpar file.
    shape : tuple of ints
        Shape of NMR data.

    Returns
    --------
    torder : {'r', 'f', 'o'}
        File ordering for using in :py:func:`read` or :py:func:`write`.

    """
    ndim = len(shape)

    # 1 and 2D files are flat
    if ndim < 3:
        return 'f'  # flat

    if "array" not in dic:
        warn("no array in dictionary, torder set to 'r'")
        return 'r'

    # extract the array list
    al = dic['array']['values'][0].split(',')

    # remove one dimension if non-phase parameter is present in array list
    if False in [s.startswith('phase') for s in al]:
        ndim = ndim - 1
        if ndim < 3:
            return 'f'  # flat

    if ndim == 3:
        if "phase" in al and "phase2" in al:
            if al.index("phase") > al.index("phase2"):
                return 'r'  # regular
            else:
                return 'o'  # opposite
        else:
            warn("missing phase order, torder set to 'r'")
            return 'r'

    if ndim == 4:
        if "phase" in al and "phase2" in al and "phase3" in al:
            if al.index("phase") > al.index("phase2") > al.index("phase3"):
                return 'r'  # regular
            if al.index("phase") < al.index("phase2") < al.index("phase3"):
                return 'o'  # opposite
            else:
                warn("bad phase order, torder set to 'r'")
                return 'r'
        else:
            warn("missing phase order, torder set to 'r'")
            return 'r'

    warn("No trace ordering for" + str(ndim) +
         "dimensional data, torder set to 'r'")
    return 'r'


def torder2i2t(torder):
    """
    Convert torder to an index2trace function
    """
    # if torder is a function, return it
    if inspect.isfunction(torder):
        return torder
    elif torder in ('flat', 'f'):
        return fileiobase.index2trace_flat
    elif torder in ('opposite', 'o'):
        return fileiobase.index2trace_opp
    elif torder in ('regular', 'r'):
        return fileiobase.index2trace_reg
    else:
        raise ValueError("unknown torder" + str(torder))


def torder2t2i(torder):
    """
    Convert torder to a trace2index functions
    """
    if inspect.isfunction(torder):
        return torder
    elif torder in ('flat', 'f'):
        return fileiobase.trace2index_flat
    elif torder in ('opposite', 'o'):
        return fileiobase.trace2index_opp
    elif torder in ('regular', 'r'):
        return fileiobase.trace2index_reg
    else:
        raise ValueError("unknown torder" + str(torder))


def reorder_data(data, shape, torder):
    """
    Reorder raw data from file.

    Parameters
    ----------
    data : 2D ndarray
        Raw data as ordered in binary file.
    shape : tuple of ints
        Shape of the NMR data.
    torder : {'f', 'r', 'o' or a Python function}
        Trace ordering . See :py:func:`read` for details.

    Returns
    -------
    tdata : ndarray
        Array of NMR data.

    Notes
    -----
    Minimal error checking is done to see if data and shape contain the same
    number of values.  This should be done before calling this function.

    """
    # take care of flat files...
    if torder in ('flat', 'f'):
        try:
            data = data.reshape(shape)
        except ValueError:
            warn(str(data.shape) + "cannot be shaped into" + str(shape))
        return data

    # all other cases
    # make an empty array to hold reordered data
    ndata = np.empty(shape, dtype=data.dtype)

    # index2tuple converter
    i2t = torder2i2t(torder)

    # loop over all non-direct dimension indices
    for tup in np.ndindex(shape[:-1]):
        # determine the corresponding trace
        ntrace = i2t(shape[:-1], tup)
        # write the trace to the index
        ndata[tup] = data[ntrace]
    return ndata


def order_data(data, torder):
    """
    Order NMR data for writing to file.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    torder : {'f', 'r', 'o'}
        Trace ordering.  See :py:func:`read` for details.

    Returns
    -------
    raw_data : 2D ndarray
        2D array ordered for writing to Agilent/Varian binary file.

    """
    # determine the shape of the on disk data matrix
    ntraces = reduce(operator.mul, data.shape[:-1])
    nshape = (ntraces, data.shape[-1])

    # take care of flat files
    if torder in ('flat', 'f'):
        return data.reshape(nshape)

    # make an emprt array to hold the 2D disk formatted data matrix
    ndata = np.empty(nshape, dtype=data.dtype)

    # index2tuple converter
    t2i = torder2t2i(torder)

    # loop over all non-direct dimension indices
    for ntrace in range(ntraces):
        tup = t2i(data.shape[:-1], ntrace)
        ndata[ntrace] = data[tup]
    return ndata


#######################
# fid reading/writing #
#######################


def read_fid(filename, shape=None, torder='flat', as_2d=False,
             read_blockhead=False):
    """
    Read a Agilent/Varian binary (fid) file.

    Parameters
    ----------

    filename : str
        Filename of Agilent/Varian binary file (fid) to read.
    shape : tuple of ints, optional
        Shape of the binary data. If not provided data is returned as a 2D
        array. Required if more than one trace per block (non-standard).
    torder : {'f', 'n', 'o'}
        Trace order. See :py:func:`read` for details.
    as_2d : bool, optional
        True to return the data as a 2D array, ignoring the shape and torder
        parameters.
    read_blockhead : bool, optional
        True to read the Agilent/Varian blockheaders(s) into the returned
        dictionary. False ignores them.

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian binary file parameters.
    data : ndarray
        Array of NMR data.

    See Also
    --------
    read_fid_lowmem : Read a Agilent/Varian binary file using minimal amounts
        of memory.
    read : Read Agilent/Varian files from a directory.

    """
    # open the file
    f = open(filename, 'rb')

    # read the fileheader
    dic = fileheader2dic(get_fileheader(f))

    # if ntraces is not 1 use _ntraces version
    if dic["ntraces"] != 1:
        return read_fid_ntraces(filename, shape, torder, as_2d, read_blockhead)

    # data parameters
    dt = find_dtype(dic)
    nblocks = dic["nblocks"]
    pts = dic["np"]
    nbheaders = dic["nbheaders"]

    # read the data
    if read_blockhead:
        bdic, data = get_nblocks(f, nblocks, pts, nbheaders, dt,
                                 read_blockhead)
        dic["blockheader"] = bdic
    else:
        data = get_nblocks(f, nblocks, pts, nbheaders, dt, read_blockhead)
    f.close()

    # uninterleave the real and imaginary data
    data = uninterleave_data(data)

    # return as 2D is requested
    if as_2d:
        return dic, data

    # return raw data if no shape provided.
    if shape is None:
        warn("unknown shape, returning unshaped data")
        return dic, data

    # check for 1D
    if data.shape[0] == 1:
        return dic, np.squeeze(data)

    # try to reshape
    # reorder 3D/4D data
    if len(shape) >= 3:
        try:
            return dic, reorder_data(data, shape, torder)
        except:
            warn("data cannot be re-ordered, returning raw 2D data\n" +
                 "Provided shape: " + str(shape) + " torder: " + str(torder))
            return dic, data

    try:
        data = data.reshape(shape)
    except ValueError:
        warn(str(data.shape) + "cannot be shaped into" + str(shape))
        return dic, data

    return dic, data


def read_fid_lowmem(filename, shape=None, torder='flat', as_2d=False,
                    read_blockhead=False):
    """
    Read a Agilent/Varian binary (fid) file using minimal amounts of memory.

    Parameters
    ----------
    filename : str
        Filename of Agilent/Varian binary file (fid) to read.
    shape : tuple of ints, optional
        Shape of the binary data. If not provided data is returned as a 2D
        array. Required if more than one trace per block (non-standard).
    torder : {'f', 'n', 'o'}
        Trace order. See :py:func:`read` for details.
    read_blockhead : bool, optional
        Not used, retained for compatibility with :py:func:`read_fid`

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian binary file parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    See Also
    --------
    read_fid : Read a Agilent/Varian binary file.
    read_lowmem : Read Agilent/Varian files from a directory using minimal
        amounts of memory.

    """
    # open the file
    f = open(filename, 'rb')

    # read the fileheader
    dic = fileheader2dic(get_fileheader(f))
    f.close()

    if dic["ntraces"] != 1:
        raise NotImplementedError

    i2tfunc = torder2i2t(torder)
    data = fid_nd(filename, i2tfunc, shape)
    return dic, data


def read_fid_ntraces(filename, shape=None, torder='flat', as_2d=False,
                     read_blockhead=False):
    """
    Read a Agilent/Varian binary (fid) file possibility having multiple
    traces per block.

    Parameters
    ----------
    filename : str
        Filename of Agilent/Varian binary file (fid) to read.
    shape : tuple of ints, optional
        Shape of the binary data. If not provided data is returned as a 2D
        array. Required if more than one trace per block (non-standard).
    torder : {'f', 'n', 'o'}
        Trace order. See :py:func:`read` for details.
    as_2d : bool, optional
        True to return the data as a 2D array, ignoring the shape and torder
        parameters.
    read_blockhead : bool, optional
        True to read the Agilent/Varian blockheaders(s) into the returned
        dictionary. False ignores them.

    Returns
    -------
    dic : dict
        Dictionary of Agilent/Varian binary file parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    See Also
    --------
    read_fid : Read a Agilent/Varian binary file with one trace per block.
    read_fid_lowmem : Read a Agilent/Varian binary file with one trace per
        block using minimal amounts of memory.


    """
    # open the file
    f = open(filename, 'rb')

    # read the fileheader
    dic = fileheader2dic(get_fileheader(f))

    # data parameters
    dt = find_dtype(dic)
    nblocks = dic["nblocks"]
    pts = dic["np"]
    nbheaders = dic["nbheaders"]
    ntraces = dic["ntraces"]

    # read the data
    if read_blockhead:
        bdic, data = get_nblocks_ntraces(f, nblocks, ntraces, pts,
                                         nbheaders, dt, read_blockhead)
        dic["blockheader"] = bdic
    else:
        data = get_nblocks_ntraces(f, nblocks, ntraces, pts, nbheaders, dt,
                                   read_blockhead)
    f.close()

    # uninterleave the real and imaginary data
    data = uninterleave_data(data)

    # if 2D array requested, return unshaped
    if as_2d:
        return dic, data

    # check for 1D
    if data.shape[0] == 1:
        return dic, np.squeeze(data)

    # try to reshape
    if shape is None:
        warn("unknown shape, returning unshaped data")
        return dic, data

    # reorder 3D/4D data
    if len(shape) >= 3:
        return dic, reorder_data(data, shape, torder)

    try:
        data = data.reshape(shape)
    except ValueError:
        warn(str(data.shape) + "cannot be shaped into" + str(shape))
        return dic, data

    return dic, data


def write_fid(filename, dic, data, torder='flat', repack=False, correct=True,
              overwrite=False):
    """
    Write a Agilent/Varian binary (fid) file.

    Parameters
    ----------
    filename : str
        Filename of Agilent/Varian binary file to write to.
    dic : dict
        Dictionary of Agilent/Varian file parameters.
    data : ndarray
        Array of NMR data.
    torder : {'f', 'r', 'o'}
        Trace ordering.  See :py:func:`read` for details.
    repack : bool, optional
        True to repack file and block headers. False leave as is.
    correct : bool, optional
        True (the default) will correct mis-sized np and nblocks values in dic.
        False will write out the incorrect values to the file header and
        record the full data set, the resulting file will not be readable by
        nmrglue.
    overwrite : bool, optional
        Set True to overwrite an existing file, False will raise a Warning if
        the file exists.

    See Also
    --------
    write_fid_lowmem : Write a Agilent/Varian binary file using minimal
        amounts of memory
    write : Write Agilent/Varian files to a directory.

    """
    data = np.array(data)

    # convert 1D data to 2D
    if data.ndim == 1:
        data = data.reshape((1, data.shape[0]))

    # reform 3D+ data
    if data.ndim >= 3:
        data = order_data(data, torder)

    # verify data and dic shapes
    if data.shape[1] != (dic['np'] // 2):
        warn("data and np size mismatch")
        if correct:
            dic['np'] = int(data.shape[1] * 2)
    if data.shape[0] != dic["nblocks"]:
        warn("data and block size mismatch")
        if correct:
            dic['nblocks'] = int(data.shape[0])

    # open file for writing
    f = fileiobase.open_towrite(filename, overwrite=overwrite)

    if repack:
        dic = repack_fileheader(dic)

    # write the fileheader to file
    put_fileheader(f, dic2fileheader(dic))

    # determine data type
    dt = find_dtype(dic)

    if "blockheader" in dic and len(dic["blockheader"]) == data.shape[0]:
        for i in range(data.shape[0]):
            if repack:
                bh = dic2blockheader(repack_blockheader(dic["blockheader"][0]))
            else:
                bh = dic2blockheader(dic["blockheader"][i])

            trace = np.array(interleave_data(data[i]), dtype=dt)
            put_block(f, trace, dic["nbheaders"], bh)

    else:   # create a generic blockheader
        bh = dic2blockheader(make_blockheader(dic, 1))
        for i in range(data.shape[0]):
            bh[2] = int(i + 1)
            trace = np.array(interleave_data(data[i]), dtype=dt)
            put_block(f, trace, dic["nbheaders"], bh)

    f.close()


def write_fid_lowmem(filename, dic, data, torder='f', repack=False,
                     overwrite=False):
    """
    Write a Agilent/Varian binary (fid) file using minimal amounts of memory.

    File is written trace by trace with each trace read from data before
    writing to reduce memory usage.

    Parameters
    ----------
    filename : str
        Filename of Agilent/Varian binary file to write to.
    dic : dict
        Dictionary of Agilent/Varian file parameters.
    data : array_like
        Array of NMR data.
    torder : {'f', 'r', 'o'}
        Trace ordering.  See :py:func:`read` for details.
    repack : bool, optional
        True to repack file and block headers. False leave as is.
    correct : bool, optional
        True (the default) will correct mis-sized np and nblocks values in dic.
        False will write out the incorrect values to the file header and
        record the full data set, the resulting file will not be readable by
        nmrglue.
    overwrite : bool, optional
        Set True to overwrite an existing file, False will raise a Warning if
        the file exists.

    See Also
    --------
    write_fid : Write a Agilent/Varian binary.
    write : Write Agilent/Varian files to a directory.

    """
    # convert 1D data to 2D
    if data.ndim == 1:
        data = data.reshape((1, data.shape[0]))

    t2i = torder2t2i(torder)

    # verify data and dic shapes
    if data.shape[-1] != (dic["np"] // 2):
        warn("data and np size mismatch")
        if correct:
            dic['np'] = int(data.shape[1] * 2)
    nblocks = int(reduce(operator.mul, data.shape[:-1]))
    if nblocks != dic["nblocks"]:
        warn("data and block size mismatch")
        if correct:
            dic['nblocks'] = nblocks

    # open file for writing
    f = fileiobase.open_towrite(filename, overwrite=overwrite)

    if repack:
        dic = repack_fileheader(dic)

    # write the fileheader to file
    put_fileheader(f, dic2fileheader(dic))

    # determine data type
    dt = find_dtype(dic)

    if "blockheader" in dic and len(dic["blockheader"]) == dic["nblocks"]:
        for ntrace in range(nblocks):
            if repack:
                bh = dic2blockheader(repack_blockheader(dic["blockheader"][0]))
            else:
                bh = dic2blockheader(dic["blockheader"][0])

            tup = t2i(data.shape[:-1], ntrace)
            trace = np.array(interleave_data(data[tup]), dtype=dt)
            put_block(f, trace, dic["nbheaders"], bh)

    else:   # create a generic blockheader
        bh = dic2blockheader(make_blockheader(dic, 1))
        for ntrace in range(nblocks):
            bh[2] = int(ntrace + 1)
            tup = t2i(data.shape[:-1], ntrace)
            trace = np.array(interleave_data(data[tup]), dtype=dt)
            put_block(f, trace, dic["nbheaders"], bh)
    f.close()


#####################
# get/put functions #
#####################


def get_nblocks(f, nblocks, pts, nbheaders, dt, read_blockhead):
    """
    Read multiple blocks from a Agilent/Varian binary file.

    Parameters
    ----------
    f : file object
        File object of open Agilent/Varian binary file.
    nblock : int
        Number of blocks to read.
    pts : int
        Number of points per trace.
    nbheaders : int
        Number of block header in each block.
    dt : dtype
        Data type of data in binary file (real type).
    read_blockheaders : bool
        True to read the blockheader(s) and return them in the returned
        dictionary.  False ignores blockheader, reading over them.

    Returns
    -------
    dic : dict, optional
        Dictionary of blockheaders, only returned if read_blockheaders is True.
    data : 2D ndarray of size (nblocks, pts)
        Array of read blocks.

    """
    # create an empty array to hold data
    data = np.empty((nblocks, pts), dtype=dt)
    if read_blockhead:
        bdic = [0] * nblocks

    # read the data
    for i in range(nblocks):
        if read_blockhead:
            bdic[i], data[i] = get_block(f, pts, nbheaders, dt, read_blockhead)
        else:
            data[i] = get_block(f, pts, nbheaders, dt, read_blockhead)

    # return the data
    if read_blockhead:
        return bdic, data
    else:
        return data


def get_block(f, pts, nbheaders, dt, read_blockhead=False):
    """
    Read a single block from Agilent/Varian binary file.

    Parameters
    ----------
    f : file object
        File object of open Agilent/Varian binary file.
    pts : int
        Number of points per trace.
    nbheaders : int
        Number of block header in each block.
    dt : dtype
        Data type of data in binary file (real type).
    read_blockheaders : bool, optional.
        True to read the blockheader(s) and return them in the returned
        dictionary.  False (default) ignores blockheader, reading over them.

    Returns
    -------
    dic : dict, optional
        Dictionary of blockheaders, only returned if read_blockheaders is True.
    data : ndarray of size
        Array of read block.

    """
    if read_blockhead is False:  # do not return blockheaders
        for i in range(nbheaders):
            skip_blockheader(f)
        trace = get_trace(f, pts, dt)
        return trace

    else:   # read the block headers
        dic = dict()
        if nbheaders >= 1:
            dic.update(blockheader2dic(get_blockheader(f)))
        if nbheaders >= 2:
            dic["hyperhead"] = hyperheader2dic(get_hyperheader(f))
        if nbheaders >= 3:
            for i in range(2, nbheaders):
                skip_blockheader(f)
        # read the data
        trace = get_trace(f, pts, dt)

        return dic, trace


def get_nblocks_ntraces(f, nblocks, ntraces, pts, nbheaders, dt,
                        read_blockhead):
    """
    Read multiple blocks from a Agilent/Varian binary file which may have
    multiple traces per block.

    Parameters
    ----------
    f : file object
        File object of open Agilent/Varian binary file.
    nblock : int
        Number of blocks to read.
    ntrace : int
        Number of traces per block
    pts : int
        Number of points per trace.
    nbheaders : int
        Number of block header in each block.
    dt : dtype
        Data type of data in binary file (real type).
    read_blockheaders : bool
        True to read the blockheader(s) and return them in the returned
        dictionary.  False ignores blockheader, reading over them.

    Returns
    -------
    dic : dict, optional
        Dictionary of blockheaders, only returned if read_blockheaders is True.
    data : 2D ndarray of size (nblocks * ntraces, pts)
        Array of read blocks.

    """
    # create an empty array to hold data
    data = np.empty((nblocks * ntraces, pts), dtype=dt)
    if read_blockhead:
        bdic = [0] * nblocks

    # read the data
    for i in range(nblocks):
        if read_blockhead:
            bdic[i], bdata = get_block_ntraces(f, ntraces, pts, nbheaders, dt,
                                               True)
            data[i * ntraces:(i + 1) * ntraces] = bdata
        else:
            bdata = get_block_ntraces(f, ntraces, pts, nbheaders, dt, False)
            data[i * ntraces:(i + 1) * ntraces] = bdata
    if read_blockhead:
        return bdic, data
    else:
        return data


def get_block_ntraces(f, ntraces, pts, nbheaders, dt, read_blockhead=False):
    """
    Read a single block from Agilent/Varian binary file which may have
    multiple traces per block.

    Parameters
    ----------
    f : file object
        File object of open Agilent/Varian binary file.
    ntrace : int
        Number of traces per block
    pts : int
        Number of points per trace.
    nbheaders : int
        Number of block header in each block.
    dt : dtype
        Data type of data in binary file (real type).
    read_blockheaders : bool
        True to read the blockheader(s) and return them in the returned
        dictionary.  False ignores blockheader, reading over them.

    Returns
    -------
    dic : dict, optional
        Dictionary of blockheaders, only returned if read_blockheaders is True.
    data : 2D ndarray of shape (ntraces, pts)
        Array of read blocks.

    """

    if read_blockhead is False:  # do not return blockheaders
        for i in range(nbheaders):
            skip_blockheader(f)
        trace = get_trace(f, pts * ntraces, dt)
        return trace.reshape(ntraces, pts)
    else:   # read the blockheaders
        dic = dict()
        # read the headers
        if nbheaders >= 1:
            dic.update(blockheader2dic(get_blockheader(f)))
        if nbheaders >= 2:
            dic["hyperhead"] = hyperheader2dic(get_hyperheader(f))
        if nbheaders >= 3:
            for i in range(2, nbheaders):
                skip_blockheader(f)
        # read the data
        trace = get_trace(file, pts * ntraces, dt)
        return dic, trace.reshape(ntraces, pts)


def get_trace(f, pts, dt):
    """
    Read trace of pts points of dtype dt from Agilent/Varian binary file

    Endianness should be handled by dt.
    """
    bsize = pts * dt.itemsize  # number of bytes in trace
    return np.frombuffer(f.read(bsize), dt)


def get_fileheader(f):
    """
    Unpack file header parameters into a list.

    Reads the 32-byte file header from file and unpacks into a list.  Endianness
    is corrected as needed.

    Returned list contents:

    =   =========   ======================================
    N   Variable    Description
    =   =========   ======================================
    0   nblocks     data blocks in file
    1   ntraces     traces per block
    2   np          elements per trace
    3   ebytes      bytes per element
    4   tbytes      bytes per trace
    5   bbytes      bytes per block
    6   vers_id     software version, file_id status bits
    7   status      status of whole file
    8   nbheaders   number of block headers per block (1)
    =   =========   ======================================

    """
    # header is packed big-endian as 6 longs, 2 shorts, 1 long
    return struct.unpack('>6lhhl', f.read(32))


def get_blockheader(f):
    """
    Unpack block header parameters into a list.

    Reads the 28-byte block header from f and unpacks into a list.  Endianness
    is corrected as needed.

    Returned list contents:

    =   ========    ========================
    N   Variable    Description
    =   ========    ========================
    0   scale       scaling factor
    1   status      status of data in block
    2   index       block index
    3   mode        block mode
    4   ctcount     ct value of FID
    5   lpval       left phase
    6   rpval       right phase
    7   lvl         level drift correction
    8   tlt         tilt drift correction
    =   ========    ========================

    """
    # block header is packed big-endian as 4 shorts, 1 long, 4 floats
    return struct.unpack('>4hl4f', f.read(28))


def skip_blockheader(f):
    """
    Read a block header but do not unpack.

    This is a replacement for get_blockheader.  It skips f ahead 28 bytes.
    """
    f.read(28)


def get_hyperheader(file):
    """
    Unpack hypercomplex header parameters to a list.

    Reads the 28-bytes block header from file and unpacks into a list. Endianness
    is corrected as needed.

    Returned list contents:

    =   ========    ================
    N   Variable    Description
    =   ========    ================
    0   s_spare1    Not Used
    1   status      status of block
    2   s_spare2    Not Used
    3   s_spart3    Not Used
    4   l_spare1    Not Used
    5   lpval1      2D left phase
    6   rpval1      2D right phase
    7   f_spare1    Not Used
    8   f_spare2    Not Used
    =   ========    ================

    """
    # hypercomplex header is packed big-endian as 4 shorts, 1 long, 4 floats
    return struct.unpack('>4hl4f', file.read(28))


def put_block(f, trace, nbheaders, bh, hh=False):
    """
    Put blockheader(s) and the trace to file.

    Parameters
    ----------
    f : file object
        Open file object to write to.
    trace : ndarray
        Trace to write to current block.
    nbheaders : int
        Number of block headers.
    bh : list
        Blockheader list.
    hh : list, optional
        Hyperheader list. Required when nbheaders >= 2.

    Notes
    -----
    When nbheaders > 2, additional headers are written as all zeros.

    """
    # write the block headers
    if nbheaders >= 1:
        put_blockheader(f, bh)

    if nbheaders >= 2:
        if hh is False:
            raise Exception("Hyperheader required")
        put_hyperheader(f, hh)

    # write any additional blockheaders as 0s if needed
    for i in range(nbheaders - 2):
        put_blockheader(f, [0] * 9)

    # write the trace
    put_trace(f, trace)


def put_trace(f, trace):
    """
    Write a trace to file f.
    """
    f.write(trace.tobytes())


def put_fileheader(f, fh):
    """
    Write a fileheader list to file (32-bytes written).

    Parameters
    ----------
    f : file object
        Open file object to write to.
    fh : list with 9 elements
        Fileheader list.

    """
    f.write(struct.pack('>6lhhl', *fh))


def put_blockheader(f, bh):
    """
    Write a blockheader list to file (28-bytes written)

    Parameters
    ----------
    f : file object
        Open file object to write to.
    bh : list with 9 elements
        Blockheaders list.

    """
    f.write(struct.pack('>4hl4f', *bh))


def put_hyperheader(f, hh):
    """
    Write hyperheader list to file (28-bytes written)

    Parameters
    ----------
    f : file object
        Open file object to write to.
    hh : list with 9 elements
        Hyperheader list.

    """
    f.write(struct.pack('>4hl4f', *hh))


#########################
# dictionary conversion #
#########################


def hyperheader2dic(head):
    """
    Convert a hypercomplex block header into a Python dictionary.
    """
    dic = dict()
    dic["s_spare1"] = head[0]
    dic["status"] = head[1]
    dic["s_spare2"] = head[2]
    dic["s_spare3"] = head[3]
    dic["l_spare1"] = head[4]
    dic["lpval1"] = head[5]
    dic["rpval1"] = head[6]
    dic["f_spare1"] = head[7]
    dic["f_spare2"] = head[8]

    # unpack the status bits
    dic["UHYPERCOMPLEX"] = (dic["status"] & 0x2) // 0x2

    return dic


def repack_hyperheader(dic):
    """
    Repack a hyperheader dictionary bit flag parameters into status.
    """
    dic["status"] = dic["UHYPERCOMPLEX"] * 0x2
    return dic


def dic2hyperheader(dic):
    """
    Convert a Python dictionary into a hypercomplex block header list.

    Does not repack status from bit flags.
    """
    head = [0] * 9
    head[0] = dic["s_spare1"]
    head[1] = dic["status"]
    head[2] = dic["s_spare2"]
    head[3] = dic["s_spare3"]
    head[4] = dic["l_spare1"]
    head[5] = dic["lpval1"]
    head[6] = dic["rpval1"]
    head[7] = dic["f_spare1"]
    head[8] = dic["f_spare2"]
    return head


def make_blockheader(filedic=False, index=1):
    """
    Make a generic blockheader dictionary with a given block index.

    filedic can be provided for status flags, if not provided creates
    header for float32 data

    """
    dic = dict()

    if filedic is False:
        # common status flags
        dic["S_DATA"] = filedic["S_DATA"]
        dic["S_SPEC"] = filedic["S_SPEC"]
        dic["S_32"] = filedic["S_32"]
        dic["S_FLOAT"] = filedic["S_FLOAT"]
        dic["S_COMPLEX"] = filedic["S_COMPLEX"]
        dic["S_HYPERCOMPLEX"] = filedic["S_HYPERCOMPLEX"]
    else:
        # common status flags
        dic["S_DATA"] = 1
        dic["S_SPEC"] = 0
        dic["S_32"] = 0
        dic["S_FLOAT"] = 1
        dic["S_COMPLEX"] = 0
        dic["S_HYPERCOMPLEX"] = 0

    # blockheader specific status flags
    dic["MORE_BLOCKS"] = 1
    dic["NP_CMPLX"] = 0
    dic["NF_CMPLX"] = 0
    dic["NI_CMPLX"] = 0
    dic["NI2_CMPLX"] = 0

    # mode flags
    dic["NP_PHMODE"] = 0
    dic["NP_AVMODE"] = 0
    dic["NP_PWRMODE"] = 0
    dic["NF_PHMODE"] = 0
    dic["NF_AVMODE"] = 0
    dic["NF_PWRMODE"] = 0
    dic["NI_PHMODE"] = 0
    dic["NI_AVMODE"] = 0
    dic["NI_PWRMODE"] = 0
    dic["NI2_PHMODE"] = 0
    dic["NI2_AVMODE"] = 0
    dic["NI2_PWRMODE"] = 0

    # main parameters
    dic["scale"] = 0
    dic["status"] = 0
    dic["index"] = index
    dic["mode"] = 0
    dic["ctcount"] = 1
    dic["lpval"] = 0.0
    dic["rpval"] = 0.0
    dic["lvl"] = 0.0
    dic["tlt"] = 0.0

    repack_blockheader(dic)  # repack the header
    return dic


def blockheader2dic(head):
    """
    Convert a block header list into a Python dictionary.
    """
    dic = dict()

    dic["scale"] = head[0]
    dic["status"] = head[1]
    dic["index"] = head[2]
    dic["mode"] = head[3]
    dic["ctcount"] = head[4]
    dic["lpval"] = head[5]
    dic["rpval"] = head[6]
    dic["lvl"] = head[7]
    dic["tlt"] = head[8]

    # unpack the status parameters
    dic["S_DATA"] = (dic["status"] & 0x1) // 0x1
    dic["S_SPEC"] = (dic["status"] & 0x2) // 0x2
    dic["S_32"] = (dic["status"] & 0x4) // 0x4
    dic["S_FLOAT"] = (dic["status"] & 0x8) // 0x8
    dic["S_COMPLEX"] = (dic["status"] & 0x10) // 0x10
    dic["S_HYPERCOMPLEX"] = (dic["status"] & 0x20) // 0x20

    dic["MORE_BLOCKS"] = (dic["status"] & 0x80) // 0x80
    dic["NP_CMPLX"] = (dic["status"] & 0x100) // 0x100
    dic["NF_CMPLX"] = (dic["status"] & 0x200) // 0x200
    dic["NI_CMPLX"] = (dic["status"] & 0x400) // 0x400
    dic["NI2_CMPLX"] = (dic["status"] & 0x800) // 0x800

    # unpack the mode parameter
    dic["NP_PHMODE"] = (dic["mode"] & 0x1) // 0x1
    dic["NP_AVMODE"] = (dic["mode"] & 0x2) // 0x2
    dic["NP_PWRMODE"] = (dic["mode"] & 0x4) // 0x4

    dic["NF_PHMODE"] = (dic["mode"] & 0x10) // 0x10
    dic["NF_AVMODE"] = (dic["mode"] & 0x20) // 0x20
    dic["NF_PWRMODE"] = (dic["mode"] & 0x40) // 0x40

    dic["NI_PHMODE"] = (dic["mode"] & 0x100) // 0x100
    dic["NI_AVMODE"] = (dic["mode"] & 0x200) // 0x200
    dic["NI_PWRMODE"] = (dic["mode"] & 0x400) // 0x400

    dic["NI2_PHMODE"] = (dic["mode"] & 0x1000) // 0x1000
    dic["NI2_AVMODE"] = (dic["mode"] & 0x2000) // 0x2000
    dic["NI2_PWRMODE"] = (dic["mode"] & 0x4000) // 0x4000

    return dic


def repack_blockheader(dic):
    """
    Repack blockheader dic bit flag parameters into status and mode.
    """

    dic["status"] = (dic["S_DATA"] * 0x1 + dic["S_SPEC"] * 0x2 +
                     dic["S_32"] * 0x4 + dic["S_FLOAT"] * 0x8 +
                     dic["S_COMPLEX"] * 0x10 + dic["S_HYPERCOMPLEX"] * 0x20 +
                     dic["MORE_BLOCKS"] * 0x80 + dic["NP_CMPLX"] * 0x100 +
                     dic["NF_CMPLX"] * 0x200 + dic["NI_CMPLX"] * 0x400 +
                     dic["NI2_CMPLX"] * 0x800)

    dic["mode"] = (dic["NP_PHMODE"] * 0x1 + dic["NP_AVMODE"] * 0x2 +
                   dic["NP_PWRMODE"] * 0x4 + dic["NF_PHMODE"] * 0x10 +
                   dic["NF_AVMODE"] * 0x20 + dic["NF_PWRMODE"] * 0x40 +
                   dic["NI_PHMODE"] * 0x100 + dic["NI_AVMODE"] * 0x200 +
                   dic["NI_PWRMODE"] * 0x400 + dic["NI2_PHMODE"] * 0x1000 +
                   dic["NI2_AVMODE"] * 0x2000 + dic["NI2_PWRMODE"] * 0x4000)
    return dic


def dic2blockheader(dic):
    """
    Convert a python dictionary into block header list.

    Does not repack status and mode from bit flags.
    """
    head = [0] * 9
    head[0] = dic["scale"]
    head[1] = dic["status"]
    head[2] = dic["index"]
    head[3] = dic["mode"]
    head[4] = dic["ctcount"]
    head[5] = dic["lpval"]
    head[6] = dic["rpval"]
    head[7] = dic["lvl"]
    head[8] = dic["tlt"]
    return head


def fileheader2dic(head):
    """
    Convert fileheader list into a Python dictionary
    """
    dic = dict()
    dic["nblocks"] = head[0]
    dic["ntraces"] = head[1]
    dic["np"] = head[2]
    dic["ebytes"] = head[3]
    dic["tbytes"] = head[4]
    dic["bbytes"] = head[5]
    dic["vers_id"] = head[6]
    dic["status"] = head[7]
    dic["nbheaders"] = head[8]

    # unpack the status parameter
    dic["S_DATA"] = (dic["status"] & 0x1) // 0x1
    dic["S_SPEC"] = (dic["status"] & 0x2) // 0x2
    dic["S_32"] = (dic["status"] & 0x4) // 0x4
    dic["S_FLOAT"] = (dic["status"] & 0x8) // 0x8
    dic["S_COMPLEX"] = (dic["status"] & 0x10) // 0x10
    dic["S_HYPERCOMPLEX"] = (dic["status"] & 0x20) // 0x20
    dic["S_ACQPAR"] = (dic["status"] & 0x80) // 0x80
    dic["S_SECND"] = (dic["status"] & 0x100) // 0x100
    dic["S_TRANSF"] = (dic["status"] & 0x200) // 0x200
    dic["S_NP"] = (dic["status"] & 0x800) // 0x800
    dic["S_NF"] = (dic["status"] & 0x1000) // 0x1000
    dic["S_NI"] = (dic["status"] & 0x2000) // 0x2000
    dic["S_NI2"] = (dic["status"] & 0x4000) // 0x4000

    return dic


def repack_fileheader(dic):
    """
    Repack blockheader dic bit flag parameters into status and mode.
    """
    dic["status"] = (dic["S_DATA"] * 0x1 + dic["S_SPEC"] * 0x2 +
                     dic["S_32"] * 0x4 + dic["S_FLOAT"] * 0x8 +
                     dic["S_COMPLEX"] * 0x10 + dic["S_HYPERCOMPLEX"] * 0x20 +
                     dic["S_ACQPAR"] * 0x80 + dic["S_SECND"] * 0x100 +
                     dic["S_TRANSF"] * 0x200 + dic["S_NP"] * 0x800 +
                     dic["S_NF"] * 0x1000 + dic["S_NI"] * 0x2000 +
                     dic["S_NI2"] * 0x4000)
    return dic


def dic2fileheader(dic):
    """
    Convert a Python dictionary into a fileheader list

    Does not repack status from bit flags
    """

    head = [0] * 9
    head[0] = dic["nblocks"]
    head[1] = dic["ntraces"]
    head[2] = dic["np"]
    head[3] = dic["ebytes"]
    head[4] = dic["tbytes"]
    head[5] = dic["bbytes"]
    head[6] = dic["vers_id"]
    head[7] = dic["status"]
    head[8] = dic["nbheaders"]

    return head

##################
# misc functions #
##################


def find_shape(pdic):
    """
    Determine the shape of a Agilent/Varian file from the procpar dictionary
    """
    # Varian files are typically either from imaging experiments or
    # from NMR experiments.  For both the direct dimension (R+I) is stored in
    # the np parameters in the procpar file.  In addition an array may be
    # present.

    # Imaging experiments have a 'seqcon' parameter present in the procpar and
    # have indirect dimension are set by nv, nv2, and nv3.

    # NMR experiments has indirect dimension shapes (R|I) are stored in
    # ni,ni2,ni3 with the phase arrays stored as phase,phase2,phase3

    # Thanks to the VeSPA project (http://scion.duhs.duke.edu/vespa/) for
    # information on imaging experiments

    try:
        shape = [int(pdic["np"]["values"][0]) // 2]
    except:
        # when shape finding fails issue warning and return None
        warn("shape cannot be determined.")
        return None

    # check for an array, we will deal only with the inner most array.
    # When multiple parameters are arrayed the ordering of indirect dimension
    # must also be known, which typically is not provided, so we will leave
    # this edge case to the user.
    if "array" in pdic:
        array_name = pdic['array']['values'][0]
        array_name = array_name.split(',')[-1]   # keep only innermost array
        if array_name.startswith('phase') is False and array_name in pdic:
            shape.insert(0, len(pdic[array_name]["values"]))

    # imaging files have a seqcon parameter
    if 'seqcon' in pdic:
        if "nv" in pdic:
            s = max(int(pdic["nv"]["values"][0]), 1)
            if s > 1:
                shape.insert(0, s)
        if "nv2" in pdic:
            s = max(int(pdic["nv2"]["values"][0]), 1)
            if s > 1:
                shape.insert(0, s)

        if "nv3" in pdic:
            s = max(int(pdic["nv3"]["values"][0]), 1)
            if s > 1:
                shape.insert(0, s)
        return tuple(shape)

    # assume we have NMR data
    if "ni" in pdic:
        multi = 2       # assume R+I in cases where no phase parameter.
        if "phase" in pdic:
            multi = len(pdic["phase"]["values"])
        s = max(int(pdic["ni"]["values"][0]), 1)
        shape.insert(0, s * multi)

    if "ni2" in pdic:
        multi = 2
        if "phase2" in pdic:
            multi = len(pdic["phase2"]["values"])
        s = max(int(pdic["ni2"]["values"][0]), 1)
        shape.insert(0, s * multi)

    if "ni3" in pdic:
        multi = 2
        if "phase3" in pdic:
            multi = len(pdic["phase"]["values"])
        s = max(int(pdic["ni3"]["values"][0]), 1)
        shape.insert(0, s * multi)

    return tuple(shape)


def find_cdtype(dic):
    """
    Find the complex dtype from a Agilent/Varian dictionary
    """
    if dic["S_FLOAT"] == 1:
        return np.dtype("complex64")
    else:
        if dic["S_32"] == 1:
            return np.dtype("complex128")
        else:
            return np.dtype("complex64")


def find_dtype(dic):
    """
    Find the real dtype from a dictionary
    """
    if dic["S_FLOAT"] == 1:
        return np.dtype('>f4')  # float32
    else:
        if dic["S_32"] == 1:
            return np.dtype('>i4')  # int32
        else:
            return np.dtype('>i2')  # int16


def uninterleave_data(data):
    """
    Unpack interleaved real, imag data

    ==========  ============
    data dtype  Return dtype
    ==========  ============
    int16       'complex64'
    float32     'complex64'
    int32       'complex128'
    ==========  ============

    """
    # determine the output dtype
    rdt = data.dtype.name

    if rdt in ('int16', "float32"):
        cdt = "complex64"
    elif rdt == 'int32':
        cdt = "complex128"
    else:
        raise ValueError("unknown dtype")

    return data[..., ::2] + np.array(data[..., 1::2] * 1.j, dtype=cdt)


def interleave_data(data_in):
    """
    Interleave real, imag data

    Does not check if resulting dtype is a valid Agilent/Varian dtype
    """
    size = list(data_in.shape)
    size[-1] = size[-1] * 2
    data_out = np.empty(size, dtype=data_in.real.dtype)
    data_out[..., ::2] = data_in.real
    data_out[..., 1::2] = data_in.imag
    return data_out

###########################
# procpar reading/writing #
###########################

# procpar functions


def read_procpar(filename):
    """
    Read a procpar file returning a dictionary of procpar parameters.
    """
    with open(filename, 'rb') as f:
        dic = dict()
        length = os.stat(filename).st_size

        # test to see if end of file
        while f.tell() != length:
            p = get_parameter(f)
            dic[p["name"]] = p
        return dic


def get_parameter(f):
    """
    Reads a procpar parameter from a file object.

    Returns a dictionary with the attributes of the parameter.

    """
    dic = dict()

    # read and decode the first line
    line = f.readline().decode().split()

    dic["name"] = line[0]
    dic["subtype"] = line[1]
    dic["basictype"] = line[2]
    dic["maxvalue"] = line[3]
    dic["minvalue"] = line[4]
    dic["stepsize"] = line[5]
    dic["Ggroup"] = line[6]
    dic["Dgroup"] = line[7]
    dic["protection"] = line[8]
    dic["active"] = line[9]
    dic["intptr"] = line[10]

    # read in the values of the parameter
    line = f.readline().decode()
    num = int(line.split()[0])
    values = []

    if dic["basictype"] == "1":     # real values, only one line
        values = line.split()[1:]
    elif dic["basictype"] == "2":   # strings, may have multiple lines
        values.append(line.split("\"")[1])  # split on "s
        for i in range(num - 1):
            values.append(f.readline().decode().split("\"")[1])

    dic["values"] = values
    line = f.readline().decode()

    # read and decode the enumerables
    dic["enumerable"] = line.split()[0]

    if dic["enumerable"] != "0":
        if dic["basictype"] == "1":     # reals
            dic["enumerables"] = line.split()[1:]
        elif dic["basictype"] == "2":   # strings
            dic["enumerables"] = line.split("\"")[1::2]

    return dic


def write_procpar(filename, dic, overwrite=False):
    """
    Write a Agilent/Varian procpar file from a dictionary
    """
    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite=overwrite, mode='w')

    for key in dic.keys():  # loop over the parameters

        d = dic[key]
        # print(out the first line)
        print(d["name"], d["subtype"], d["basictype"],
              d["maxvalue"], d["minvalue"], d["stepsize"],
              d["Ggroup"], d["Dgroup"], d["protection"],
              d["active"], d["intptr"], file=f)

        # print out the next line(s) (and more if strings)
        if d["basictype"] == "1":   # real values, one line
            print(len(d["values"]), end=' ', file=f)  # don't end the line
            for value in d["values"]:
                print(value, end=' ', file=f)    # still not yet
            print("", file=f)   # now end the line

        elif d["basictype"] == "2":     # strings, may have multiple lines
            print(len(d["values"]), end=' ', file=f)  # don't end the line
            for value in d["values"]:
                # now end the line (for each string)
                print('"' + value + '"', file=f)

        # print out the last line
        print(d["enumerable"], end=' ', file=f)

        if d["enumerable"] != "0":
            for e in d["enumerables"]:
                if d["basictype"] == "1":  # reals
                    print(e, end=' ', file=f)
                elif d["basictype"] == "2":  # strings
                    print('"' + e + '"', end=' ', file=f)
        print("", file=f)   # end the enumerable line

    f.close()


subtypes = ["undefined", "real", "string", "delay", "flag", "frequency",
            "pulse", "integer"]

basictypes = ["undefined", "real", "string"]

############################################
# low memory numpy.ndarray emulating class #
############################################


class fid_nd(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of Agilent/Varian fid files which must have one trace per block.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of Agilent/Varian binary file.
    i2t_func : function
        Python function which maps an index to a trace.
    fshape : tuple of ints, optional
        Shape of data in file, if None will be assumed to be 2D data.
    order : tuple
        Ordering of axes compared to file. None will results in (0, 1, 2, ...)
        ordering.

    """
    def __init__(self, filename, i2t_func, fshape=None, order=None):
        """
        Create and set up object.
        """
        # read the file dictionary
        f = open(filename, 'rb')
        dic = fileheader2dic(get_fileheader(f))
        f.close()

        # check fshape
        if fshape is None:
            # by default open as 2D with nblocks,np/2 fshape
            fshape = (dic['nblocks'], int(dic["np"] / 2))
        else:
            # check that last dimension has np/2 points
            if fshape[-1] != int(dic["np"] / 2):
                s = "last dimension should have size %i" % (int(dic["np"] / 2))
                raise ValueError(s)
            # product of all but last dim should be number of blocks
            if reduce(operator.mul, fshape[:-1]) != dic['nblocks']:
                s = "number of traces in file does not match fshape"
                raise ValueError(s)

        # check order
        if order is None:
            order = range(len(fshape))

        # finalize
        self.fdtype = find_dtype(dic)
        self.pts = dic["np"]
        self.nbh = dic["nbheaders"]
        self.bbytes = dic["bbytes"]
        self.dtype = find_cdtype(dic)
        self.filename = filename
        self.i2t = i2t_func
        self.fshape = fshape
        self.order = order
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy.
        """
        n = fid_nd(self.filename, self.i2t, self.fshape, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values

        slices is a well formatted tuple of slices
        """
        # separate the last slice from the first slices
        lslice = slices[-1]
        fslice = slices[:-1]

        # and the same for fshape
        lfshape = self.fshape[-1]
        ffshape = self.fshape[:-1]

        # find the output size and make a in/out nd iterator
        osize, nd_iter = fileiobase.size_and_ndtofrom_iter(ffshape, fslice)
        osize.append(len(range(lfshape)[lslice]))

        # create an empty array to store the selected slices
        out = np.empty(tuple(osize), dtype=self.dtype)

        with open(self.filename, 'rb') as f:

            # read in the data trace by trace
            for out_index, in_index in nd_iter:

                # determine the trace number from the index
                ntrace = self.i2t(ffshape, in_index)

                # seek to the correct place in the file
                f.seek(ntrace * self.bbytes + 32)

                # retrieve trace and save to output
                trace = get_block(f, self.pts, self.nbh, self.fdtype, False)
                trace = uninterleave_data(trace)
                out[out_index] = trace[lslice]

        return out
