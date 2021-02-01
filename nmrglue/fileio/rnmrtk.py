"""
Fuctions for reading and writing Rowland NMR Toolkit (RNMRTK) files
"""

from __future__ import division

__developer_info__ = """
Information of the Rowland NMR Toolkit file format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

import numpy as np
from warnings import warn

from . import fileiobase

###################
# unit conversion #
###################


def make_uc(dic, data, dim=-1):
    """
    Creat a unit conversion object

    Parameters
    ----------
    dic : dict
        Dictionary of RNMRTK parameters.
    data : ndarray
        Array of NMR data.
    dim : int, optional
        Demension number to create unit conversion object for. Default is for
        the last dimension.

    Returns
    -------
    uc : unit conversion object.
        Unit conversion object for given dimension.

    """
    if dim < 0:     # negative dimensions
        dim = data.ndim + dim
    size = data.shape[dim]  # R|I

    ddim = find_dic_dim(dic, dim)
    cplx = {'R': False, 'C': True}[dic['nptype'][ddim]]
    sw = dic['sw'][ddim]   # Hz
    obs = dic['sf'][ddim]   # MHz
    car = dic['ppm'][ddim] * dic[obs]   # Hz
    return fileiobase.unit_conversion(size, cplx, sw, obs, car)


#################
# data creation #
#################


def create_data(data):
    """
    Create a RNMRTK data array (recast into float32 or complex64)
    """
    if np.iscomplexobj(data):
        return np.array(data, dtype="complex64")
    else:
        return np.array(data, dtype="float32")

########################
# universal dictionary #
########################


def guess_udic(dic, data):
    """
    Guess parameters of a universal dictionary from a dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of RNMRTK parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.

    """
    # create an empty universal dictionary
    ndim = dic['ndim']
    udic = fileiobase.create_blank_udic(ndim)

    # fill in parameters from RNMRTK dictionary for each dimension
    for iudim in range(ndim):

        # find the corresponding dimension in the RNMRTK parameter dictionary
        idim = find_dic_dim(dic, iudim)

        udic[iudim]['encoding'] = dic['quad'][idim].lower()
        udic[iudim]['car'] = dic['ppm'][idim] * dic['sf'][idim]
        udic[iudim]['obs'] = dic['sf'][idim]
        udic[iudim]['sw'] = dic['sw'][idim]
        udic[iudim]['size'] = dic['npts'][idim]

        # set quadrature and correct size
        if dic['nptype'][idim] == 'C':
            if iudim != ndim - 1:
                udic[iudim]['size'] *= 2    # don't double size of last dim
            udic[iudim]['complex'] = True
        else:
            udic[iudim]['complex'] = False

        # set label to T1 or F1, etc
        udic[iudim]['label'] = dic['dom'][idim] + str(idim + 1)

        # set time or frequency domain
        if dic['dom'][idim].upper() == 'T':
            udic[iudim]['freq'] = False
            udic[iudim]['time'] = True
        else:
            udic[iudim]['freq'] = True
            udic[iudim]['time'] = False

    return udic


def create_dic(udic, dim_order=None):
    """
    Create a RNMRTK dictionary from a universal dictionary.

    Parameters
    ----------
    udic : dict
        Universal dictionary of spectral parameters.
    dim_order : list, optional
        List mapping axis numbers in the universal dictionary to the order
        in which they will appear in the RNMRTK dictionary.  If None, the
        default [0, 1, 2, ...] is used.

    Returns
    --------
    dic : dict
        Dictionary of RNMRTK parameters.

    """

    # create the RNMRTK dictionary and fill with some default values
    dic = {}
    dic['comment'] = ''
    dic['ndim'] = ndim = int(udic['ndim'])
    dic['format'] = np.dtype('float32').str

    if dim_order is None:
        dim_order = range(ndim)  # default to 0, 1, 2, ...

    # set various parameters from the universal dictionary
    dic['dom'] = [['F', 'T'][udic[i]['time']] for i in dim_order]
    dic['nptype'] = [['R', 'C'][udic[i]['complex']] for i in dim_order]
    dic['ppm'] = [udic[i]['car'] / udic[i]['obs'] for i in dim_order]
    dic['sf'] = [udic[i]['obs'] for i in dim_order]
    dic['sw'] = [udic[i]['sw'] for i in dim_order]
    dic['npts'] = [udic[i]['size'] for i in dim_order]
    dic['quad'] = [udic[i]['encoding'].lower() for i in dim_order]

    # xfirst and xstep are freq domains values, correct later for time domain
    dic['xfirst'] = [-0.5 * i for i in dic['sw']]
    dic['xstep'] = [udic[i]['sw'] / udic[i]['size'] for i in dim_order]

    # these we guess on, they may be incorrect
    dic['cphase'] = [0.0] * ndim
    dic['lphase'] = [0.0] * ndim
    dic['nacq'] = [udic[i]['size'] for i in dim_order]

    # make small corrections as needed
    rnmrtk_quads = ['states', 'states-tppi', 'tppi', 'tppi-redfield']
    for i in range(ndim):

        # fix quadrature if not a valid RNMRTK quadrature
        if dic['quad'][i] not in rnmrtk_quads:
            dic['quad'][i] = 'states'

        # fix parameters if time domain data
        if dic['dom'][i] == 'T':    # time domain data
            dic['xfirst'][i] = 0.0
            if dic['quad'][i] in ['states']:
                # states time domain data
                dic['xstep'][i] = 1. / dic['sw'][i]
            else:
                # tppi time domain data
                dic['xstep'][i] = 0.5 / dic['sw'][i]

        # half the number of points if dimension is complex
        if dic['nptype'][i] == 'C':
            dic['npts'][i] //= 2

    # determine and set layout
    size = [udic[i]['size'] for i in range(ndim)]
    domains = [dic['dom'][i] + str(i + 1) for i in range(ndim)]

    # correct size of last dimension if complex
    if dic['nptype'][dim_order[-1]] == 'C':
        dic['npts'][dim_order[-1]] *= 2
        size[-1] *= 2
        if dic['dom'][dim_order[-1]] == 'F':
            dic['xstep'][dim_order[-1]] /= 2.

    dic['layout'] = (size, domains)

    return dic

#######################
# Reading and Writing #
#######################


def read(filename, par_file=None):
    """
    Read RNMRTK files.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file to read (.sec).
    par_file : str or None, optional
        Filename of RNMRTK parameter file. If None (default) a the last four
        characters of `file` are changed to .par.

    Returns
    -------
    dic : dic
        Dictionary of RNMRTK parameters.
    data : ndarray
        Array of NMR data.

    Notes
    -----
    The dictionary parameters are ordered opposite the data layout, that is to
    say the the FIRST parameter in each list corresponds to the LAST axis in
    the data array.

    See Also
    --------
    read_lowmem : Read RNMRTK files with minimal memory usage.
    write : Write RNMRTK files.

    """
    # determine par_file name if not given
    if par_file is None:
        par_file = filename[:-4] + ".par"
    dic = read_par(par_file)

    # determine sec file parameters from parameter dictionary
    dtype = dic["format"]
    shape = dic["layout"][0]
    cplex = {'R': False, 'C': True}[dic['nptype'][-1]]

    # read in the data
    data = read_sec(filename, dtype, shape, cplex)
    return dic, data


def read_lowmem(filename, par_file=None):
    """
    Read RNMRTK files with minimal memory usage

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file to read (.sec).
    par_file : str or None, optional
        Filename of RNMRTK parameter file. If None (default) a the last four
        characters of `file` are changed to .par.

    Returns
    -------
    dic : dic
        Dictionary of RNMRTK parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    Notes
    -----
    The dictionary parameters are ordered opposite the data layout, that is to
    say the the FIRST parameter in each list corresponds to the LAST axis in
    the data array.

    See Also
    --------
    read : Read RNMRTK files.
    write : Write RNMRTK files.

    """
    # determine par_file name if not given
    if par_file is None:
        par_file = filename[:-4] + ".par"
    dic = read_par(par_file)

    # determine shape, complexity and endiness from dictionary
    fshape = list(dic["layout"][0])
    cplex = {'R': False, 'C': True}[dic['nptype'][-1]]
    if cplex:
        fshape[-1] //= 2
    big = {'<': False, '>': True}[dic['format'][0]]
    data = rnmrtk_nd(filename, fshape, cplex, big)
    return dic, data


def write(filename, dic, data, par_file=None, overwrite=False):
    """
    Write RNMRTK files.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file to write to (.sec).
    dic : dict
        Dictionary of RNMRTK parameters.
    data : ndarray
        Array of NMR data.
    par_file : str or None, optional
        Filename of RNMRTK parameter file. If None (default) a the last four
        characters of `file` are changed to .par.
    overwrite : bool, optional
        True to overwrite existing files. False will raises a Warning if the
        file exists.

    See Also
    --------
    write_lowmem : Write RNMRTK files using minimal amounts of memory.
    read : Read RNMRTK files.

    """
    # determine par_file name if not given
    if par_file is None:
        par_file = filename[:-4] + ".par"
    write_par(par_file, dic, overwrite)
    dtype = dic["format"]
    write_sec(filename, data, dtype, overwrite)


def write_lowmem(filename, dic, data, par_file=None, overwrite=False):
    """
    Write RNMRTK files using minimal amounts of memory (trace by trace).

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file to write to (.sec).
    dic : dict
        Dictionary of RNMRTK parameters.
    data : array_like
        Array of NMR data.
    par_file : str or None, optional
        Filename of RNMRTK parameter file. If None (default) a the last four
        characters of `file` are changed to .par.
    overwrite : bool, optional
        True to overwrite existing files. False will raises a Warning if the
        file exists.

    See Also
    --------
    write : Write RNMRTK files using minimal amounts of memory.
    read_lowmem : Read RNMRTK files using minimal amounts of memory.

    """
    # determine par_file name if not given
    if par_file is None:
        par_file = filename[:-4] + ".par"
    write_par(par_file, dic, overwrite)

    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite=overwrite, mode='wb')

    # write out the file trace by trace
    for tup in np.ndindex(data.shape[:-1]):
        put_trace(f, data[tup])
    f.close()
    return

#######################
# sec reading/writing #
#######################


def write_sec(filename, data, dtype='f4', overwrite=False):
    """
    Write a RNMRTK .sec file.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file to write to (.sec).
    data : array_like
        Array of NMR data.
    dtype : dtype
        Data type to convert data to before writing to disk.
    overwrite : bool, optional
        True to overwrite existing files. False will raises a Warning if the
        file exists.

    See Also
    --------
    write : Write RNMRTK files.

    """
    # open file
    f = fileiobase.open_towrite(filename, overwrite, mode='wb')

    # interleave read/imag if needed
    if np.iscomplexobj(data):
        data = interleave_data(data)

    # write data and close file
    f.write(data.astype(dtype).tobytes())
    f.close()
    return


def read_sec(filename, dtype, shape, cplex):
    """
    Read a RNMRTK parameter .par file.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK (.sec) file to read .
    dtype : dtype
        Type of data in file, typically 'float32'.
    shape : tuple
        Shape of data.
    cplex : bool
        True if the last (fast) dimension is complex. False is real only.

    Returns
    -------
    data : ndarray
        Array of NMR data.

    """
    data = get_data(filename, dtype)
    data = data.reshape(shape)
    if cplex:
        data = uninterleave_data(data)
    return data

##########################
# data get/put functions #
##########################


def get_data(filename, dtype):
    """
    Get spectral data from a RNMRTK file.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file (.sec) to get data from.
    dtype : dtype
        Type of data in file, typically 'float32'

    Returns
    -------
    rdata : ndarray
        Raw NMR data, unshaped and typically not complex.

    """
    return np.fromfile(filename, dtype)


def get_trace(f, num_points, big):
    """
    Get a trace from an open RNMRTK file.

    Parameters
    -----------
    f : file object
        Open file object to read from.
    num_points : int
        Number of points in trace (R+I)
    big : bool
        True for data that is big-endian, False for little-endian.

    Returns
    -------
    trace : ndarray
        Raw trace of NMR data.

    """
    if big:
        bsize = num_points * np.dtype('>f4').itemsize
        return np.frombuffer(f.read(bsize), dtype='>f4')
    else:
        bsize = num_points * np.dtype('<f4').itemsize
        return np.frombuffer(f.read(bsize), dtype='<f4')


def put_trace(f, trace):
    """
    Put a trace to an open RNMRTK file.

    Parameters
    ----------
    f : file object
        Open file object to read from.
    trace : ndarray
        Raw trace of NMR data, may be complex64.

    """
    f.write(trace.view('float32').tobytes())


def uninterleave_data(data):
    """
    Remove interleaving of real. imag data in last dimension of data.
    """
    return data.view('complex64')


def interleave_data(data):
    """
    Interleave real, imag data in data
    """
    return data.view('float32')
    # size = list(data.shape)
    # size[-1] = size[-1]*2
    # data_out = np.empty(size,dtype="float32")
    # data_out[...,::2] = data.real
    # data_out[...,1::2] = data.imag
    # return data_out

######################
# low-memory objects #
######################


class rnmrtk_nd(fileiobase.data_nd):
    """
    Emulate a ndarray objects without loading data into memory for low memory
    reading of RNMRTK files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of RNMRTK file (.sec) to read.
    fshape : tuple of ints
        Shape of data in file.
    cplex : bool
        True if the last (fast) axis is complex.
    big : bool
        True for big-endian data, False for little-endian.
    order : tuple
        Ordering of axes against file. None for (0, 1, 2, ...).

    """

    def __init__(self, filename, fshape, cplex, big, order=None):
        """
        Create and set up
        """
        # check and set order
        if order is None:
            order = range(len(fshape))
        self.order = order

        # set additional parameters
        self.fshape = fshape    # shape on disk
        self.cplex = cplex
        self.filename = filename
        self.big = big

        if self.cplex:
            self.dtype = np.dtype('complex64')
        else:
            self.dtype = np.dtype('float32')

        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = rnmrtk_nd(self.filename, self.filename, self.fshape, self.cplex,
                      self.big, order)
        return n

    def __fgetitem__(self, slices):
        """
        Return ndarray of selected values.

        slices is a well formatted tuple of slices

        """
        # separate the last slice from the leading slices
        lslice = slices[-1]
        fslice = slices[:-1]

        # and the same for fshape
        lfshape = self.fshape[-1]
        ffshape = self.fshape[:-1]

        # find the output size and make an to/from nd_iterator
        osize, nd_iter = fileiobase.size_and_ndtofrom_iter(ffshape, fslice)
        osize.append(len(range(lfshape)[lslice]))

        # create an empty array to store the selected slices
        out = np.empty(tuple(osize), dtype=self.dtype)

        # opent the file for reading
        f = open(self.filename, 'rb')

        # read in the data trace by trace
        for out_index, in_index in nd_iter:

            # determine the trace number from the index
            ntrace = fileiobase.index2trace_flat(ffshape, in_index)

            # seek to the correct place in the file
            if self.cplex:
                ts = ntrace * lfshape * 2 * 4
                f.seek(ts)
                trace = get_trace(f, lfshape * 2, self.big)
                trace = uninterleave_data(trace)
            else:
                ts = ntrace * lfshape * 4
                f.seek(ts)
                trace = get_trace(f, lfshape, self.big)

            # put the trace into the output array
            out[out_index] = trace[lslice]

        # close the file and return
        f.close()
        return out

##################################
# Parameter dictionary utilities #
##################################


def find_dic_dim(dic, dim):
    """
    Find dimension in dictionary which corresponds to array dimension.

    Parameters
    ----------
    dic : dict
        Dictionary of RNMRTK parameters.
    dim : int, non-negative
        Dimension of data array.

    Returns
    -------
    ddim : int
        Dimension in dic which corresponds to array dimension, dim.

    """
    dic_dims = [int(i[1]) - 1 for i in dic['layout'][1]]
    return dic_dims.index(dim)


def find_array_dim(dic, ddim):
    """
    Find array dimension which corresponds to dictionary dimension.

    Parameters
    ----------
    dic : dict
        Dictionary of RNMRTK parameters.
    ddim : int, non-negative
        Dimension in dictionary.

    Returns
    -------
    dim : int
        Dimension in array which corresponds to dictionary dimension, ddim.

    """
    dic_dims = [int(i[1]) for i in dic['layout'][1]]
    return dic_dims[ddim]

############################
# parameter file functions #
############################


def read_par(filename):
    """
    Parse a RNMRTK parameter (.par) file.

    Parameters
    ----------
    file : str
        Filename of RNMRTK parameter file (.par) to read

    Returns
    -------
    dic : dict
        Dictionary of RNMRTK parameters.

    """
    dic = {}
    with open(filename, 'r') as f:
        for line in f:
            if len(line.split()) >= 2:
                parse_par_line(line, dic)

    # check that order and layout match, if they do remove from dictionary
    if dic['order'] != [int(i[1]) for i in dic['layout'][1]]:
        warn('Dom order and layout order do not match')
    else:
        dic.pop('order')

    return dic


def write_par(par_file, dic, overwrite):
    """
    Write a RNMRTK parameter file (.par).

    Parameters
    -----------
    par_file : str
        Filename of RNMRTK parameter file (.par) to write.
    dic : dict
        Dictionary of NMR parameters.
    overwrite : bool
        Set True to overwrite existing files, False will raise a Warning if the
        file exists.

    """
    # open file for writing
    f = fileiobase.open_towrite(par_file, overwrite, mode='w')

    # write comment line
    f.write('Comment \'' + dic['comment'] + '\'\n')

    # Dom line, set from layout
    l = "Dom " + " ".join(dic['layout'][1])
    f.write(l + "\n")

    # N line
    s = ["%14i %c" % (t) for t in zip(dic['npts'], dic['nptype'])]
    l = "N".ljust(8) + "".join(s)
    f.write(l + "\n")

    # write out additional lines Command    Value lines
    order = ['Sw', 'Xfirst', 'Xstep', 'Cphase', 'Lphase', 'Sf', 'Ppm', 'Nacq']
    codes = {'Sw': '%16.3f', 'Xfirst': '%16.5f', 'Xstep': '%16.5G',
             'Cphase': '%16.3f', 'Lphase': '%16.3f', 'Sf': '%16.2f',
             'Ppm': '%16.3f', 'Nacq': '%16i'}

    for lc in order:
        t = [codes[lc] % i for i in dic[lc.lower()]]
        l = lc.ljust(8) + "".join(t)
        f.write(l + "\n")

    # Quad line
    quad_dic = {'states': 'States', 'states-tppi': 'States-TPPI',
                'tppi': 'TPPI', 'tppi-redfield': 'TPPI-Redfield'}
    t = ["%16s" % (quad_dic[i]) for i in dic['quad']]
    l = "Quad".ljust(8) + "".join(t)
    f.write(l + "\n")

    # write format line
    if dic['format'][0] == '<':
        f.write('Format  Little-endian  IEEE-Float\n')
    else:
        f.write('Format  Big-endian     IEEE-Float\n')

    l = "Layout  " + " ".join([j + ":" + str(i) for i, j in
                               zip(*dic['layout'])])
    f.write(l + "\n")
    f.close()

    return


def parse_par_line(line, dic):
    """
    Parse a line from a RNMRTK parameter file (.par).
    """
    c, pl = line.split()[0], line.split()[1:]
    c = c.upper()

    if c == 'COMMENT':
        dic['comment'] = pl[0].strip('\'')

    elif c == 'DOM':
        dom = [s[0] for s in pl]        # dom as it appears in the file
        dic['ndim'] = ndim = len(pl)
        dic['order'] = order = [int(s[1]) for s in pl]  # dimension order
        # dom in accending order (to match other parameter)
        dic['dom'] = [dom[order.index(i)] for i in range(1, ndim + 1)]

    elif c == 'N':
        dic['npts'] = [int(i) for i in pl[::2]]
        dic['nptype'] = list(pl[1::2])

    elif c in ['SW', 'XFIRST', 'XSTEP', 'CPHASE', 'LPHASE', 'SF', 'PPM']:
        dic[c.lower()] = [float(i) for i in pl]
    elif c == "NACQ":
        dic['nacq'] = [int(i) for i in pl]
    elif c == "QUAD":
        dic['quad'] = [str(s).lower() for s in pl]

    # format assumes IEEE-Float type, only checks endiness
    elif c == 'FORMAT':
        if pl[0].upper() == "LITTLE-ENDIAN":
            dic['format'] = '<f4'
        else:
            dic['format'] = '>f4'
    elif c == 'LAYOUT':
        size = [int(p.split(":")[1]) for p in pl]
        domains = [p.split(":")[0] for p in pl]
        dic['layout'] = size, domains
    return
