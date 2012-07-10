""" 
Fuctions for reading and writing Rowland NMR Toolkit (RNMRTK) files 
"""

__developer_info__ = """
Information of the Rowland NMR Toolkit file format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

import numpy as np

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
    if dim < 0:   # negative dimensions
        dim = data.ndim + dim
    size = data.shape[dim]  # R|I
    cplx = {1:False, 2:True}[dic['nptype'][dim]]
    sw = dic['sw'][dim]   # Hz
    obs = dic['sf'][dim]   # MHz
    car = dic['zero_freq_ppm'][dim] * dic[obs]   # Hz
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
    udic = fileiobase.create_black_udic(data.ndim)

    # update the default values
    # TODO: write this 
    return udic

def create_dic(udic):
    """
    Create a RNMRTK dictionary from a universal dictionary.

    Parameters
    ----------
    udic : dict
        Universal dictionary of spectral parameters.

    Returns
    --------
    dic : dict
        Dictionary of RNMRTK parameters.

    """
    dic = make_empty_dic()

    # TODO
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
    cplex = {1:False, 2:True}[dic['nptype'][0]]

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
    fshape = dic["layout"][0]
    cplex = {1:False, 2:True}[dic['nptype'][0]]
    big = {'<':False, '>':True}[dic['format'][0]]
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
    f = fileiobase.open_towrite(filename, overwrite=overwrite)
   
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
    f = fileiobase.open_towrite(filename, overwrite)

    # interleave read/imag if needed
    if np.iscomplexobj(data):
        data = interleave_data(data)

    # write data and close file
    f.write(data.astype(dtype).tostring())
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
    f.write(trace.view('float32').tostring())

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
    #size = list(data.shape)
    #size[-1] = size[-1]*2
    #data_out = np.empty(size,dtype="float32")
    #data_out[...,::2] = data.real
    #data_out[...,1::2] = data.imag
    #return data_out

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
        if order == None:
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
        # seperate the last slice from the leading slices
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


############################
# parameter file functions #
############################

def find_shape(dic):
    """ 
    Determine the data shape (R+I for all dims) from a RNMRTK dictionary.
    """
    s = [dic['npts'][i] * dic['nptype'][i] for i in range(dic['ndim'])]
    return tuple(s[::-1])

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
    dic = make_empty_dic()
    f = open(filename, 'rb')
    for line in f:
        if len(line.split()) >= 2:
            parse_par_line(line, dic)
    # reorder dictionary lists per layout
    ndim = dic['ndim']
    plist = [int(s[1]) - 1 for s in dic['layout'][1][::-1]] + range(ndim, 4)
    permute_dic(dic, plist)
    return dic

def permute_dic(dic, plist):
    """
    Permute all parameters in dictionary according to plist (4 element list)
    """
    lists_keys = ['dom', 'ldim', 'nacq', 'npts', 'nptype', 'p0', 'p1', 'quad',
                    'sf', 'sw', 'xfirst', 'xstep', 'zero_freq_ppm']
    for key in lists_keys:
        dic[key] = [dic[key][i] for i in plist]
    return 

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
    # permute the dictionary so that the ordering is 1,2,...0,0
    plist = [0, 1, 2, 3]
    for i, v in enumerate(dic['ldim']):
        if v != 0:
            plist[i] = v - 1
    permute_dic(dic, plist)

    # open file for writing
    f = fileiobase.open_towrite(par_file, overwrite)
    
    # write comment line
    f.write('Comment \'' + dic['comment'] + '\'\n')
    
    # Dom line, set from layout
    ndim = dic['ndim'] 
    #l = "Dom "+" ".join([dic['dom'][i]+str(i+1) for i in range(ndim)])
    l = "Dom " + " ".join(dic['layout'][1])
    f.write(l + "\n")
    
    # N line
    npts = dic['npts']
    str_nptype = [['', 'R', 'C'][i] for i in dic['nptype']]
    t = [ "%14i %c" % (npts[i], str_nptype[i]) for i in range(ndim)]
    l = "N".ljust(8) + "".join(t)
    f.write(l + "\n")

    # write out additional lines Command    Value lines
    order = ['Sw', 'Xfirst', 'Xstep', 'Cphase', 'Lphase', 'Sf', 'Ppm', 'Nacq']
    codes = {'Sw':'%16.3f', 'Xfirst':'%16.5f', 'Xstep':'%16.5G',
             'Cphase':'%16.3f', 'Lphase':'%16.3f', 'Sf':'%16.2f', 
             'Ppm':'%16.3f', 'Nacq':'%16i'}
     
    c2d = {'Sw':'sw', 'Xfirst':'xfirst', 'Xstep':'xstep', 'Cphase':'p0',
           'Lphase':'p1', 'Sf':'sf', 'Ppm':'zero_freq_ppm', 'Nacq':'nacq'}

    for lc in order:
        t = [codes[lc] % (dic[c2d[lc]][i]) for i in range(ndim)]
        l = lc.ljust(8) + "".join(t)
        f.write(l + "\n")
    
    # Quad line
    str_quad = [{0:'States', 1:'States-TPPI'}[i] for i in dic['quad']]
    t = ["%16s" % (str_quad[i]) for i in range(ndim)]
    l = "Quad".ljust(8) + "".join(t)
    f.write(l + "\n")

    # write format line
    if dic['format'][0] == '<':
        f.write('Format  Little-endian  IEEE-Float\n')
    else:
        f.write('Format  Big-endian     IEEE-Float\n')

    l = "Layout  " + " ".join([ j+":"+str(i) for i, j in zip(*dic['layout'])])
    f.write(l + "\n")
    f.close()

    # permute the dictionary back
    permute_dic(dic, plist)
    return

def parse_par_line(line, dic):
    """
    Parse a line from a RNMRTK parameter file (.par).
    """
    c, pl = line.split()[0], line.split()[1:]
    c = c.upper()

    if c == 'COMMENT':
        dic['comment'] = pl[0].strip('\'')
    
    if c == 'DOM':
        for i, p in enumerate(pl):
            dic['dom'][int(p[1]) - 1] = p[0]
        # also update ndim and ldim keys
        dic['ndim'] = ndim = len(pl)
        dic['ldim'][:ndim] = range(1, ndim + 1) # since we have not permuted

    if c == 'N':
        for i, p in enumerate(pl[::2]):
            dic['npts'][i] = int(p)
        for i, p in enumerate(pl[1::2]):
            dic['nptype'][i] = {'C':2, 'R':1}[p]

    float_p = { 'SW':'sw', 'XFIRST':'xfirst', 'XSTEP':'xstep', 'CPHASE':'p0',
                'LPHASE':'p1', 'SF':'sf', 'PPM':'zero_freq_ppm'}

    if c in float_p.keys():
        for i, p in enumerate(pl):
            dic[float_p[c]][i] = float(p)
    elif c == "NACQ":
        for i, p in enumerate(pl):
            dic['nacq'][i] = int(p)
    elif c == "QUAD":
        for i, p in enumerate(pl):
            dic['quad'][i] = {'States':0, 'States-TPPI':1}[p]

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

def make_empty_dic():
    """
    Make an empty RNMRTK parameter dictionary.
    """
    dic = { 'pbufptr':0,
            'cbufptr':0,
            'npages':0,
            'pbuf':0,
            'ndim':0,
            'ldim':[0, 0, 0, 0],
            'npts':[1, 1, 1, 1],
            'nptype':[1, 1, 1, 1],
            'nacq':[0, 0, 0, 0],
            'sw':[4000.0, 4000.0, 4000.0, 4000.0],
            'xfirst':[0.0, 0.0, 0.0, 0.0],
            'xstep':[1.0/4000.0, 1/4000.0, 1/4000.0, 1/4000.0],
            'p0':[0.0, 0.0, 0.0, 0.0],
            'p1':[0.0, 0.0, 0.0, 0.0],
            'zero_freq_ppm':[4.77, 4.77, 4.77, 4.77],
            'sf':[399.65, 399.65, 399.65, 399.65], 
            'quad':[0, 0, 0, 0],
            'dom':['T', 'T', 'T', 'T'],
            'comment':'',
            'format':'<f4',
            'layout':([], [])
            }
    return dic
