"""
Functions for reading and writing spectral data to the glue format which
stores spectral data in the Hierarchical Data Format (HDF5).

glue files are HDF5 files with the spectral data stored in a dataset names
'spectrum' and any parameters stored in the dataset attributes.  At minimum
the parameter dictionary must contain a ndim key with the dimentionality of
the data and a dictionry for each axis numbered (0,1,2...) with the following
keys:

========    =====================================================
Key         Description
========    =====================================================
car         Carrier freqiency in Hz
complex     True for complex data, False for magnitude
encoding    type of encoding, 'states', 'tppi', 'direct', etc
freq        True for frequency domain data, False otherwise
label       Axis label, ('13C', etc)
obs         Observation frequency in MHz
size        Dimension size (R|I for direct dim, R+I for indirect)
sw          Spectral width in Hz
time        True for time domain data, False otherwise
========    =====================================================

"""

import numpy as np
import h5py
from . import fileiobase


# unit conversion functions
def make_uc(dic, data, dim=-1):
    """
    make a unit conversion object
    """
    if dim == -1:
        dim = data.ndim - 1  # last dimention

    size = dic[dim]["size"]
    cplex = dic[dim]["complex"]
    sw = dic[dim]["sw"]
    obs = dic[dim]["obs"]
    car = dic[dim]["car"]

    return fileiobase.unit_conversion(size, cplex, sw, obs, car)


# dictionary/data creation functions (null functions)
def create_data(data):
    """
    Create glue data array (return data)
    """
    return data


def guess_udic(dic, data):
    """
    Guess parameters of a universal dictionary from dic,data pair
    """
    return dic


def create_dic(udic):
    """
    Create a glue dictionary from a universal dictionary
    """
    return dic


# file reading/writing functions
def read(filename, dataset="spectrum"):
    """
    Read a glue file
    """
    dic, data = read_lowmem(filename, dataset)
    return dic, data[:]


def read_lowmem(filename, dataset="spectrum"):
    """
    Read a glue file using mimimal memory usage
    """
    f = h5py.File(filename, 'r')
    dic = get_dic(f, dataset)
    data = wrap_data(f[dataset])
    return dic, data


def write(filename, dic, data, dataset="spectrum", overwrite=False):
    """
    Write dic,data pair to a HDF5 file
    """
    # create the file
    f = fileiobase.open_towrite(filename, overwrite=overwrite)
    f.close()
    f = h5py.File(filename, 'w')

    # write the dictionary and data to the file
    f.create_dataset(dataset, data=data)
    put_dic(f, dic, dataset)

    f.close()
    return


# data wrapping functions
def wrap_data(data):
    """
    wrap h5py.highlevel.Dataset objects into a more numpy like object
    """
    if len(data.shape) == 1:
        return data
    elif len(data.shape) == 2:
        return glue_2d(data)
    elif len(data.shape) == 3:
        return glue_3d(data)
    else:
        return data


# dictionary get/put
def get_dic(f, dataset="spectrum"):
    """
    Get a dictionary from dataset in a HDF5 File
    """

    # select the data set
    dset = f[dataset]

    dic = {}
    # loop over the attributes
    for key, value in dset.attrs.iteritems():

        if "_" in key:
            # we have an axis key
            axis, subkey = key.split("_", 1)
            axis = int(axis)
            if axis not in dic:
                dic[axis] = {}
            dic[axis][subkey] = value
        else:
            dic[key] = value

    return dic


def put_dic(f, dic, dataset="spectrum"):
    """
    Put a dictionary to the dataset in a HDF5 File
    """

    # select the data set
    dset = f[dataset]

    for key, value in dic.iteritems():

        # axis dictionaries
        if type(key) == int and type(value) == dict:
            axis = key
            for axiskey, axisvalue in value.iteritems():
                fullkey = str(axis) + "_" + axiskey
                dset.attrs[fullkey] = axisvalue

        # simple value
        else:
            dset.attrs[key] = value
    return


# glue_* objects
class glue_2d(fileiobase.data_nd):
    """
    glue_2d emulates numpy.ndarray objects without loading data into memory

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes functions create a new fid_2d object with the
      new axes ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self, Dataset, order=["y", "x"]):
        """
        Create and set up object
        """
        # set up array attributes
        self.Dataset = Dataset
        self.lenX = int(Dataset.shape[1])
        self.lenY = int(Dataset.shape[0])
        a = [self.lenY, self.lenX]
        self.order = order
        self.shape = tuple([a[order.index(k)] for k in ["y", "x"]])
        self.ndim = 2
        self.dtype = Dataset.dtype

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = glue_3d(self.Dataset, order)
        return n

    def __fgetitem__(self, (sY, sX)):
        """
        Return ndarray of selected values

        (sY,sX) is a well formated tuple of slices
        """
        return self.Dataset[sY, sX]


class glue_3d(fileiobase.data_nd):
    """
    glue_3d emulates numpy.ndarray objects without loading data into memory

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes functions create a new fid_2d object with the
      new axes ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self, Dataset, order=["z", "y", "x"]):
        """
        Create and set up object
        """
        # set up array attributes
        self.Dataset = Dataset
        self.lenX = int(Dataset.shape[2])
        self.lenY = int(Dataset.shape[1])
        self.lenZ = int(Dataset.shape[0])
        a = [self.lenZ, self.lenY, self.lenX]
        self.order = order
        self.shape = tuple([a[order.index(k)] for k in ["z", "y", "x"]])
        self.ndim = 3
        self.dtype = Dataset.dtype

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = glue_3d(self.Dataset, order)
        return n

    def __fgetitem__(self, (sZ, sY, sX)):
        """
        Return ndarray of selected values

        (sZ, sY, sX) is a well formated tuple of slices
        """
        return self.Dataset[sZ, sY, sX]
