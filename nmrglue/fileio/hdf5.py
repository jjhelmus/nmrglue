"""
Functions for reading and writing spectral data to the Hierarchical Data Format
(HDF5).  

"""

import numpy as np
import h5py
import fileiobase


def read(filename,dataset="spectra"):
    """
    Read a HDF5 file
    """
    f = h5py.File(filename,'r')
    dic = get_dic(f,dataset)
    data = f[dataset]
    return dic,data

def write(filename,dic,data,dataset="spectra",overwrite=False):
    """
    Write dic,data pair to a HDF5 file
    """
    # create the file
    f = fileiobase.open_towrite(filename,overwrite=overwrite)
    f.close()
    f = h5py.File(filename,'w')

    # write the dictionary and data to the file
    f.create_dataset(dataset,data=data)
    put_dic(f,dic,dataset)

    f.close()
    return


# dictionary get/put

def get_dic(f,dataset="spectra"):
    """
    Get a dictionary from dataset in a HDF5 File
    """
    
    # select the data set
    dset = f[dataset]

    dic = {}
    # loop over the attributes
    for key,value in dset.attrs.iteritems():
        
        if "_" in key:
            # we have an axis key
            axis,subkey = key.split("_",1)
            axis = int(axis)
            if axis not in dic:
                dic[axis] = {}
            dic[axis][subkey] = value
        else:
            dic[key] = value

    return dic

def put_dic(f,dic,dataset="spectra"):
    """
    Put a dictionary to the dataset in a HDF5 File
    """

    # select the data set
    dset = f[dataset]

    for key,value in dic.iteritems():
        
        # axis dictionaries
        if type(key) == int and type(value) == dict:
            axis = key
            for axiskey,axisvalue in value.iteritems():
                fullkey = str(axis)+"_"+axiskey
                dset.attrs[fullkey] = axisvalue

        # simple value
        else:
            dset.attrs[key] = value
    return



class TDataset(object):
    """ 
    h5py object which supports transposes and swap axes methods
    """

    def __init__(self,Dataset,order=None):

        self.Dataset = Dataset
        
        # set array like properties
        self.dshape = Dataset.shape
        self.dtype = Dataset.dtype 
        self.ndim = len(Dataset.shape)

        # order specifies how the axes match up to the Dataset axes
        # self.or
        if order == None:     # set order as normal if not specified
            self.order = tuple(range(self.ndim))
        else:
            self.order = order

        # shape from order and Dataset shape
        self.shape = tuple([self.Dataset.shape[i] for i in self.order]) 


    def __getitem__(self, args):
        """ 
        Read a slice from the HDF5 dataset.

        Takes slices and recarray-style field names (more than one is
        allowed!) in any order.  Obeys basic NumPy rules, including
        broadcasting.

        Also supports:

        * Boolean "mask" array indexing
        * Advanced dataspace selection via the "selections" module
        """

        args = args if isinstance(args, tuple) else (args,)

        # Sort field indices from the rest of the args.
        names = tuple(x for x in args if isinstance(x, str))
        args = tuple(x for x in args if not isinstance(x, str))

        # Create NumPy datatype for read, using only the named fields
        # as specified by the user.
        basetype = self.Dataset.dtype
        if len(names) == 0:
            new_dtype = basetype
        elif basetype.names is None:
            raise ValueError("Field names only allowed for compound types")
        else:
            for name in names:
                if not name in basetype.names:
                    raise ValueError("Field %s does not appear in this type." % name)
                new_dtype = numpy.dtype([(name, basetype.fields[name][0]) for name in names])

        # Perform the selection in the TDataset order.
        sel = h5py.selections.select(self.shape, args)

        if sel.nselect == 0:
            return numpy.ndarray((0,), dtype=new_dtype)

        if not isinstance(sel,h5py.selections.SimpleSelection):
            raise TypeError("Only rectangular (regular) selections supported")

        # unpack the selection
        start,count,step,scalar = sel._sel
        
        # reorder these into Dataset order
        dstart = tuple([start[i] for i in self.order])
        dcount = tuple([count[i] for i in self.order])
        dstep  = tuple([step[i] for i in self.order])
        dscalar = tuple([scalar[i] for i in self.order])

        # determind Dataset order mshape
        dmshape = tuple(x for x, y in zip(dcount, dscalar) if not y)

        # create a SpaceID object representing Dataset order selection
        id = h5py.h5s.create_simple(self.Dataset.shape, 
                    (h5py.h5s.UNLIMITED,)*len(self.Dataset.shape))
        id.select_hyperslab(dstart,dcount,dstep)

        # Create the output array to hold Dataset ordered selection
        arr = np.ndarray(dmshape, new_dtype, order='C')

        # This is necessary because in the case of array types, NumPy
        # discards the array information at the top level.
        mtype = h5py.h5t.py_create(new_dtype)

        # Perfom the actual read
        mspace = h5py.h5s.create_simple(dmshape)
        fspace = id
        self.Dataset.id.read(mspace, fspace, arr, mtype)

        # swap axes to put into TDataset ordering
        arr = arr.transpose(self.order)

        # Patch up the output for NumPy
        if arr.shape == ():
            arr = np.asscalar(arr)
        return arr
