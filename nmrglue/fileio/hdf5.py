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

class JDataset(h5py.Dataset):
    """
    Dataset subclass which behaves more like a Numpy ndarray
    """

    pass
