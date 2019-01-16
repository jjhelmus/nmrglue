""" Unit tests for nmrglue/fileio/bruker.py module """

import os

import numpy as np
import nmrglue as ng
import shutil
import tempfile

# Test data.
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'bruker')

def test_read_pdata():
    """Reading processed bruker 1D data"""

    # read processed data
    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, '1', 'pdata', '1'))

    # check that the data is correct
    assert np.all(data == [-36275840.0, -34775104.0])

    # check dictionaries are correct
    assert len(dic.keys()) == 2
    assert 'procs' in dic.keys()
    assert 'acqus' in dic.keys()
    assert len(dic['procs'].keys()) == 131
    assert len(dic['acqus'].keys()) == 298


def test_reorder_submatrix():
    """reordering submatrix back and forth"""

    # make a dummy matrix
    data = np.arange(16, dtype='float64').reshape(4, 4)
    
    # correctly reordered matrix
    rdata = np.array([[ 0,  1,  4,  5],
                      [ 2,  3,  6,  7],
                      [ 8,  9, 12, 13],
                      [10, 11, 14, 15]], dtype='float64')

    # reorder from the submatrix form
    r1data = ng.bruker.reorder_submatrix(data, shape=(4, 4), 
                                         submatrix_shape=(2, 2), reverse=False)

    # reorder to the submatrix form
    r2data = ng.bruker.reorder_submatrix(r1data, shape=(4, 4), 
                                         submatrix_shape=(2, 2), reverse=True)
    # checks
    assert np.all(rdata == r1data)
    assert np.all(data == r2data)
    

def test_write_pdata():
    """ Writing a processed Bruker dataset """

    dic, data = ng.bruker.read_pdata(os.path.join(DATA_DIR, '1', 'pdata', '1'))
    
    # write to a temperory file
    td = tempfile.mkdtemp('.')
    ng.bruker.write_pdata(td, dic, data, write_procs=True, pdata_folder=10)

    assert os.path.isdir(os.path.join(td, 'pdata', '10'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', 'procs'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', 'proc'))
    assert os.path.isfile(os.path.join(td, 'pdata', '10', '1r'))

    rdic, rdata = ng.bruker.read_pdata(os.path.join(td, 'pdata', '10'))

    assert np.all(data == rdata)
    assert rdic['procs'].keys() == dic['procs'].keys()
    shutil.rmtree(td)
