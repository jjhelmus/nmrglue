"""
Functions for reading and writing Bruker binary (ser/fid) files, Bruker 
JCAMP-DX parameter (acqus) files, and Bruker pulse program (pulseprogram) 
files.

Bruker binary files (ser/fid) store data as an array of int32s whose endiness 
is determinded by the parameter BYTORDA (1 = big endian, 0 = little endian).
Typically the direct dimension is digitally filtered. The exact method of
removing this filter is unknown but an approximation is avaliable.

Bruker JCAMP-DX files (acqus, etc) are text file which are described by the 
`JCAMP-DX standard <http://www.jcamp-dx.org/>`_.  Bruker parameters are 
prefixed with a '$'.

Bruker pulseprogram files are text files described in various Bruker manuals.
Of special important are lines which describe exterenal variable assignments 
(surrounded by "'s), loops (begin with lo), phases (contain ip of dp) or 
increments (contain id, dd, ipu or dpu).  These lines are parsed when reading
the file.

"""

import numpy as np
import os
import fileiobase
from nmrglue.process import proc_base

# data creation

def create_data(data):
    """ 
    Create a bruker data array (recast into a complex128 or int32)
    """
    if np.iscomplexobj(data):
        return np.array(data,dtype='complex128')
    else:
        return np.array(data,dtype='int32')


# universal dictionary functions

def guess_udic(dic,data):
    """ 
    Guess parameters of universal dictionary from dic,data pair
    """
    # XXX if pprog, acqus are in dic use them 
    
    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in xrange(data.ndim):
        udic[i]["size"] = data.shape[i]

    return udic

def create_dic(udic):
    """ 
    Create a bruker dictionary from a universal dictionary
    """
    # determind the size in bytes
    if udic[udic["ndim"]-1]["complex"]:
        bytes = 8
    else:
        bytes = 4

    for k in xrange(udic["ndim"]): 
        bytes*=udic[k]["size"]

    dic= {"FILE_SIZE":bytes}

    return dic


# Global read/write function and related utilities

def read(dir=".",bin_file=None,acqus_files=None,pprog_file=None,shape=None,
         cplex=None,big=None):
    """ 
    Read Bruker files in directory

    Parameters:

    * dir           Directory to read from.
    * bin_file      Filename of binary file in directory.
    * acqus_files   List of filename(s) of acqus parameter files in directory.
    * pprog_file    Filename of pulseprogram in directory.
    * shape         Shape of resulting data (tuple).
    * cplex         Complexity of direct dimention (True/False).
    * big           Endianness of binary file. Set to True for big-endian, 
                    False for little-endian and None to determind automatically

    Only the dir parameter must be defined, others will be determined 
    automatically if not specified.

    """

    if os.path.isdir(dir) != True:
        raise IOError,"directory %s does not exist"%(dir)

    # determind parameter automatically
    if bin_file == None:
        if os.path.isfile(os.path.join(dir,"fid")):
            bin_file = "fid"
        elif os.path.isfile(os.path.join(dir,"ser")):
            bin_file = "ser"
        else:
            raise IOError,"no Bruker binary file could be found in %s"%(dir)

    if acqus_files == None:
        acqus_files = []
        for f in ["acqus","acqu2s","acqu3s","acqu4s"]:
            if os.path.isfile(os.path.join(dir,f)):
                acqus_files.append(f)

    if pprog_file == None:
        pprog_file = "pulseprogram"

    # create an empty dictionary
    dic = dict()

    # read the acqus_files and add to the dictionary
    for f in acqus_files:
        dic[f] = read_jcamp(os.path.join(dir,f))

    # read the pulse program and add to the dictionary
    dic["pprog"] = read_pprog(os.path.join(dir,pprog_file))

    # determind file size and add to the dictionary
    dic["FILE_SIZE"] = os.stat(os.path.join(dir,bin_file)).st_size

    # determind shape and complexity for direct dim if needed
    if shape == None or cplex == None:
        gshape,gcplex = guess_shape(dic)
    if shape == None:
        shape = gshape
    if cplex == None:
        cplex = gcplex
   
    # determind endianness (assume little-endian unless BYTORDA is 1)
    if big == None:
        big = False # default value
        if "acqus" in dic and "BYTORDA" in dic["acqus"]:
            if dic["acqus"]["BYTORDA"] == 1:
                big = True
            else:
                big = False

    # read the binary file
    f = os.path.join(dir,bin_file)
    null,data = read_binary(f,shape=shape,cplex=cplex,big=big)
    return dic,data


def read_lowmem(dir=".",bin_file=None,acqus_files=None,pprog_file=None,
                shape=None,cplex=None,big=None):
    """ 
    Read Bruker files using minimal amounts of memory 

    Parameters:

    * dir           Directory to read from.
    * bin_file      Filename of binary file in directory.
    * acqus_files   List of filename(s) of acqus parameter files in directory.
    * pprog_file    Filename of pulseprogram in directory.
    * shape         Shape of resulting data (tuple).
    * cplex         Complexity of direct dimention (True/False).
    * big           Endianness of binary file. Set to True for big-endian, 
                    False for little-endian and None to determind automatically

    Only the dir parameter must be defined, others will be determined 
    automatically if not specified.

    """

    if os.path.isdir(dir) != True:
        raise IOError,"directory %s does not exist"%(dir)

    # determind parameter automatically
    if bin_file == None:
        if os.path.isfile(os.path.join(dir,"fid")):
            bin_file = "fid"
        elif os.path.isfile(os.path.join(dir,"ser")):
            bin_file = "ser"
        else:
            raise IOError,"no Bruker binary file could be found in %s"%(dir)

    if acqus_files == None:
        acqus_files = []
        for f in ["acqus","acqu2s","acqu3s","acqu4s"]:
            if os.path.isfile(os.path.join(dir,f)):
                acqus_files.append(f)

    if pprog_file == None:
        pprog_file = "pulseprogram"

    # create an empty dictionary
    dic = dict()

    # read the acqus_files and add to the dictionary
    for f in acqus_files:
        dic[f] = read_jcamp(os.path.join(dir,f))

    # read the pulse program and add to the dictionary
    dic["pprog"] = read_pprog(os.path.join(dir,pprog_file))

    # determind file size and add to the dictionary
    dic["FILE_SIZE"] = os.stat(os.path.join(dir,bin_file)).st_size

    # determind shape and complexity for direct dim if needed
    if shape == None or cplex == None:
        gshape,gcplex = guess_shape(dic)
    if shape == None:
        shape = gshape
    if cplex == None:
        cplex = gcplex
   
    # determind endianness (assume little-endian unless BYTORDA is 1)
    if big == None:
        big = False # default value
        if "acqus" in dic and "BYTORDA" in dic["acqus"]:
            if dic["acqus"]["BYTORDA"] == 1:
                big = True
            else:
                big = False

    # read the binary file
    f = os.path.join(dir,bin_file)
    null,data = read_binary_lowmem(f,shape=shape,cplex=cplex,big=big)
    return dic,data



def write(dir,dic,data,bin_file=None,acqus_files=None,pprog_file=None,
    overwrite=False,big=True):
    """ 
    Write Bruker files 

    Parameters:

    * dir           Directory to write to.
    * dic           dictionary holding acqus_files and pprog_file parameters.
    * data          array of data
    * bin_file      Filename of binary file to write to in directory
    * acqus_files   Filename(s) of acqus files in directory to write to.
    * pprog_file    Filename of pulseprogram in directory.
    * overwrite     True to overwrite files, False to warn.
    * big           Endiness to write binary data with,
                    bigendian=True, little=False

    If any of bin_file,acqus_files or pprog_file are None the associated 
    file(s) will be determined automatically

    """

    # determind parameters automatically
    if bin_file == None:
        if data.ndim == 1:
            bin_file = "fid"
        else:
            bin_file = "ser"

    if acqus_files == None:
        acq = ["acqus","acqu2s","acqu3s","acqu4s"]
        acqus_files = [k for k in acq if dic.has_key(k)]

    if pprog_file == None:
        pprog_file = "pulseprogram"


    # write out the acqus files
    for f in acqus_files:
        write_jcamp(dic[f],os.path.join(dir,f),overwrite=overwrite)

    # write out the pulse program
    write_pprog(os.path.join(dir,pprog_file),dic["pprog"],overwrite=overwrite)

    # write out the binary data
    bin_full = os.path.join(dir,bin_file)
    write_binary(bin_full,dic,data,big=big,overwrite=overwrite)

    return


def guess_shape(dic):
    """ 
    Determind data shape and complexity from dictionary
    
    Returns: (shape,cplex)

    * shape   Tuple representing shape of the binary file.
    * cplex   Complexity of direct dimension.

    When dictionary does not contain enough information warning will be issued
    and (1),True is returned
    
    """

    # Clunkly and need some error checking for defaults, etc

    # check to see if we have needed dictionary keys
    if "pprog" not in dic or "acqus" not in dic or "FILE_SIZE" not in dic:
        print "Warning: Cannot determind shape do to missing dictionary keys"
        return (1),True

    # unpack the dictionary
    pd = dic["pprog"]
    acqus = dic["acqus"]

    # pprog error checks
    if "loop" not in pd or "incr" not in pd:
        print "Warning: Cannot determind shape do to bad pprog dictionary"
        return (1),True

    # unpack the pulse program dictionary
    lp = pd["loop"]
    ic = pd["incr"]

    # acqus dictionary checks
    if "AQ_mod" not in acqus or "TD" not in acqus:
        print "Warning: Cannot determind shape do to bad acqus dictionary"
        return (1),True

    # determind the minimum dimentionality from acqus files
    mindim=1
    if dic.has_key('acqu2s'):
        mindim=2
    if dic.has_key('acqu3s'):
        mindim=3
    if dic.has_key('acqu4s'):
        mindim=4

    # figure out the dim size predicted from size of loop
    # if greater than mindim set as ndim
    loopdim = [0,2,2,2,3,3,4,4][len(lp)]
    if loopdim>mindim:
        ndim = loopdim
    else:
        ndim = mindim
    
    # create a empty shape list
    shape = []

    # Determind X dimension size (round TD up to nearest 256)
    td = acqus["TD"]
    x = int( np.ceil(td/256.)*256 )

    # if Real (0) or Sequential (2) don't divide by two 
    if acqus["AQ_mod"] == 0 or acqus["AQ_mod"] == 2:
        shape.append(x)
        cplex = False
    else:
        shape.append(x/2)
        cplex = True

    # Determind Y dimension size if needed
    if ndim >= 2:
        y = 1       # default value
        if dic.has_key("acqu2s") and dic["acqu2s"].has_key("TD"):
            y = dic["acqu2s"]["TD"]
    
        if len(lp)==2 or len(lp)==4 or len(lp)==6:  # even number of loops
            if lp[0]==2 and len(ic[0])==0 and len(ic[1])!=0:
                y = 2*lp[1]
        
        if len(lp)%2==1:   # odd number of loops
            if lp[1]==2 and len(ic[0])==0 and len(ic[1])==0 and len(ic[2])!=0:
                y = 2*lp[2]

        shape.append(y)

    # Determind Z dimension size if needed
    if ndim >= 3:

        z = dic["FILE_SIZE"] / (x*y*4)

        if len(lp)==4 or len(lp)==6:
            if lp[2]==2 and len(ic[2])==0 and len(ic[3])!=0:
                z = 2*lp[3]
        
        if len(lp)==5 or len(lp)==7:
            if lp[3]==2 and len(ic[0])==0 and len(ic[3])==0 and len(ic[4])!=0:
                z = 2*lp[4]

        shape.append(z)

    # Determind A dimension size if needed
    if ndim >=4:
        
        a = dic["FILE_SIZE"] / (x*y*16*4)

        if len(lp)==6 and lp[4]==2 and len(ic[4])==0 and len(ic[5])!=0:
                a = 2*lp[5]
        
        if len(lp)==7:
            if lp[5]==2 and len(ic[0])==0 and len(ic[5])==0 and len(ic[6])!=0:
                a = 2*lp[6]

        shape.append(a)


    return tuple(shape[::-1]),cplex


# Bruker binary (fid/ser) reading and writing


def read_binary(filename,shape=(1),cplex=True,big=True):
    """ 
    Read Bruker binary data from file and return dic,data pair

    Parameters:

    * filename  Filename of Bruker binary file
    * shape     Tuple describing shape of resulting file
    * cplex     Flag indicating if direct dimension is complex
    * big       Endianness of binary file, True for big-endian, False for 
                little-endian

    If data cannot be reshaped 1D version of data will be returned

    """
    
    # open the file and get the data
    f = open(filename)
    data = get_data(f,big=big)

    # complexify if needed
    if cplex:
        data = complexify_data(data)

    # create dictionary
    dic = {"FILE_SIZE":os.stat(filename).st_size}

    # reshape if possible
    try:
        return dic,data.reshape(shape)

    except ValueError:
        print "Warning:",data.shape,"cannot be shaped into",shape
        return dic,data


def read_binary_lowmem(filename,shape=(1),cplex=True,big=True):
    """ 
    Read Bruker binary data from file using minimal memory.

    Parameters:

    * filename  Filename of Bruker binary file
    * shape     Tuple describing shape of resulting file
    * cplex     Flag indicating if direct dimension is complex
    * big       Endianness of binary file, True for big-endian, False for 
                little-endian

    Will issue warning if shape/cplex does not seem correct given the file
    size.

    """

    # create dictionary
    dic = {"FILE_SIZE":os.stat(filename).st_size}

    if len(shape) == 2:
        return dic,ser_2d(filename,shape=shape,cplex=cplex,big=big)
    elif len(shape) == 3:
        return dic,ser_3d(filename,shape=shape,cplex=cplex,big=big)
    elif len(shape) == 4:   # this does not work yet...
        return dic,ser_4d(filename,shape=shape,cplex=cplex,big=big)
    else:
        raise ValueError,"incorrect shape dimensionality, must be 2, 3 or 4"


def write_binary(filename,dic,data,overwrite=False,big=True):
    """ 
    Write Bruker binary data to file

    Parameters:

    * filename      Filename to write to.
    * dic           dictionary holding acqus_files and pprog_file parameters.
    * data          array of data
    * overwrite     True to overwrite files, False to warn.
    * big           Endiness to write binary data with,
                    bigendian=True, little=False

    """

    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    if type(data) == np.ndarray:
        if np.iscomplexobj(data):
            put_data(f,uncomplexify_data(data),big)
        else:
            put_data(f,data,big)
        
        f.close()
        return

    # we have a more complex data object, assume complex
    if data.ndim == 1:
        trace = data[:]
        if np.iscomplexobj(data):
            put_data(f,uncomplexify_data(trace),big)
        else:
            put_data(f,trace,big)

    elif data.ndim == 2:
        for trace in data:
            if np.iscomplexobj(data):
                put_data(f,uncomplexify_data(trace),big)
            else:
                put_data(f,trace,big)

    elif data.ndim == 3:
        for i in xrange(data.shape[0]):
            for j in xrange(data.shape[1]):
                trace = data[i,j]
                if np.iscomplexobj(trace):
                    put_data(f,uncomplexify_data(trace),big)
                else:
                    put_data(f,trace,big)

    else:   
        raise ValueError("unknown dimensionality")

    f.close()
    return


# ser_* objects

class ser_2d(fileiobase.data_2d):
    """
    set_2d emulates a numpy.ndarray object without loading data into memory

    * slicing operations returns ndarray objects.
    * can iterate over with expected results
    * transpose and swapaxes functions creates a new ser_2d object with the 
       new axes ordering
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,shape,cplex,big,order=["y","x"]):
        """create and set up a ser_2d object"""
        
        if len(shape) != 2:
            raise ValueError,"shape parameter is not a 2 tuple"

        # unpack shape
        lenY,lenX = shape

        # check if size correct
        
        if cplex == True:
            if os.stat(filename).st_size != lenY*lenX*4*2:
                raise IndexError,"invalid dimensions"
        else:
            if os.stat(filename).st_size != lenY*lenX*4:
                raise IndexError,"invalid dimensions"

        # open the file
        self.filename = filename
        self.f = open(filename)

        # sizes 
        self.lenX = lenX
        self.lenY = lenY

        # order
        self.order = order

        # shape based on order
        a = [self.lenY,self.lenX]
        self.shape = ( a[order.index("y")], a[order.index("x")])
        del(a)

        # complexity and endianness
        self.cplex = cplex
        self.big = big

        # dtype
        if self.cplex:
            self.dtype = np.dtype("complex128")
        else:
            self.dtype = np.dtype("int32")

        # dimensionality
        self.ndim = 2


    def __fcopy__(self,order):
        
        n = ser_2d(self.filename,(self.lenY,self.lenX),self.cplex,self.big,
                   order)
        return n

    def __fgetitem__(self,(sY,sX)):
        """returns ndarray of selected values

        (sY,sX) is a well formatted tuple of slices

        """

        # make the empty data set
        lenY = len(range(self.lenY)[sY])
        lenX = len(range(self.lenX)[sX])
        if self.cplex:
            data = np.empty( (lenY,lenX),dtype="complex128")
        else:
            if self.big:
                data = np.empty( (lenY,lenX),dtype=">i4")
            else:
                data = np.empty( (lenY,lenX),dtype="<i4")

        # read the data
        for jY,iY in enumerate(range(self.lenY)[sY]):
            
            if self.cplex:
                self.f.seek(iY*self.lenX*2*4)   # 2x4-byte int32 per cplex pt. 
                trace = get_trace(self.f,self.lenX*2,self.big)
                trace = complexify_data(trace)[sX]
            else:
                self.f.seek(iY*self.lenX*4)     # 1x4-byte int32 per real pt.
                trace = get_trace(self.f,self.lenX,self.big)

            data[jY] = trace

        return data


class ser_3d(fileiobase.data_3d):
    """
    set_3d emulates a numpy.ndarray object without loading data into memory

    * slicing operations returns ndarray objects.
    * can iterate over with expected results
    * transpose and swapaxes functions creates a new ser_3d object with the 
       new axes ordering
    * has ndim, shape, and dtype attributes.
    
    """

    def __init__(self,filename,shape,cplex,big,order=["z","y","x"]):
        """ create and set up a ser_3d object """

        if len(shape) !=3:
            raise ValueError,"shape parameter is not a 3 tuple"

        # unpack shape
        lenZ,lenY,lenX = shape

        # check if size correct
        if cplex == True:
            if os.stat(filename).st_size != lenZ*lenY*lenX*4*2:
                raise IndexError,"invalid dimentions"
        else:
            if os.stat(filename).st_size != lenZ*lenY*lenX*4:
                raise IndexError,"invalid dimentions"

        # open the file
        self.filename = filename
        self.f = open(filename)

        # sizes 
        self.lenX = lenX
        self.lenY = lenY
        self.lenZ = lenZ

        # order
        self.order = order

        # shape based on order
        a = [self.lenZ,self.lenY,self.lenX]
        self.shape = ( a[order.index("z")],a[order.index("y")],
                       a[order.index("x")])
        del(a)

        # complexity and endianness
        self.cplex = cplex
        self.big = big

        # dtype
        if self.cplex:
            self.dtype = np.dtype("complex128")
        else:
            self.dtype = np.dtype("int32")

        # dimensionality
        self.ndim = 3

    def __fcopy__(self,order):
        
        n = ser_3d(self.filename,(self.lenZ,self.lenY,self.lenX),self.cplex,
                   self.big,order)
        return n

    def __fgetitem__(self,(sZ,sY,sX)):
        """ returns ndarray of selected values

        (sZ,sY,sX) is a well formatted tuple of slices

        """

        # make the empty data set
        lenZ = len(range(self.lenZ)[sZ])
        lenY = len(range(self.lenY)[sY])
        lenX = len(range(self.lenX)[sX])

        if self.cplex:
            data = np.empty( (lenZ,lenY,lenX),dtype="complex128")
        else:
            if self.big:
                data = np.empty( (lenZ,lenY,lenX),dtype=">i4")
            else:
                data = np.empty( (lenZ,lenY,lenX),dtype="<i4")

        # read the data
        for jZ,iZ in enumerate(range(self.lenZ)[sZ]): 
            for jY,iY in enumerate(range(self.lenY)[sY]):
            
                if self.cplex:
                    # each complex point is 2 4-byte int32 
                    ts = iZ*self.lenY*self.lenX*2*4+iY*self.lenX*2*4
                    self.f.seek(ts)
                    trace = get_trace(self.f,self.lenX*2,self.big)
                    trace = complexify_data(trace)[sX]
                else:
                    # each real point is a 4-byte int32
                    ts = iZ*self.lenY*self.lenX*4+iY*self.lenX*4
                    self.f.seek(ts)
                    trace = get_trace(self.f,self.lenX,self.big)

                data[jZ,jY] = trace

        return data



# binary get/put functions

def get_data(file,big):
    """ 
    Get binary data from file object with given endiness
    """
    if big == True:
        return np.frombuffer(file.read(),dtype='>i4')
    else:
        return np.frombuffer(file.read(),dtype='<i4')

def put_data(file,data,big=True):
    """ 
    Put data to file object with given endiness
    """

    if big:
        file.write(data.astype('>i4').tostring())
    else:
        file.write(data.astype('<i4').tostring())
    
    return

def get_trace(file,num_points,big):
    """ 
    Get trace of num_points from file with given endiness
    """
    if big == True:
        bsize = num_points*np.dtype('>i4').itemsize
        return np.frombuffer(file.read(bsize),dtype='>i4')
    else:
        bsize = num_points*np.dtype('<i4').itemsize
        return np.frombuffer(file.read(bsize),dtype='<i4')


# data manipulation functions

def complexify_data(data):
	""" 
    Complexify data packed real,imag data
	"""
	return data[...,::2] + data[...,1::2]*1.j


def uncomplexify_data(data_in):
    """ 
    Uncomplexify data (pack real,imag) into a int32 array
    """
    size = list(data_in.shape)
    size[-1] = size[-1]*2

    data_out = np.empty(size,dtype="int32")

    data_out[...,::2]  = data_in.real
    data_out[...,1::2] = data_in.imag

    return data_out


# digital filter functions

# Table of points to frequency shift Bruker data to remove digital filter
# (Phase is 360 degrees * num_pts)
# This table is an 'un-rounded' version base on the table by
# W.M. Westler and F. Abildgaard's offline processing note, online at:
# http://www.boc.chem.uu.nl/static/local/prospectnd/dmx_digital_filters.html
# and the updated table with additional entries at:
# http://sbtools.uchc.edu/help/nmr/nmr_toolkit/bruker_dsp_table.asp

# The rounding in the above tables appear to be based on k / (2*DECIM)
# for example 2 : 44.75   = 44 + 3/4
#             4 : 66.625  = 66 + 5/8
#             8 : 68.563 ~= 68 + 9/16 = 68.5625
# Using this the un-rounded table was created by checking possible unrounded
# fracions which would round to those in the original table.

bruker_dsp_table = {
    10: { 
        2    : 44.75,
        3    : 33.5,
        4    : 66.625,
        6    : 59.083333333333333,
        8    : 68.5625,
        12   : 60.375,
        16   : 69.53125,
        24   : 61.020833333333333,
        32   : 70.015625,
        48   : 61.34375,
        64   : 70.2578125,
        96   : 61.505208333333333,
        128  : 70.37890625,
        192  : 61.5859375,
        256  : 70.439453125,
        384  : 61.626302083333333,
        512  : 70.4697265625,
        768  : 61.646484375,
        1024 : 70.48486328125,
        1536 : 61.656575520833333,
        2048 : 70.492431640625,
        },
    11: {
        2    : 46.,
        3    : 36.5,
        4    : 48.,
        6    : 50.166666666666667,
        8    : 53.25,
        12   : 69.5,
        16   : 72.25,
        24   : 70.166666666666667,
        32   : 72.75,
        48   : 70.5,
        64   : 73.,
        96   : 70.666666666666667,
        128  : 72.5,
        192  : 71.333333333333333,
        256  : 72.25,
        384  : 71.666666666666667,
        512  : 72.125,
        768  : 71.833333333333333,
        1024 : 72.0625,
        1536 : 71.916666666666667,
        2048 : 72.03125
        },
    12: {
        2    : 46. ,
        3    : 36.5,
        4    : 48.,
        6    : 50.166666666666667,
        8    : 53.25,
        12   : 69.5,
        16   : 71.625,
        24   : 70.166666666666667,
        32   : 72.125,
        48   : 70.5,
        64   : 72.375,
        96   : 70.666666666666667,
        128  : 72.5,
        192  : 71.333333333333333,
        256  : 72.25,
        384  : 71.666666666666667,
        512  : 72.125,
        768  : 71.833333333333333,
        1024 : 72.0625,
        1536 : 71.916666666666667,
        2048 : 72.03125
        },
    13: {
        2    : 2.75, 
        3    : 2.8333333333333333,
        4    : 2.875,
        6    : 2.9166666666666667,
        8    : 2.9375,
        12   : 2.9583333333333333,
        16   : 2.96875,
        24   : 2.9791666666666667,
        32   : 2.984375,
        48   : 2.9895833333333333,
        64   : 2.9921875,
        96   : 2.9947916666666667
        } 
    }


def remove_digital_filter(dic,data):
    """
    Remove the digial filter from Bruker data.

    Use rm_dig_filter to specify decim, dspfvs, and grpdly paraters.

    Parameters:

    * dic   Dictionary of Bruker parameters
    * data  array of data.

    Returns: array with digital filter removed

    """
    if 'acqus' not in dic:
        raise ValueError("dictionary does not contain acqus parameters") 
    
    if 'DECIM' not in dic['acqus']:
        raise ValueError("dictionary does not contain DECIM parameter")
    decim = dic['acqus']['DECIM']
    
    if 'DSPFVS' not in dic['acqus']:
        raise ValueError("dictionary does not contain DSPFVS parameter")
    dspfvs = dic['acqus']['DSPFVS']

    if 'GRPDLY' not in dic['acqus']:
        grpdly = 0
    else:
        grpdly = dic['acqus']['GRPDLY']

    return rm_dig_filter(data,decim,dspfvs,grpdly)


    
def rm_dig_filter(data,decim,dspfvs,grpdly=0):
    """
    Remove the digital filter from Bruker data.

    Use remove_digital_filter to find parameters from Bruker dictionary.

    Parameters:
    
    * data      Array of data.
    * decim     Decimation rate (Bruker DECIM parameter).
    * dspfvs    Firmware version (Bruker DSPFVS parameter).
    * grpdly    Group delay, when available.  (Bruker GRPDLY parameter).

    grpdly is not always needed, but when provided decim and dspfvs are 
    ignored.

    Returns: array with digital filter removed.

    """
    
    # This algorithm gives results similar but not exactly the same
    # as NMRPipe.  It was worked out by examining sample FID converted using
    # NMRPipe against spectra shifted with nmrglue's processing functions.  
    # When a frequency shifting with a fft first (fft->first order phase->ifft)
    # the middle of the fid nearly matches NMRPipe's and the difference at the
    # beginning is simply the end of the spectra reversed.  A few points at 
    # the end of the spectra are skipped entirely. 
    # -jjh 2010.12.01

    # The algorithm is as follows:
    # 1. FFT the data
    # 2. Apply a negative first order phase to the data.  The phase is 
    #    determined by the GRPDLY parameter or found in the DSPFVS/DECIM 
    #    loopup table.
    # 3. Inverse FFT  
    # (these first three steps are a frequency shift with a FFT first, fsh2)
    # 4. Round the applied first order phase up by two integers. For example
    #    71.4 -> 73, 67.8 -> 69, and 48 -> 50, this is the number of points
    #    removed from the end of the fid.
    # 5. If the size of the removed portion is greater than 6, remove the first
    #    6 points, reverse the remaining points, and add then to the beginning
    #    of the spectra.  If less that 6 points were removed, leave the FID 
    #    alone.
          
    
    if grpdly > 0:  # use group delay value if provided (not 0 or -1)
        phase = grpdly
    
    # determind the phase correction
    else:
        if dspfvs >= 14:    # DSPFVS greater than 14 give no phase correction.
            phase = 0.
        else:   # loop up the phase in the table
            if dspfvs not in bruker_dsp_table:
                raise ValueError("dspfvs not in lookup table")
            if decim not in bruker_dsp_table[dspfvs]:
                raise ValueError("decim not in lookup table")
            phase = bruker_dsp_table[dspfvs][decim]

    # and the number of points to remove (skip) and add to the beginning
    skip = int(np.floor(phase+2.))  # round up two integers
    add = int(max(skip-6,0))        # 6 less, or 0

    # DEBUG 
    #print "phase: %f, skip: %i add: %i"%(phase,skip,add)

    # frequency shift
    pdata = proc_base.fsh2(data,-phase)
    
    # add points at the end of the specta to beginning
    pdata[...,:add] = pdata[...,:add]+pdata[...,:-(add+1):-1]
    # remove points at end of spectra
    return pdata[...,:-skip]


# JCAMP-DX functions

def read_jcamp(filename):
    """ 
    Read a Bruker JCAMP-DX file into a dictionary

    Note: This is not a fully functional JCAMP-DX reader, it is only intended
    to read Bruker acqus (and similar) files

    Creates two special dictionary keys _coreheader and _comments

    Bruker parameter "$FOO" are extracted into strings, floats or lists
    and assigned to dic["FOO"]

    """
    dic = {"_coreheader":[],"_comments":[]} # create empty dictionary
    f = open(filename,'r')

    # loop until EOF
    while len(f.read(1)):
        
        f.seek(-1,os.SEEK_CUR)  # rewind 1 byte
        line = f.readline().rstrip()    # read a line

        if line[:6] == "##END=":
            #print "End of file"
            break

        elif line[:2] == "$$":
            dic["_comments"].append(line)

        elif line[:2] == "##" and line[2] != "$":
            dic["_coreheader"].append(line)

        elif line[:3] == "##$":
            key,value = parse_jcamp_line(line,f)
            dic[key] = value

        else:
            print "Warning: Extraneous line:",line

    return dic

def parse_jcamp_line(line,f):
    """ 
    Parse a single JCAMP-DX line
    
    Extract the Bruker parameter name and value from a line from a JCAMP-DX
    file.  This may entail reading additional lines from the fileobj f if the 
    parameter value extends over multiple lines.
    
    """

    # extract key= text from line
    key = line[3:line.index("=")]
    text = line[line.index("=")+1:].lstrip()

    if "<" in text:   # string
        while ">" not in text:      # grab additional text until ">" in string
            text = text+"\n"+f.readline().rstrip() 
        value = text.replace("<","").replace(">","")

    elif "(" in text: # array
        num = int(line[line.index("..")+2:line.index(")")])+1
        value = []
        rline = line[line.index(")")+1:]

        # extract value from remainer of line
        for t in rline.split():
            if "." in t or "e" in t:
                value.append(float(t))
            else:
                value.append(int(t))
        
        # parse additional lines as necessary
        while len(value) < num:
            nline = f.readline().rstrip()
            for t in nline.split():
                if "." in t or "e" in t:
                    value.append(float(t))
                else:
                    value.append(int(t))

    elif text == "yes":
        value = True

    elif text == "no":
        value = False

    else:   # simple value
        if "." in text or "e" in text: 
            value = float(text)
        else:
            value = int(text)

    return key,value


def write_jcamp(dic,filename,overwrite=False):
    """ 
    Write a Bruker JCAMP-DX file from a dictionary

    Written file will differ slightly from bruker's JCAMP-DX files in that all
    multi-value parameters will be written on multiple lines. Bruker is 
    inconsistent on what is written to a single line and what is not.  
    In addition line breaks may be slightly different but will always be 
    within JCAMP-DX specification.  Finally long floating point values
    may loose precision when writing.

    For example:

        ##$QS= (0..7)83 83 83 83 83 83 83 22
        
        will be written as
        ##$QS= (0..7)
        83 83 83 83 83 83 83 22

    """
 
    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    # create a copy of the dictionary
    d = dict(dic)

    # remove the comments and core header key from dictionary
    comments = d.pop("_comments")
    corehdr = d.pop("_coreheader")

    # write out the core headers
    for line in corehdr:
        f.write(line)
        f.write("\n")

    # write out the comments
    for line in comments:
        f.write(line)
        f.write("\n")

    keys = d.keys()
    keys.sort()

    # write out each key,value pair
    for key in keys:
       write_jcamp_pair(f,key,d[key])

    # write ##END= and close the file
    
    f.write("##END=")
    f.close()


def write_jcamp_pair(f,key,value):
    """ 
    Write out a line of a JCAMP file

    a 'line' might actually be more than one line for arrays
    """

    # the parameter name and such
    line = "##$"+key+"= "

    if type(value) == float or type(value) == int:  # simple numbers
        line = line + str(value)

    elif type(value) == str:        # string
        line = line+"<"+value+">"

    elif type(value) == bool:   # yes or no
        if value:
            line = line+"yes"
        else:
            line = line+"no"
        
    elif type(value) == list:   # lists
        # write out the current line
        line = line+"(0.."+str(len(value)-1)+")"
        f.write(line)
        f.write("\n")
        line = ""

        # loop over elements in value printing out lines when
        # they reach > 70 characters or the next value would cause
        # the line to go over 80 characters
        for v in value:     
            if len(line) > 70:
                f.write(line)
                f.write("\n")
                line = ""

            to_add = str(v)

            if len(line+" "+to_add) > 80:
                f.write(line)
                f.write("\n")
                line = ""

            if line != "":
                line = line+to_add+" "
            else:
                line = to_add+" "

    # write out the line and a newline character
    f.write(line)
    f.write("\n")

    return

# pulse program read/writing functions

def read_pprog(filename):
    """ 
    Parse a Bruker pulse program (pulseprogram) file for information

    Parameter:

    * filename  name of pulseprogram file to read from

    Returns a dictionary with following keys:

    * var       dictionary of variables assigned in pulseprogram
    * incr      list of lists containing increment times
    * loop      list of loop multipliers
    * phase     list of lists containing phase elements
    * ph_extra  list of lists containing comments at the end of phase lines

    The incr,phase and ph_extra lists match up with loop list.  For example 
    incr[0],phase[0] and ph_extra[0] are all increment and phase commands 
    with comments which occur during loop 0 which has loop[0] steps.

    """
    
    # open the file
    f = open(filename,'r')

    # initilize lists and dictionaries
    var = dict()
    loop = []
    incr = [ [] ]
    phase = [ [] ]
    ph_extra = [ [] ]

    # loop over lines in pulseprogram looking for loops, increment, 
    # assigments and phase commands
    for line in f:

        # split line into comment and text and strip leading/trailing spaces
        if ";" in line:
            comment = line[line.index(";"):]
            text = line[:line.index(";")].strip()
        else:
            comment = ""
            text = line.strip()

        # remove label from text when first word is all digits or 
        # has "," as the last element
        if len(text.split())!=0:
            s = text.split()[0]
            if s.isdigit() or s[-1] == ",":
                text = text[len(s):].strip()

        # skip blank lines and include lines
        if text == "" or text[0] == "#":
            #print line,"--Blank, Comment or Include"
            continue

        # see if we have quotes and have an assigment 
        # syntax "foo=bar"
        # add foo:bar to var dictionary
        if "\"" in text:
            if "=" in line:
                # strip quotes, split on = and add to var dictionary
                text = text.strip("\"")
                key,value = text.split("=")
                var[key]=value
                #print line,"--Assignment"
                continue
            else:
                #print line,"--Statement"
                continue

        # loops begin with lo
        # syntax is: lo to N time M
        # add M to loop list
        if text[0:2] == "lo":
            loop.append(text.split()[4])
            incr.append([])
            phase.append([])
            ph_extra.append([])
            #print line,"--Loop"
            continue

        # increment statement have id, dd, ipu or dpu
        # syntax foo {id/dd/ipu/dpu}N
        # store N to incr list
        if ("id" in text) or ("dd" in text):
            incr[len(loop)].append(int(text.split()[1][2:]))
            #print line,"--Increment"
            continue

        if ("ipu" in text) or ("dpu" in text):
            incr[len(loop)].append(int(text.split()[1][3:]))
            #print line,"--Increment"
            continue

        # phase statement have ip or dp 
        # syntax fpp {ip/dp}N extra
        # store N to phase list and extra to ph_extra list
        if ("ip" in text) or ("dp" in text):
            phase[len(loop)].append(int(text.split()[1][2:]))

            # find the first space after "ip" and read past there
            last = text.find(" ",text.index("ip"))
            if last == -1:
                ph_extra[len(loop)].append("")
            else:
                ph_extra[len(loop)].append(text[last:].strip())
            #print line,"--Phase"
            continue

        #print line,"--Unimportant"

    f.close()

    # remove the last empty incr, phase and ph_extra lists
    incr.pop()
    phase.pop()
    ph_extra.pop()

    # convert loop to numbers if possible
    for i,t in enumerate(loop):

        if t.isdigit():
            loop[i] = int(t)
        else:
            if var.has_key(t) and var[t].isdigit():
               loop[i] = int(var[t]) 

    # create the output dictionary
    dic = {"var":var,"incr":incr,"loop":loop,"phase":phase,"ph_extra":ph_extra}

    return dic

def write_pprog(filename,dic,overwrite=False):
    """ 
    Write a minimal Bruker pulse program 
    
    DO NOT TRY TO RUN THE RESULTING PULSE PROGRAM

    This pulse program should return the same dictionary when read using 
    read_pprog, nothing else.  The pulse program will be nonsense.
    
    """
    
    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite) 

    # write a comment
    f.write("; Minimal Bruker pulseprogram created by write_pprog\n")

    # write our the variables
    for k,v in dic["var"].iteritems():
        f.write("\""+k+"="+v+"\"\n")
    
    # write out each loop
    for i,steps in enumerate(dic["loop"]):
        
        # write our the increments
        for v in dic["incr"][i]:
            f.write("d01 id"+str(v)+"\n")
            
        # write out the phases
        for v,w in zip(dic["phase"][i],dic["ph_extra"][i]):
            f.write("d01 ip"+str(v)+" "+str(w)+"\n")

        f.write("lo to 0 times "+str(steps)+"\n")

    # close the file
    f.close()
        
    return
