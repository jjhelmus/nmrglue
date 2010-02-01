"""
Functions for reading and writing varian binary (fid) files and varian 
parameter (procpar) files.

Both varian file formats are well documented in:

* Varain MR News 2005-04-18 Importing Binary VnmrJ / VNMR FIDs into Third Party Software and VnmrJ / VNMR FID Data Format
* VnmrJ User Programming - Chapter 5: Parameters and Data

These are available (as of 08/2009) online from 
`Varian <http://varianinc.com>`_.

"""

import numpy as np
import struct
import os
import fileiobase

# data creation

def create_data(data):
    """ 
    Create a varian data array (recast into complex64 array)
    """
    return np.array(data,dtype="complex64")


# universal dictionary functions

def guess_udic(dic,data):
    """ 
    Guess parameter of universal dictionary from dic,data pair
    """
    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in xrange(data.ndim):
        udic[i]["size"] = data.shape[i]
    
    return udic


def create_dic(udic):
    """ 
    Create a varian dictionary from a universal dictionary
    """
    ddic = udic[udic["ndim"]-1]   # direct dimension dictionary
    
    # number of points in direct dimension is R+I
    if ddic["complex"]:
        num_points = ddic["size"]*2
    else:
        num_points = ddic["size"]

    # number of blocks is product of R+I of all other dimensions 
    num_blocks = 1
    for f in [udic[k]["size"] for k in xrange(udic["ndim"]-1)]:
        num_blocks = num_blocks*f

    # create the varian dictionary
    dic = dict()

    # default flag bits
    dic["S_DATA"]   = 1
    dic["S_SPEC"]   = 0
    dic["S_32"]     = 0
    dic["S_FLOAT"]  = 1
    dic["S_COMPLEX"]= 0
    dic["S_HYPERCOMPLEX"] = 0

    dic["S_ACQPAR"] = 1
    dic["S_SECND"]  = 0
    dic["S_TRANSF"] = 0
    dic["S_NP"]     = 0
    dic["S_NF"]     = 0
    dic["S_NI"]     = 0
    dic["S_NI2"]    = 0

    # parameters
    dic["nblocks"]   = num_blocks
    dic["ntraces"]   = 1
    dic["np"]        = num_points
    dic["ebytes"]    = 4
    dic["tbytes"]    = 4*num_points
    dic["bbytes"]    = dic["tbytes"]+28
    dic["vers_id"]   = 0
    dic["status"]    = 201  # 6th bit set
    dic["nbheaders"] = 1

    return dic


# varian binary (fid) file reading/writing

def write_fid(filename,dic,data,overwrite=False,repack=False):
    """ 
    Write a varian binary (fid) file

    Parameters:

    * filename  name of file to write to
    * dic       varian dictionary
    * data      data to write
    * overwrite True/False to overwrite existing file
    * repack    True/False to repack file and block headers

    When present dic['blockheader'] will be used for blockheaders.  When 
    absent a standard blockheader will be used.

    """

    if data.ndim == 1:
        write_fid_1D(filename,dic,data,overwrite=False,repack=False)
    elif data.ndim == 2:
        write_fid_2D(filename,dic,data,overwrite=False,repack=False)
    elif data.ndim == 3:
        write_fid_3D(filename,dic,data,overwrite=False,repack=False)
    else:
        raise ValueError("unsupported dimensionality") 
    
    return

def write_fid_3D(filename,dic,data,overwrite=False,repack=False):
    """ 
    Write a 3D varian binary (fid) file
    """
    
    # error checking
    if data.shape[2] != (dic["np"]/2):
        print StandardError("Warning: data and np size mismatch")
    if data.shape[0]*data.shape[1] != dic["nblocks"]:
        print StandardError("Warning: data and block size mismatch")

    # open file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    if repack:
        dic = repack_fileheader(dic)

    # write the fileheader to file
    put_fileheader(f,dic2fileheader(dic))
 
    # determind data type
    dt = find_dtype(dic)

    lenY = data.shape[1]
    lenZ = data.shape[0]

    if dic.has_key("blockheader") and len(dic["blockheader"])==data.shape[0]:
        for i in xrange(dic["nblocks"]):
            if repack:
                bh=dic2blockheader(repack_blockheader(dic["blockheader"][0]))
            else:
                bh = dic2blockheader(dic["blockheader"][0])
            
            y = int( np.floor(i/2.)%lenY)
            z = int( np.floor( (i-y*2)/lenY ) + i%2 )
            trace = np.array(interleave_data(data[z,y]),dtype=dt)
            put_block(f,dic,trace,bh)
        pass    # end of for loop

    else:   # create a generic blockheader
        bh = dic2blockheader(make_blockheader(dic,1))
        for i in xrange(dic["nblocks"]):
            bh[2] = int(i+1)    # set the blockheader index correctly
            y = int( np.floor(i/2.)%lenY)
            z = int( np.floor( (i-y*2)/lenY ) + i%2 )
            trace = np.array(interleave_data(data[z,y]),dtype=dt)
            put_block(f,dic,trace,bh)
        pass # end of for loop

    f.close()

    return

def write_fid_2D(filename,dic,data,overwrite=False,repack=False):
    """ 
    Write a 2D varian binary (fid) file
    """
    
    # error checking
    if data.shape[1] != (dic["np"]/2):
        print StandardError("Warning: data and np size mismatch")
    if data.shape[0] != dic["nblocks"]:
        print StandardError("Warning: data and block size mismatch")

    # open file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    if repack:
        dic = repack_fileheader(dic)

    # write the fileheader to file
    put_fileheader(f,dic2fileheader(dic))
 
    # determind data type
    dt = find_dtype(dic)

    if dic.has_key("blockheader") and len(dic["blockheader"])==data.shape[0]:
        for i in xrange(dic["nblocks"]):
            if repack:
                bh=dic2blockheader(repack_blockheader(dic["blockheader"][0]))
            else:
                bh = dic2blockheader(dic["blockheader"][i])

            trace = np.array(interleave_data(data[i]),dtype=dt)
            put_block(f,dic,trace,bh)
        pass    # end of for loop

    else:   # create a generic blockheader
        bh = dic2blockheader(make_blockheader(dic,1))
        for i in xrange(dic["nblocks"]):
            bh[2] = int(i+1)
            trace = np.array(interleave_data(data[i]),dtype=dt)
            put_block(f,dic,trace,bh)
        pass # end of for loop

    f.close()

    return
       
def write_fid_1D(filename,dic,data,overwrite=False,repack=False):
    """ 
    Write a 1D varian binary (fid) file
    """

    # error checking
    if data.shape[0] != (dic["np"]/2):
        print StandardError("Warning: data and np size mismatch")

    # open file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    if repack:
        dic = repack_fileheader(dic)

    # write the fileheader to file
    put_fileheader(f,dic2fileheader(dic))

    
    # determind data type
    dt = find_dtype(dic)

    if dic.has_key("blockheader"):  # dictionary has a blockheader use it
        if repack:
            bh = dic2blockheader(repack_blockheader(dic["blockheader"][0]))
        else:
            bh = dic2blockheader(dic["blockheader"][0])
        trace = np.array(interleave_data(data),dtype=dt)
        put_block(f,dic,trace,bh)
    else:
        bh = dic2blockheader(make_blockheader(dic,1))
        trace = np.array(interleave_data(data),dtype=dt)
        put_block(f,dic,trace,bh)

    f.close()
    return

def read_fid(filename,read_blockhead=False,(lenZ,lenY)=(None,None)):
    """ 
    Read a varian binary (fid) file returning a dic,data pair.

    When read_blockhead is True will read blockheader(s).

    For 3D data sets the length of the two indirect dimensions must be provided
    as lenZ and lenY.  If not provided the data set is open as if it was a 2D
    data set.  read_blockhead is ignored in 3D data sets.

    The data dtype will be complex64 or complex128 so that there is no loss of 
    precision in the conversions from int16/int32.

    """

    # open the file
    f = open(filename)

    # read the fileheader
    dic = fileheader2dic(get_fileheader(f))

    # check for 3D data sets
    if lenZ!=None and lenY!=None:
        f.close()
        # set up a low memory data object and then request all items
        dic,data = read_fid_lowmem_3D(filename,(lenZ,lenY))
        return dic,data[:,:,:]
        
    # read the data
    if read_blockhead:
        bdic,data = get_nblocks(f,dic,dic["nblocks"],read_blockhead)
        dic["blockheader"] = bdic
    else:
        data = get_nblocks(f,dic,dic["nblocks"],read_blockhead)

    f.close()
    data = uninterleave_data(data)

    return dic,np.squeeze(data)

def read_fid_lowmem(filename,(lenZ,lenY)=(None,None)):
    """
    Read a varian binary (fid) file with mimimal memory usage returning a 
    dic,data pair

    data will be an fid_2d of fid_3d object which can be sliced and transposed
    like a numpy array.

    For 3D data sets the length of the two indirect dimensions must be provided
    as lenZ and lenY.  If not provided the data set is open as if it was a 2D
    data set.

    """

    # peak at fileheader to determind dimensionality
    f = open(filename)
    dic = fileheader2dic(get_fileheader(f))
    f.close()

    if dic["nblocks"] == 1: # 1D data, no lowmem function
        return read_fid(filename)
    
    if lenZ == None or lenY == None:    # 2D data
        return read_fid_lowmem_2D(filename)    
    else:   # 3D data
        return read_fid_lowmem_3D(filename)

def read_fid_lowmem_2D(filename):
    """ 
    Read a varian binary (fid) file as a 2D array with mimimal memory usage.
    """
    
    # create the fid_3d and data
    data = fid_2d(filename)
    dic = dict(data.dic)

    return dic,data


def read_fid_lowmem_3D(filename,(lenZ,lenY)):
    """ 
    Read a varian fid file as a 3D array with mimimal memory usage.
    """
    
    # create the fid_3d and data
    data = fid_3d(filename,(lenZ,lenY))
    dic = dict(data.dic)

    return dic,data


# fid_2d and fid_3d classes

class fid_2d(fileiobase.data_2d):
    """
    fid_2d emulates numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes functions create a new fid_2d object with the
      new axes ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order = ["y","x"]):
        """ set up a fid_2d object"""

        # open the file
        self.filename = filename
        self.f = open(filename)
        
        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(self.f))
        
        # sizes
        self.lenX = self.dic["np"]/2
        self.lenY = self.dic["nblocks"]

        # order  
        self.order = order

        # shape based on order
        a = [self.lenY,self.lenX]
        self.shape = ( a[order.index("y")], a[order.index("x")] )
        del(a)

        # dtype
        self.dtype = find_cdtype(self.dic)
        self.ndim = 2

    def __fcopy__(self,order):
        """ create a copy"""

        n = fid_2d(self.filename,order)
        return n
    

    def __fgetitem__(self,(sY,sX)):
        """ 
        Returns ndarray of selected values

        (sY,sX) are a well formatted 2-tuple of slices

        """
        # make the empty data array
        lensY = len(range(self.lenY)[sY])
        lensX = len(range(self.lenX)[sX])

        dt = find_cdtype(self.dic)
        data = np.empty( (lensY,lensX) , dtype=dt)

        # read in the data
        for jY,iY in enumerate(range(self.lenY)[sY]):
                
            ntrace = iY
            self.f.seek(ntrace*self.dic["bbytes"]+32)
            trace = get_block(self.f,self.dic,False)
            trace = uninterleave_data(trace)[sX]

            data[jY] = trace

        return data


class fid_3d(fileiobase.data_3d):
    """
    fid_3d emulates numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes functions create a new fid_3d object with the
      new axes ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,(lenZ,lenY),order = ["z","y","x"]):
        """set up a fid_3d fobject"""

        # open the file
        self.filename = filename
        self.f = open(filename)
        
        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(self.f))
        
        if self.dic['nblocks'] != lenZ*lenY:
            raise ValueError,"lenZ*lenY does not match size %s" % \
            self.dic['nblocks']

        # sizes
        self.lenX = self.dic["np"]/2
        self.lenY = lenY
        self.lenZ = lenZ

        # order  
        self.order = order

        # shape based on order

        a = [self.lenZ,self.lenY,self.lenX]
        self.shape = (a[order.index("z")],a[order.index("y")],
                      a[order.index("x")] )
        
        # dtype
        self.dtype = find_cdtype(self.dic)
        self.ndim = 3

        del(a)

    def __fcopy__(self,order):
        """ create a copy"""

        n = fid_3d(self.filename,(self.lenZ,self.lenY),order)
        return n
    

    def __fgetitem__(self,(sZ,sY,sX)):
        """ 
        Return ndarray of selected values
            
        (sZ,sY,sX) is a well formateed 3-tuple of slices

        """
    
        # make the empty data array
        lensZ = len(range(self.lenZ)[sZ])
        lensY = len(range(self.lenY)[sY])
        lensX = len(range(self.lenX)[sX])

        dt = find_cdtype(self.dic)
        data = np.empty( (lensZ,lensY,lensX) , dtype=dt)

        # read in the data
        for jZ,iZ in enumerate(range(self.lenZ)[sZ]):
            for jY,iY in enumerate(range(self.lenY)[sY]):
               
                ntrace = int(np.floor(iZ/2.)*self.lenY*2+iY*2+(iZ%2))
                #ntrace = iZ*self.lenY+iY   # linear blocks (not common)
                self.f.seek(ntrace*self.dic["bbytes"]+32)
                trace = get_block(self.f,self.dic,False)
                trace = uninterleave_data(trace)[sX]

                data[jZ,jY] = trace

        return data

# data type functions

def find_cdtype(dic):
    """ 
    Find the data complex dtype from a dictionary
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
    Find the data dtype from a dictionary 
    """

    if dic["S_FLOAT"] == 1:
        return np.dtype('>f4') # float32
    else:
        if dic["S_32"] == 1:
            return np.dtype('>i4') # int32
        else:
            return np.dtype('>i2') # int16


# block get/put functions

def get_block(file,filedic,read_blockhead=False):
    """ 
    Read a block from file described by filedic dictionary

    When read_blockhead is True will read blockheader(s) and returns 
    block_dic,trace.  When read_blockhead is False returns trace.

    """
    
    # find the dtype
    dt = find_dtype(filedic)

    # Do not return blockheaders
    if read_blockhead == False: # Do not return blockheaders
        for i in xrange(filedic["nbheaders"]):
            skip_blockheader(file)
        trace = get_trace(file,filedic["np"],dt)
        return trace
    
    # read the blockheaders
    else:
        dic = dict()
        # read the headers
        if filedic["nbheaders"] >= 1:
            dic.update(blockheader2dic(get_blockheader(file)))
        if filedic["nbheaders"] >= 2:
            dic["hyperhead"] = hyperheader2dic(get_hyperheader(file))
        if filedic["nbheaders"] >= 3:
            for i in xrange(2,filedic["nbheaders"]):
                skip_blockheader(file)
        # read the data
        trace = get_trace(file,filedic["np"],dt)
        
        return dic,trace

def get_nblocks(file,filedic,n,read_blockhead=False):
    """ Read n blocks from file described by filedic dictionary

    When read_blockhead is True will read blockheader(s) and returns
    block_dic,data.  When read_blockheader is False returns data.

    """

    # create an empty array to hold data
    dt = find_dtype(filedic)
    data = np.empty( (filedic["nblocks"],filedic["np"]) ,dtype=dt)

    if read_blockhead:
        bdic = [0]*filedic["nblocks"]

    # read the data
    for i in xrange(filedic["nblocks"]):
        if read_blockhead:
            bdic[i],data[i] = get_block(file,filedic,read_blockhead)
        else:
            data[i] = get_block(file,filedic,read_blockhead)

    if read_blockhead:
        return bdic,data
    else:
        return data

def put_block(file,dic,trace,bh,hh=False):
    """ 
    Put blockheader(s) and trace to file

    Parameters:

    * file  file object
    * dic   fileheader dictionary
    * trace trace to write to block
    * bh    blockheader list
    * hh    hyperheader list (if required)

    When dic["nblocks"] > 2 additional blocks are written as all zeros.

    """
    # write the block headers
    if dic["nbheaders"] >= 1:
        put_blockheader(file,bh)

    if dic["nbheaders"] >= 2:
        if hh == False:
            raise Exception,"Hyperheader required"
        put_hyperheader(file,hh)
    
    # write any additional blockheaders as 0s if needed
    for i in xrange(dic["nbheaders"]-2):
        put_blockheader(file,[0]*9)

    # write the trace
    put_trace(file,trace)

    return

def put_nblocks(file,dic,n,data,bh,hh=False):
    """ 
    Put nblocks to file

    Parameters:

    * file  file object.
    * dic   fileheader dictionary.
    * n     number of blocks to write.
    * data  sliceable object which returns traces.
    * bh    first blockheader list.
    * hh    hyperheader list (if required).

    The index in the blockheader will be automatically incremeted.

    """
    n = int(n)
    
    # check that data is large enough
    if len(data) < n:
        raise ValueError,"data must be at lease of length %i"%(n)

    # store the first block header index
    bh_i = int(bh[2])

    for i in xrange(n): 
        bh[2] = bh_i+i  # increment the block header index
        put_block(file,dic,data[i],bh,hh)

    return


# header get/put functions

def get_fileheader(file):
	"""
    Unpack file header parameters into list.

    Reads the 32-byte file header from file and unpacks into a list.  Endiness
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
	return struct.unpack('>6lhhl',file.read(32))

def put_fileheader(file,fh):
    """ 
    Write fileheader list to file (32-bytes written)

    Parameters:

    * file  file object
    * fh    fileheader list (length 9)

    """

    file.write( struct.pack('>6lhhl',*fh) )
    return

def get_blockheader(file):
	"""
    Unpack block header parameters into a list.
	
    Reads the 28-byte block header from file and unpacks into a list.  Endiness
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
	return struct.unpack('>4hl4f',file.read(28))

def skip_blockheader(file):
	""" 
    Read a block header but do not unpack

    This is a replacement for get_blockheader.  It skips file ahead 28 bytes.
	"""
	file.read(28)
	return

def put_blockheader(file,bh):
    """ 
    Write blockheader list to file (28-bytes written)

    Parameters:

    * file  file object
    * bh    blockheaders list, length 9

    """
    file.write( struct.pack('>4hl4f',*bh) )
    return

def get_hyperheader(file):
	""" 
    Unpack hypercomplex header parameters to a list.

    Reads the 28-bytes block header from file and unpacks into a list. Endiness
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
	return struct.unpack('>4hl4f',file.read(28))


def put_hyperheader(file,hh):
    """ Write hyperheader list to file (28-bytes written)
    
    Parameters:

    * file  file object
    * hh    hyperheader list, length 9

    """
    file.write( struct.pack('>4hl4f',*bh) )
    return


def get_trace(file,num_points,dt):
    """ Read trace of num_points of dtype dt from file

    Does not correct Endiness as dt should handle this

    dt should be:

    * np.dtype('>f4') when S_FLOAT == 1
    * np.dtype('>i2') when S_32 != 1 and S_FLOAT == 0
    * np.dtype('>i4') when S_32 == 1 and S_FLOAT == 0

    """

    bsize = num_points*dt.itemsize  # number of bytes in trace
    return np.frombuffer(file.read(bsize),dt)


# trace put functions

def put_trace(file,trace):
    """ 
    Write trace to file 
    """
    
    file.write(trace.tostring())
    return


# data conversion functions

def uninterleave_data(data_in):
    """ 
    Unpack interleaved real,imag data
    
    ==================  ============
    data_in dtype.name  Return dtype
    ==================  ============
    int16               'complex64'
    float32             'complex64'
    int32               'complex128'
    ==================  ============
    """
    # determind the output dtype
    rdt = data_in.dtype.name

    if rdt == 'int16' or rdt == "float32":
        cdt = "complex64"
    elif rdt == 'int32':
        cdt = "complex128"
    else:
        raise ValueError,"unknown dtype"

    return data_in[...,::2]+np.array(data_in[...,1::2]*1.j,dtype=cdt)

def interleave_data(data_in):
    """ 
    Interleave real, imag data

    Does not check if resulting dtype is a valid varian dtype

    """

    size = list(data_in.shape)
    size[-1] = size[-1]*2
    data_out = np.empty(size,dtype=data_in.real.dtype)

    data_out[...,::2] = data_in.real
    data_out[...,1::2] = data_in.imag

    return data_out

def sign_adj_2Ddata(data):
    """ 
    Sign adjust 2D data

    Use this function to match NMRPipe data generated from var2pipe
    
    Negate imag part of even fids, real part of odd fids 
    """

    if data.ndim == 1:
        data.imag = -data.imag
    else:
        data.imag[::2] = -data.imag[::2]
        data.real[1::2] = -data.real[1::2]

    return data


# Header manipulation functions

def hyperheader2dic(head):
	""" 
    Convert a hypercomplex block header into a python dictionary.
	"""
	dic = dict()
	
	dic["s_spare1"]	= head[0]
	dic["status"]	= head[1]
	dic["s_spare2"]	= head[2]
	dic["s_spare3"]	= head[3]
	dic["l_spare1"]	= head[4]
	dic["lpval1"]	= head[5]
	dic["rpval1"]	= head[6]
	dic["f_spare1"]	= head[7]
	dic["f_spare2"]	= head[8]

	#unpack the status bits
	dic["UHYPERCOMPLEX"] = (dic["status"] & 0x2) / 0x2

	return dic

def repack_hyperheader(dic):
    """ 
    Repack hyperheader dictionary bit flag parameters into status.
    """
    dic["status"] = dic["UHYPERCOMPLEX"]*0x2
    return dic


def dic2hyperheader(head):
	"""	
    Convert a python dictionary into a hypercomplex block header list.

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


def make_blockheader(filedic=False,index=1):
    """ 
    Make a generic blockheader dictionary.

    filedic can be provided for status flags, if not provided creates 
    header for float32 data

    """

    dic = dict()

    if filedic != False:
        # common status flags
        dic["S_DATA"]  = filedic["S_DATA"] 
        dic["S_SPEC"]  = filedic["S_SPEC"]
        dic["S_32"]    = filedic["S_32"]
        dic["S_FLOAT"] = filedic["S_FLOAT"]
        dic["S_COMPLEX"] = filedic["S_COMPLEX"]
        dic["S_HYPERCOMPLEX"] = filedic["S_HYPERCOMPLEX"]
    else:
        # common status flags
        dic["S_DATA"]  = 1
        dic["S_SPEC"]  = 0
        dic["S_32"]    = 0
        dic["S_FLOAT"] = 1
        dic["S_COMPLEX"] = 0 
        dic["S_HYPERCOMPLEX"] = 0
    
    # blockheader specific status flags
    dic["MORE_BLOCKS"] = 1
    dic["NP_CMPLX"]    = 0
    dic["NF_CMPLX"]    = 0
    dic["NI_CMPLX"]    = 0
    dic["NI2_CMPLX"]   = 0

    # mode flags
    dic["NP_PHMODE"]   = 0
    dic["NP_AVMODE"]   = 0
    dic["NP_PWRMODE"]  = 0
    dic["NF_PHMODE"]   = 0
    dic["NF_AVMODE"]   = 0
    dic["NF_PWRMODE"]  = 0
    dic["NI_PHMODE"]   = 0
    dic["NI_AVMODE"]   = 0
    dic["NI_PWRMODE"]  = 0
    dic["NI2_PHMODE"]  = 0
    dic["NI2_AVMODE"]  = 0
    dic["NI2_PWRMODE"] = 0

    # main parameters
    dic["scale"]   = 0
    dic["status"]  = 0
    dic["index"]   = index
    dic["mode"]    = 0
    dic["ctcount"] = 1
    dic["lpval"]   = 0.0
    dic["rpval"]   = 0.0
    dic["lvl"]     = 0.0
    dic["tlt"]     = 0.0

    repack_blockheader(dic) # repack the header

    return dic

def blockheader2dic(head):
    """ 
    Convert a block header list into a python dictionary.
    """

    dic = dict()

    dic["scale"] 	= head[0]
    dic["status"]	= head[1]
    dic["index"]	= head[2]
    dic["mode"]	    = head[3]
    dic["ctcount"]	= head[4]
    dic["lpval"]	= head[5]
    dic["rpval"]	= head[6]
    dic["lvl"]	    = head[7]
    dic["tlt"]	    = head[8]
    
    # unpack the status parameters  
    dic["S_DATA"]    = (dic["status"] & 0x1)/0x1
    dic["S_SPEC"]    = (dic["status"] & 0x2)/0x2
    dic["S_32"]      = (dic["status"] & 0x4)/0x4
    dic["S_FLOAT"]   = (dic["status"] & 0x8)/0x8
    dic["S_COMPLEX"] = (dic["status"] & 0x10)/0x10
    dic["S_HYPERCOMPLEX"]   = (dic["status"] & 0x20)/0x20

    dic["MORE_BLOCKS"]	= (dic["status"] & 0x80)/0x80
    dic["NP_CMPLX"]	 = (dic["status"] & 0x100)/0x100
    dic["NF_CMPLX"]	 = (dic["status"] & 0x200)/0x200
    dic["NI_CMPLX"]	 = (dic["status"] & 0x400)/0x400
    dic["NI2_CMPLX"] = (dic["status"] & 0x800)/0x800

    # unpack the mode parameter
    dic["NP_PHMODE"]	= (dic["mode"] & 0x1) / 0x1
    dic["NP_AVMODE"]	= (dic["mode"] & 0x2) / 0x2
    dic["NP_PWRMODE"]	= (dic["mode"] & 0x4) / 0x4

    dic["NF_PHMODE"]	= (dic["mode"] & 0x10) / 0x10
    dic["NF_AVMODE"]	= (dic["mode"] & 0x20) / 0x20
    dic["NF_PWRMODE"]	= (dic["mode"] & 0x40) / 0x40
    	
    dic["NI_PHMODE"]	= (dic["mode"] & 0x100) / 0x100
    dic["NI_AVMODE"]	= (dic["mode"] & 0x200) / 0x200
    dic["NI_PWRMODE"]	= (dic["mode"] & 0x400) / 0x400 

    dic["NI2_PHMODE"]	= (dic["mode"] & 0x1000) / 0x1000
    dic["NI2_AVMODE"]	= (dic["mode"] & 0x2000) / 0x2000
    dic["NI2_PWRMODE"]	= (dic["mode"] & 0x4000) / 0x4000

    return dic


def repack_blockheader(dic):
    """ 
    Repack blockheader dic bit flag parameters into status and mode.
    """

    dic["status"] = dic["S_DATA"]*0x1 + dic["S_SPEC"]*0x2 + dic["S_32"]*0x4 + \
                    dic["S_FLOAT"]*0x8 + dic["S_COMPLEX"]*0x10 +              \
                    dic["S_HYPERCOMPLEX"]*0x20 + dic["MORE_BLOCKS"]*0x80 +    \
                    dic["NP_CMPLX"]*0x100 + dic["NF_CMPLX"]*0x200 +           \
                    dic["NI_CMPLX"]*0x400 + dic["NI2_CMPLX"]*0x800         

    dic["mode"] = dic["NP_PHMODE"]*0x1 + dic["NP_AVMODE"]*0x2 +          \
                  dic["NP_PWRMODE"]*0x4 + dic["NF_PHMODE"]*0x10 +        \
                  dic["NF_AVMODE"]*0x20 + dic["NF_PWRMODE"]*0x40 +       \
                  dic["NI_PHMODE"]*0x100 + dic["NI_AVMODE"]*0x200 +      \
                  dic["NI_PWRMODE"]*0x400 + dic["NI2_PHMODE"]*0x1000 +   \
                  dic["NI2_AVMODE"]*0x2000 + dic["NI2_PWRMODE"]*0x4000
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
    Convert fileheader list into a python dictionary
	"""
	dic = dict()

	dic["nblocks"]	= head[0]
	dic["ntraces"]	= head[1]
	dic["np"]	= head[2]
	dic["ebytes"]	= head[3]
	dic["tbytes"]	= head[4]
	dic["bbytes"]	= head[5]
	dic["vers_id"]	= head[6]
	dic["status"]	= head[7]
	dic["nbheaders"]= head[8]

	# unpack the status parameter
	dic["S_DATA"]	= (dic["status"] & 0x1) / 0x1
	dic["S_SPEC"]	= (dic["status"] & 0x2) / 0x2
	dic["S_32"]	    = (dic["status"] & 0x4) / 0x4
	dic["S_FLOAT"]	= (dic["status"] & 0x8) / 0x8
	dic["S_COMPLEX"]= (dic["status"] & 0x10) / 0x10
	dic["S_HYPERCOMPLEX"]	= (dic["status"] & 0x20) / 0x20
	dic["S_ACQPAR"]	= (dic["status"] & 0x80) / 0x80
	dic["S_SECND"]	= (dic["status"] & 0x100) / 0x100
	dic["S_TRANSF"]	= (dic["status"] & 0x200) / 0x200
	dic["S_NP"]	    = (dic["status"] & 0x800) / 0x800
	dic["S_NF"]	    = (dic["status"] & 0x1000) / 0x1000
	dic["S_NI"]	    = (dic["status"] & 0x2000) / 0x2000
	dic["S_NI2"]	= (dic["status"] & 0x4000) / 0x4000
	
	return dic


def repack_fileheader(dic):
    """ 
    Repack blockheader dic bit flag parameters into status and mode.
    """

    dic["status"] = dic["S_DATA"]*0x1 + dic["S_SPEC"]*0x2 + dic["S_32"]*0x4 + \
                    dic["S_FLOAT"]*0x8 + dic["S_COMPLEX"]*0x10 +             \
                    dic["S_HYPERCOMPLEX"]*0x20 + dic["S_ACQPAR"]*0x80 +      \
                    dic["S_SECND"]*0x100 + dic["S_TRANSF"]*0x200 +           \
                    dic["S_NP"]*0x800 + dic["S_NF"]*0x1000 +                 \
                    dic["S_NI"]*0x2000 + dic["S_NI2"]*0x4000

    return dic

def dic2fileheader(dic):
	""" 
    Convert a python dictionary into a fileheader list

    Does not repack status from bit flags
	"""

	head = [0] * 9
 	head[0]= dic["nblocks"]
  	head[1]= dic["ntraces"]
  	head[2]= dic["np"]
  	head[3]= dic["ebytes"]
  	head[4]= dic["tbytes"]
  	head[5]= dic["bbytes"]
  	head[6]= dic["vers_id"]
  	head[7]= dic["status"]
  	head[8]= dic["nbheaders"]

	return head


# procpar functions

def read_procpar(filename):
	""" 
    Read a procpar file returning a python dictionary
	"""

	f = open(filename)

	dic = dict()

	length = os.stat(filename).st_size

	# test to see if end of file
	while f.tell() != length:
		p = get_parameter(f)
		dic[p["name"]] = p

	return dic


def write_procpar(filename,dic,overwrite=False):
    """ 
    Write a varian procpar file from python dictionary
    """

    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite=overwrite)

    for key in dic.keys():	# loop over the parameters

        d = dic[key]
        # print out the first line
        print >> f,d["name"],d["subtype"],d["basictype"],     \
                   d["maxvalue"],d["minvalue"],d["stepsize"], \
                   d["Ggroup"],d["Dgroup"],d["protection"],   \
                   d["active"],d["intptr"]

        # print out the next line(s) (and more if strings)
        if d["basictype"] == "1":	# real values, one line

            print >> f,len(d["values"]), # don't end the line
            for value in d["values"]:
                print >>f,value,	# still not yet
            print >> f,""	# now end the line

        elif d["basictype"] == "2":	# strings, may have multiple lines

            print >> f,len(d["values"]),	# don't end the line
            for value in d["values"]:
                print >> f,'"'+value+'"'	# now end the line (for each string)

        # print out the last line
        print >> f,d["enumerable"],

        if d["enumerable"] != "0":
            for e in d["enumerables"]:
                if d["basictype"] == "1": #reals
                    print >> f,e,

                elif d["basictype"] == "2": #strings
                    print >> f,'"'+e+'"',

        print >> f,""	# end the enumerable line	

    f.close()

    return

subtypes = ["undefined", "real", "string","delay","flag",
	"frequency","pulse","integer"]

basictypes = ["undefined","real","string"]

def get_parameter(file):
    """ 
    Reads a procpar parameter from a file object

    Returns a dictionary with the attributes of the parameter
    """

    dic = dict()

    # read and decode the first line
    line = file.readline().split()

    dic["name"] 		= line[0]
    dic["subtype"]		= line[1]
    dic["basictype"]	= line[2]
    dic["maxvalue"]		= line[3]
    dic["minvalue"]		= line[4]
    dic["stepsize"]		= line[5]
    dic["Ggroup"]		= line[6]
    dic["Dgroup"]		= line[7]
    dic["protection"]	= line[8]
    dic["active"]		= line[9]
    dic["intptr"]		= line[10]

    # read in the values of the parameter
    line = file.readline()
    
    num = int(line.split()[0])

    values = []

    if dic["basictype"] == "1":	# real values, only one line
        values = line.split()[1:]

    elif dic["basictype"] == "2":	# strings, may have multiple lines

        values.append(line.split("\"")[1])	# split on "s

        for i in range(num-1):
            values.append(file.readline().split("\"")[1])

    dic["values"] = values

    line = file.readline()

    # read and decode the enumerables
    dic["enumerable"] = line.split()[0] 

    if dic["enumerable"] != "0":

        if dic["basictype"] == "1":	# reals
            dic["enumerables"] = line.split()[1:]

        elif dic["basictype"] == "2":	#strings		
            dic["enumerables"] = line.split("\"")[1::2]

    return dic
