"""
Functions for reading and writing sparky (.ucsf) files.
"""
__developer_info__ = """
Information on the sparky file format can be found online at:
http://www.cgl.ucsf.edu/home/sparky/manual/files.html
and in the source file ucsffile.cc.

"""

import numpy as np
import struct
import os
import fileiobase
import datetime

# unit conversion function

def make_uc(dic,data,dim=-1):
    """ 
    Make a unit conversion object 
 
    Parameters:

    * dic   Sparky dictionary
    * data  data array
    * dim   dimention to make converter for (0,1,2,3 or -1, last) 
    """

    if dim == -1:
        dim = data.ndim - 1 # last dimention 
   
    wdic = dic["w"+str(int(1+dim))]

    size = float(wdic["npoints"])
    cplx = False
    sw   = wdic["spectral_width"]
    obs  = wdic["spectrometer_freq"]
    car  = wdic["xmtr_freq"]*obs        

    return fileiobase.unit_conversion(size,cplx,sw,obs,car)


# dictionary/data creation


def create_data(data):
    """ 
    Create a sparky data array (recast into float32 array)
    
    """
    return np.array(data,dtype="float32")


def guess_udic(dic,data):
    """ 
    Guess parameter of universal dictionary from dic,data pair
    """

    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in xrange(data.ndim):

        adic = dic["w"+str(i+1)]
        udic[i]["size"] = data.shape[i]
        udic[i]["sw"] = adic['spectral_width']
        udic[i]["obs"] = adic['spectrometer_freq']
        udic[i]["car"] = adic['xmtr_freq']*adic['spectrometer_freq']
        udic[i]["label"] = adic['nucleus']

        udic[i]["complex"] = False
        udic[i]["time"] = False
        udic[i]["freq"] = True

    return udic

def create_dic(udic,datetimeobj=datetime.datetime.now(),user='user'):
    """ 
    Create a sparky dictionary from universal dictionary
    """

    dic = dict()
   
    # determind shape of array
    shape = [udic[k]["size"] for k in xrange(udic["ndim"])]

    # populate the dictionary
    dic["ident"] = 'UCSF NMR'
    dic["naxis"] = udic["ndim"]
    dic["ncomponents"] = 1
    dic["encoding"] = 0
    dic["version"] = 2
    dic["owner"] = user
    dic["date"] =  datetimeobj.ctime()
    dic["comment"] = ''
    dic["scratch"] = ''

    # calc a good tile shape
    tshape = calc_tshape(shape)

    # total number of tiles
    ntiles = 1
    for tlen,slen in zip(tshape,shape):
        ntiles*=np.ceil(float(slen)/tlen)
   
    # points in tile
    tpoints = np.array(tshape).prod()

    # data bytes
    dbytes = tpoints*ntiles*4

    # total file size if data size plus leaders
    dic["seek_pos"] = int(dbytes+180+128*len(shape))


    # populate the dictionary with axis dictionaries
    for i,(tlen,dlen) in enumerate(zip(tshape,shape)):
        
        dic["w"+str(i+1)] = create_axisdic(udic[i],tlen,dlen)

    return dic


def create_axisdic(adic,tlen,dlen):
    """ 
    Make an sparky axis dictionary from a universal axis dictionary 

    Parameters:

    * adic  axis dictionary from universal dictionary
    * tlen  tile length
    * dlen  data length

    """

    dic = dict()

    dic["nucleus"] = adic["label"]
    dic["spectral_shift"] = 0
    dic["npoints"] = int(dlen)
    dic["size"] = int(dlen)
    dic["bsize"] = int(tlen)
    dic["spectrometer_freq"] = float(adic["obs"])
    dic["spectral_width"] = float(adic["sw"])
    dic["xmtr_freq"] = float(adic["car"])/dic["spectrometer_freq"]
    dic["zero_order"] = 0.0
    dic["first_order"] = 0.0
    dic["first_pt_scale"] = 0.0
    dic["extended"] = '\x80'        # transform bit set

    return dic


def datetime2dic(datetimeobj,dic):
    """ 
    Add time datetime object to dictionary 
    """
    dic["date"] = datetimeobj.ctime()

    return dic


def dic2datetime(dic):
    """ 
    Create a datetime object from sparky dictionary
    """
    return datetime.datetime.strptime(dic["date"],"%a %b %d %H:%M:%S %Y")


def calc_tshape(shape,kbyte_max=128):
    """ 
    Calculate tile shape from data shape
    
    shape is a tuple representing the data shape (data.shape)
    kbyte_max determined the largest tile size in Kbytes

    Algorithm divides each dimention by 2 until under kbyte_max tile size.

    """

    s = np.array(shape,dtype="int")
    i = 0
    while (s.prod()*4./1024. > kbyte_max):
        s[i] = np.floor(s[i]/2.)
        i = i+1
        if i == len(s): i = 0
    return tuple(s)


# global read/write functions

def read(filename):
    """ 
    Read a sparky file returning a dic,data pair
    """
    # open the file
    f = open(filename)

    # determind the dimentionality
    n = fileheader2dic(get_fileheader(f))["naxis"]
    f.close()

    if n == 2:
        return read_2D(filename)
    if n == 3:
        return read_3D(filename)
    
    raise ValueError,"unknown dimentionality: %s"%n


def read_lowmem(filename):
    """ 
    Read a sparky file with minimal memory returning a dic,data pair.
    """

    # open the file
    f = open(filename)

    # determind the dimentionality
    n = fileheader2dic(get_fileheader(f))["naxis"]
    f.close()

    if n == 2:
        return read_lowmem_2D(filename)
    if n == 3:
        return read_lowmem_3D(filename)
    
    raise ValueError,"unknown dimentionality: %s"%order


def write(filename,dic,data,overwrite=False):
    """ 
    Write a sparky file

    Parameters:

    * filename  Name of file to write to.
    * data      Data array.
    * dic       Sparky parameter dictionary.
    * overwrite Set to True to overwrite existing file.

    No return

    """

    n = dic["naxis"]

    if n == 2:
        return write_2D(filename,dic,data,overwrite=overwrite)
    if n == 3:
        return write_3D(filename,dic,data,overwrite=overwrite)
    
    raise ValueError,"unknown dimentionality: %s"%order

def write_lowmem(filename,dic,data,overwrite=False):
    """
    Write a sparky file tile by tile (low memory)

    Parameters:

    * filename  Name of file to write to.
    * data      Data array.
    * dic       Sparky parameter dictionary.
    * overwrite Set to True to overwrite existing file.
    
    No return

    """
    # write also writes tile by tile...
    return write(filename,dic,data,overwrite)

# dimensional reading/writing functions

def read_2D(filename):
    """ 
    Read a 2D sparky file returning a dic,data pair.
    """

    seek_pos = os.stat(filename).st_size
    f = open(filename)

    # read the file header
    dic = fileheader2dic(get_fileheader(f))

    # check for file size mismatch
    if seek_pos != dic["seek_pos"]:
        raise IOError,"Bad file size %s vs %s",(seek_pos,dic["seek_pos"])

    # read the axis headers...
    for i in xrange(dic['naxis']):
        dic["w"+str(i+1)] = axisheader2dic(get_axisheader(f))

    # read the data and untile
    lenY = dic["w1"]["npoints"]  
    lenX = dic["w2"]["npoints"]
    lentY = dic["w1"]["bsize"]
    lentX = dic["w2"]["bsize"]
    data = get_data(f)
    data = untile_data2D(data,(lentY,lentX),(lenY,lenX))

    return dic,data


def write_2D(filename,dic,data,overwrite=False):
    """ 
    Write a sparky file from 2D data
    """

    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite)

    # write the file header
    put_fileheader(f,dic2fileheader(dic))

    # write the axis headers
    put_axisheader(f,dic2axisheader(dic["w1"]))
    put_axisheader(f,dic2axisheader(dic["w2"]))

    lentX = dic["w2"]["bsize"]
    lentY = dic["w1"]["bsize"]
    t_tup = (lentY,lentX)

    ttX = np.ceil(data.shape[1] / float(lentX)) # total tiles in X dim
    ttY = np.ceil(data.shape[0] / float(lentY)) # total tiles in Y dim
    
    tt = ttX*ttY

    for i in xrange(int(tt)):
        put_data( f,find_tilen_2d( data,i,(t_tup) ) )

    f.close()
    return


def read_3D(filename):
    """ 
    Read a 3D sparky file returning a dic,data pair
    """

    seek_pos = os.stat(filename).st_size
    f = open(filename)

    # read the file header
    dic = fileheader2dic(get_fileheader(f))

    # check for file size mismatch
    if seek_pos != dic["seek_pos"]:
        raise IOError,"Bad file size %s vs %s",(seek_pos,dic["seek_pos"])

    # read the axis headers...
    for i in xrange(dic['naxis']):
        dic["w"+str(i+1)] = axisheader2dic(get_axisheader(f))

    # read the data and untile
    lenZ = dic["w1"]["npoints"]  
    lenY = dic["w2"]["npoints"]
    lenX = dic["w3"]["npoints"]
    lentZ = dic["w1"]["bsize"]
    lentY = dic["w2"]["bsize"]
    lentX = dic["w3"]["bsize"]
    data = get_data(f)
    data = untile_data3D(data,(lentZ,lentY,lentX),(lenZ,lenY,lenX))

    return dic,data


def write_3D(filename,dic,data,overwrite=False):
    """ 
    Write a sparky file from 3D data
    """

    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite)

    # write the file header
    put_fileheader(f,dic2fileheader(dic))

    # write the axis headers
    put_axisheader(f,dic2axisheader(dic["w1"]))
    put_axisheader(f,dic2axisheader(dic["w2"]))
    put_axisheader(f,dic2axisheader(dic["w3"]))

    lentX = dic["w3"]["bsize"]
    lentY = dic["w2"]["bsize"]
    lentZ = dic["w1"]["bsize"]

    t_tup = (lentZ,lentY,lentX)

    ttX = np.ceil(data.shape[2] / float(lentX)) # total tiles in X dim
    ttY = np.ceil(data.shape[1] / float(lentY)) # total tiles in Y dim
    ttZ = np.ceil(data.shape[0] / float(lentZ)) # total tiles in Z dim
    
    tt = ttX*ttY*ttZ

    for i in xrange(int(tt)):
        put_data( f,find_tilen_3d( data,i,(t_tup) ) )

    f.close()
    return


# read_lowmem functions


def read_lowmem_2D(filename):
    """ 
    Read a 2D sparky file with minimal memory usage.
    """
    seek_pos = os.stat(filename).st_size

    # create the sparky_2d file
    data = sparky_2d(filename)
    dic  = dict(data.dic)

    # check for file size mismatch
    if seek_pos != dic["seek_pos"]:
        raise IOError,"Bad file size %s vs %s",(seek_pos,dic["seek_pos"])
    return dic,data


def read_lowmem_3D(filename):
    """ 
    Read a 3D sparky file with minimal memory usage
    """
    seek_pos = os.stat(filename).st_size

    # create the sparky_3d file
    data = sparky_3d(filename)
    dic  = dict(data.dic)

    # check for file size mismatch
    if seek_pos != dic["seek_pos"]:
        raise IOError,"Bad file size %s vs %s",(seek_pos,dic["seek_pos"])

    return dic,data


# sparky_* objects

class sparky_2d(fileiobase.data_nd):
    """
    Emulates a numpy.ndarray object without loading data into memory for low
    memory reading of 2D Sparky files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order=None):
        """
        Create and set up object
        """
        # open the file
        self.filename = filename
        f = open(filename)

        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(f))

        if self.dic["naxis"] != 2:
            raise StandardError,"file is not a 2D Sparky file"

        # read in the axisheaders
        self.dic["w1"] = axisheader2dic(get_axisheader(f))
        self.dic["w2"] = axisheader2dic(get_axisheader(f))
        f.close()

        # sizes
        self.lenY = self.dic["w1"]["npoints"]
        self.lenX = self.dic["w2"]["npoints"]

        # tile sizes
        self.lentY = self.dic["w1"]["bsize"]
        self.lentX = self.dic["w2"]["bsize"]

        # check order
        if order == None:
            order = (0, 1)
        
        # finalize
        self.dtype = np.dtype("float32")
        self.order = order
        self.fshape = (self.lenY, self.lenX)
        self.__setdimandshape__()

    def __fcopy__(self,order):
        
        n = sparky_2d(self.filename,order)
        return n

    def __fgetitem__(self,slices):
        """
        Returns ndarray of selected values

        slices is a well formatted 2-tuple of slices
        """
        sY,sX = slices
        
        f = open(self.filename)

        #print sY,sX
        gY = range(self.lenY)[sY]   # list of values to take in Y
        gX = range(self.lenX)[sX]   # list of values to take in X

        # tiles to get in each dim to read
        gtY = set([np.floor(i/self.lentY) for i in gY]) # Y tile to read
        gtX = set([np.floor(i/self.lentX) for i in gX]) # X tile to read

        # create a empty output directory
        out = np.empty( (len(gY),len(gX) ),dtype=self.dtype)
        
        for iY in gtY:      # loop over Y tiles to get
            for iX in gtX:  # loop over X tiles to get
                
                # get the tile and reshape it
                ntile = iY*np.ceil(self.lenX/self.lentX)+iX
                tile = get_tilen(f,ntile,(self.lentX,self.lentY))
                tile = tile.reshape(self.lentY,self.lentX)

                # tile minimum and max values for each dim
                minX = iX*self.lentX
                maxX = (iX+1)*self.lentX

                minY = iY*self.lentY
                maxY = (iY+1)*self.lentY

                # determind what elements are needed from this tile
                XinX = [i for i in gX if maxX>i>=minX ] # values in gX
                XinT = [i-minX for i in XinX]           # tile index values
                XinO = [gX.index(i) for i in XinX]      # output indexes

                YinY = [i for i in gY if maxY>i>=minY ] # values in gX
                YinT = [i-minY for i in YinY]           # tile index values
                YinO = [gY.index(i) for i in YinY]      # output indexes

                # take elements from the tile
                ctile = tile.take(XinT,axis=1).take(YinT,axis=0)

                # DEBUGGING info
                #print "-------------------------------"
                #print "iX:",iX,"iY:",iY,"ntile:",ntile
                #print "tile.shape",tile.shape
                #print "minX:",minX,"maxX",maxX
                #print "minY:",minY,"maxY",maxY
                #print "XinX",XinX
                #print "XinT",XinT
                #print "XinO",XinO
                #print "YinY",YinY
                #print "YinT",YinT
                #print "YinO",YinO

                # put the cut tile to the out array (uses some fancy indexing)
                out[np.ix_(YinO,XinO)] = ctile

        f.close()
        return out


class sparky_3d(fileiobase.data_nd):
    """
    Emulates a numpy.ndarray object without loading data into memory for low
    memory reading of 3D Sparky files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order=None):

        # open the file
        self.filename = filename
        f = open(filename)

        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(f))
        
        if self.dic["naxis"] != 3:
            raise StandardError,"file not 3D Sparky file"

        # read in the axisheaders
        self.dic["w1"] = axisheader2dic(get_axisheader(f))
        self.dic["w2"] = axisheader2dic(get_axisheader(f))
        self.dic["w3"] = axisheader2dic(get_axisheader(f))
        f.close()

        # sizes
        self.lenZ = self.dic["w1"]["npoints"]
        self.lenY = self.dic["w2"]["npoints"]
        self.lenX = self.dic["w3"]["npoints"]

        # tile sizes
        self.lentZ = self.dic["w1"]["bsize"]
        self.lentY = self.dic["w2"]["bsize"]
        self.lentX = self.dic["w3"]["bsize"]

        # check order
        if order == None:
            order = (0, 1, 2)
        
        # finalize
        self.dtype = np.dtype("float32")
        self.order = order
        self.fshape = (self.lenZ, self.lenY, self.lenX)
        self.__setdimandshape__()

    def __fcopy__(self,order):
        
        n = sparky_3d(self.filename,order)
        return n

    def __fgetitem__(self,slices):
        """ 
        Returns ndarray of selected values

        slices is a well formateed 3-tuple of slices

        """
        sZ,sY,sX = slices
        f = open(self.filename)

        gZ = range(self.lenZ)[sZ]   # list of values to take in Z
        gY = range(self.lenY)[sY]   # list of values to take in Y
        gX = range(self.lenX)[sX]   # list of values to take in X

        # tiles to get in each dim to read
        gtZ = set([np.floor(float(i)/self.lentZ) for i in gZ]) # Z tile to read
        gtY = set([np.floor(float(i)/self.lentY) for i in gY]) # Y tile to read
        gtX = set([np.floor(float(i)/self.lentX) for i in gX]) # X tile to read

        # total tiles in each dim
        ttX = np.ceil(self.lenX/float(self.lentX)) # total tiles in X
        ttY = np.ceil(self.lenY/float(self.lentY)) # total tiles in Y
        ttZ = np.ceil(self.lenZ/float(self.lentZ)) # total tiles in Z

        tile_tup = (self.lentZ,self.lentY,self.lentX)

        # create a empty output array
        out = np.empty( ( len(gZ),len(gY),len(gX) ),dtype=self.dtype)
        
        for iZ in gtZ:          # loop over Z tiles to get
            for iY in gtY:      # loop over Y tiles to get
                for iX in gtX:  # loop over X tiles to get
                
                    # get the tile and reshape it
                    ntile = iZ*ttX*ttY + iY*ttX + iX
                    tile = get_tilen(f,ntile,tile_tup)
                    tile = tile.reshape(tile_tup)

                    # tile minimum and max values for each dim
                    minX = iX*self.lentX
                    maxX = (iX+1)*self.lentX

                    minY = iY*self.lentY
                    maxY = (iY+1)*self.lentY

                    minZ = iZ*self.lentZ
                    maxZ = (iZ+1)*self.lentZ

                    # determind what elements are needed from this tile
                    XinX = [i for i in gX if maxX>i>=minX ] # values in gX
                    XinT = [i-minX for i in XinX]           # tile index values
                    XinO = [gX.index(i) for i in XinX]      # output indexes

                    YinY = [i for i in gY if maxY>i>=minY ] # values in gX
                    YinT = [i-minY for i in YinY]           # tile index values
                    YinO = [gY.index(i) for i in YinY]      # output indexes

                    ZinZ = [i for i in gZ if maxZ>i>=minZ ] # values in gX
                    ZinT = [i-minZ for i in ZinZ]           # tile index values
                    ZinO = [gZ.index(i) for i in ZinZ]      # output indexes

                    # take elements from the tile
                    ctile = tile.take(XinT,axis=2).take(YinT,axis=1)
                    ctile = ctile.take(ZinT,axis=0)

                    # DEBUGGING info
                    #print "-------------------------------"
                    #print "iX:",iX,"iY:",iY,"iZ:",iZ,"ntile:",ntile
                    #print "ttX:",ttX,"ttY:",ttY,"ttZ",ttZ
                    #print "tile.shape",tile.shape
                    #print "minX:",minX,"maxX",maxX
                    #print "minY:",minY,"maxY",maxY
                    #print "minZ:",minZ,"maxZ",maxZ
                    #print "XinX",XinX
                    #print "XinT",XinT
                    #print "XinO",XinO
                    #print "YinY",YinY
                    #print "YinT",YinT
                    #print "YinO",YinO
                    #print "ZinZ",ZinZ
                    #print "ZinT",ZinT
                    #print "ZinO",ZinO

                    # put the cut tile to the out array 
                    out[np.ix_(ZinO,YinO,XinO)] = ctile
        f.close()
        return out


# tile and data get/put functions

def get_tilen(f,n_tile,tw_tuple):
    """ 
    Read in tile n from file object with tile sizes given by tw_tuple
    
    Current file position is loss (store before calling if desired)

    """
    # determind the size of the tile in bytes
    tsize = 4
    for i in tw_tuple:
        tsize = tsize*i

    # seek to the beginning of the tile
    f.seek(int(180+128*len(tw_tuple)+n_tile*tsize))
    
    return np.frombuffer(f.read(tsize),dtype='>f4')


def get_tile(f,num_points):
    """ 
    Read tile data from file object
    """
    bsize = num_points*4        # size in bytes

    return np.frombuffer(f.read(bsize),dtype='>f4')


def put_tile(f,tile):
    """ 
    Put tile data to file
    """
    f.write(tile.astype('>f4').tostring())
    return


def get_data(f):
    """ 
    Read all data from sparky file object
    """
    return np.frombuffer(f.read(),dtype='>f4')


def put_data(f,data):
    """ 
    Put data to file

    Does not untile data, assumes this has been done
    
    """
    f.write(data.astype('>f4').tostring())
    return


# tiling/untileing functions

def find_tilen_2d(data,ntile,(lentY,lentX)):
    """ 
    Return a single tile from untiled data 
    
    Parameters:

    * data            untiled data
    * ntile           Tile number to return
    * (lentY,lentX)   Tuple representing tile size

    Returns 1D numpy array of floats, zero filling as needed.

    """

    ttX = np.ceil(data.shape[1] / float(lentX))    # total tiles in X dim
    ttY = np.ceil(data.shape[0] / float(lentY))    # total tiles in Y dim

    # tile number in each dim
    Xt = ntile % ttX
    Yt = int(np.floor(ntile/ttX))

    # dimention limits
    Xmin = int(Xt*lentX)
    Xmax = int((Xt+1)*lentX)

    Ymin = int(Yt*lentY)
    Ymax = int((Yt+1)*lentY)

    tile = data[Ymin:Ymax,Xmin:Xmax]

    # some edge tiles might need zero filling
    # see if this is the case
    if tile.shape == (lentY,lentX):    # well sized tile
        return tile.flatten()
    else:  
        new_tile = np.zeros( (lentY,lentX),dtype="float32")
        new_tile[:tile.shape[0],:tile.shape[1]] = tile
        return new_tile.flatten()


def tile_data2d(data,(lentY,lentX)):
    """ 
    Tile sparky data into 1D numpy array

    Parameters:

    * data    Two-dimensional data array
    * lentY   Y (w1) dimention tile size
    * lentX   X (w2) dimention tile size

    Returns 1D numpy array of floats

    """
    # determind the number of tiles in data
    ttX = np.ceil(data.shape[1] / float(lentX))    # total tiles in X dim
    ttY = np.ceil(data.shape[0] / float(lentY))    # total tiles in Y dim
    tt = ttX*ttY                            # total number of tiles
    
    # calc some basic parameter
    tsize = lentX*lentY     # number of points in one tile
    t_tup = (lentY,lentX)   # tile size tuple

    # create an empty array to store file data
    out = np.empty( (tt*tsize),dtype="float32")
    
    for i in xrange(int(tt)):
        out[i*tsize:(i+1)*tsize] = find_tilen_2d(data,i,t_tup)

    return out


def untile_data2D(data,(lentY,lentX),(lenY,lenX)):
    """ 
    Reorganize tiled sparky data into 2D data

    Parameters:

    * data    1D numpy array of tile data
    * lentY   size of tile in Y (w1) dim
    * lentX   size of tile in X (w2) dim
    * lenY    size of data in Y dim
    * lenX    size of data in X dim

    Returns 2D np.array of floats with size (lenY,lenX)

    """
    # determind the number of tiles in data
    ttX = np.ceil(lenX / float(lentX))    # total tiles in X dim
    ttY = np.ceil(lenY / float(lentY))    # total tiles in Y dim
    tt = ttX*ttY

    # calc some basic parameter
    tsize = lentX*lentY     # number of points in one tile
    t_tup = (lentY,lentX)   # tile size tuple

    # create an empty array to store file data
    out = np.empty( (ttY*lentY,ttX*lentX),dtype="float32")

    for iY in xrange(int(ttY)):
        for iX in xrange(int(ttX)):
            
            minX = iX*lentX
            maxX = (iX+1)*lentX

            minY = iY*lentY
            maxY = (iY+1)*lentY

            ntile = iY*ttX + iX
            minT = ntile * tsize
            maxT = (ntile+1) * tsize

            #print "ntile",ntile
            #print "minX",minX,"maxX",maxX
            #print "minY",minY,"maxY",maxY
            #print "minT",minT,"maxT",maxT

            #print out[minY:maxY,minX:maxX].shape
            #print data[minT:maxT].reshape(t_tup).shape

            out[minY:maxY,minX:maxX] = data[minT:maxT].reshape(t_tup)

    return out[:lenY,:lenX]
            

def find_tilen_3d(data,ntile,(lentZ,lentY,lentX)):
    """ 
    Return a single tile from untiled data 
    
    Parameters:

    * data    untiled data
    * ntile   Tile number to return
    * lentZ   Tile Z (w1) size
    * lentY   Tile Y (w2) size
    * lentX   Tile X (w3) size

    Returns 1D numpy array of floats, zero filling as needed.

    """

    ttX = np.ceil(data.shape[2] / float(lentX))    # total tiles in X dim
    ttY = np.ceil(data.shape[1] / float(lentY))    # total tiles in Y dim
    ttZ = np.ceil(data.shape[0] / float(lentZ))    # total tiles in Z dim

    # tile number in each dim
    Xt = ntile % ttX
    Yt = int(np.floor(ntile/ttX)) % ttY
    Zt = int( np.floor( ntile / (ttX*ttY) ) )

    # dimention limits
    Xmin = int(Xt*lentX)
    Xmax = int((Xt+1)*lentX)

    Ymin = int(Yt*lentY)
    Ymax = int((Yt+1)*lentY)

    Zmin = int(Zt*lentZ)
    Zmax = int((Zt+1)*lentZ)
    
    tile = data[Zmin:Zmax,Ymin:Ymax,Xmin:Xmax]

    # some edge tiles might need zero filling
    # see if this is the case
    if tile.shape == (lentZ,lentY,lentX):    # well sized tile
        return tile.flatten()
    else:  
        new_tile = np.zeros( (lentZ,lentY,lentX),dtype="float32")
        new_tile[:tile.shape[0],:tile.shape[1],:tile.shape[2]] = tile
        return new_tile.flatten()


def tile_data3d(data,(lentZ,lentY,lentX)):
    """ 
    Tile sparky data into 1D numpy array

    Parameters:

    * data    Three-dimensional data array
    * lentZ   Z (w1) dimention tile size
    * lentY   Y (w2) dimention tile size
    * lentX   X (w3) dimention tile size

    Returns 1D numpy array of floats

    """
    

    # determind the number of tiles in data
    ttX = np.ceil(data.shape[2] / float(lentX))    # total tiles in X dim
    ttY = np.ceil(data.shape[1] / float(lentY))    # total tiles in Y dim
    ttZ = np.ceil(data.shape[0] / float(lentZ))    # total tiles in Z dim

    tt = ttX*ttY*ttZ    # total number of tiles
    
    # calc some basic parameter
    tsize = lentX*lentY*lentZ       # number of points in one tile
    t_tup = (lentZ,lentY,lentX)     # tile size tuple

    # create an empty array to store file data
    out = np.empty( (tt*tsize),dtype="float32")
    
    for i in xrange(int(tt)):
        out[i*tsize:(i+1)*tsize] = find_tilen_3d(data,i,t_tup)

    return out


def untile_data3D(data,(lentZ,lentY,lentX),(lenZ,lenY,lenX)):
    """ 
    Reorganize tiled sparky data into 3D data

    Parameters:

    * data    1D numpy array of tile data
    * lentZ   size of tile in Z (w1) dim
    * lentY   size of tile in Y (w2) dim
    * lentX   size of tile in X (w3) dim
    * lenZ    size of data in Z dim
    * lenY    size of data in Y dim
    * lenX    size of data in X dim

    Returns 3D np.array of floats with size (lenZ,lenY,lenX)

    """

    
    # determind the number of tiles in data
    ttX = np.ceil(lenX / float(lentX))    # total tiles in X dim
    ttY = np.ceil(lenY / float(lentY))    # total tiles in Y dim
    ttZ = np.ceil(lenZ / float(lentZ))    # total tiles in Z dim
    tt = ttX*ttY*ttZ

    # calc some basic parameter
    tsize = lentX*lentY*lentZ       # number of points in one tile
    t_tup = (lentZ,lentY,lentX)     # tile size tuple

    # create an empty array to store file data
    out = np.empty( (ttZ*lentZ,ttY*lentY,ttX*lentX),dtype="float32")

    for iZ in xrange(int(ttZ)):
        for iY in xrange(int(ttY)):
            for iX in xrange(int(ttX)):
            
                minX = iX*lentX
                maxX = (iX+1)*lentX

                minY = iY*lentY
                maxY = (iY+1)*lentY

                minZ = iZ*lentZ
                maxZ = (iZ+1)*lentZ

                ntile = iZ*ttX*ttY + iY*ttX + iX
                minT = ntile * tsize
                maxT = (ntile+1) * tsize

                out[minZ:maxZ,minY:maxY,minX:maxX] =  \
                data[minT:maxT].reshape(t_tup)

    return out[:lenZ,:lenY,:lenX]
 

# fileheader functions

def get_fileheader(f):
    """ 
    Get fileheader from file and return a list

    Reads the 180 byte file header of a sparky file

    """
    # file header as descriped in ucsffile.cc of sparky source
    # header is packed as follows: 
    # ident(10s),naxis(c),ncomponents(c),encoding(c),version(c)
    # owner(9s),date(26s),comment(80s),pad(3x),seek_pos(l),scratch(40s),
    # pad(4x)

    # note that between comment and seek_pos is a 3 byte pad
    # so that the long is @ a multiple of 4
    # also sparky always packs big-endian, hence > 

    return struct.unpack('>10s 4c 9s 26s 80s 3x l 40s 4x',f.read(180) )  


def put_fileheader(f,fl):
    """ 
    Write fileheader list to file (180-bytes)
    """

    f.write( struct.pack('>10s 4c 9s 26s 80s 3x l 40s 4x',*fl))
    return


def fileheader2dic(header):
    """ 
    Convert fileheader list into dictionary
    """

    dic = dict()

    dic["ident"]        = str(header[0]).strip('\x00')
    dic["naxis"]        = ord(header[1])
    dic["ncomponents"]  = ord(header[2])
    dic["encoding"]     = ord(header[3])
    dic["version"]      = ord(header[4])
    dic["owner"]        = str(header[5]).strip('\x00')
    dic["date"]         = str(header[6]).strip('\x00')
    dic["comment"]      = str(header[7]).strip('\x00')
    dic["seek_pos"]     = header[8]     # eof seek position
    dic["scratch"]      = str(header[9]).strip('\x00')
 
    return dic  


def dic2fileheader(dic):
    """ 
    Convert fileheader dictionary to list
    """

    fl = [0]*10
    fl[0]  = dic["ident"]  
    fl[1]  = chr(dic["naxis"])
    fl[2]  = chr(dic["ncomponents"])
    fl[3]  = chr(dic["encoding"])
    fl[4]  = chr(dic["version"])
    fl[5]  = dic["owner"]
    fl[6]  = dic["date"]
    fl[7]  = dic["comment"]
    fl[8]  = dic["seek_pos"]
    fl[9]  = dic["scratch"]

    return fl


# axisheader functions


def get_axisheader(f):
    """ 
    Get axisheader from file and return a list

    Only the first 44 bytes are examined, the NMR_PROCESSED and other header
    parameters are ignored since the current version of sparky does not use 
    them.

    """
    # axis header is described in ucsffile.cc
    # axis header is packed as follows
    # nucleus(6s),spectral_shift(h),npoints(I),size(I),bsize(I)
    # spectrometer_freq(f),spectral_width(f),xmtr_freq(f),zero_order(f),
    # first_order(f),first_pt_scale(f),ZEROS

    return struct.unpack('>6s h 3I 6f 84s',f.read(128) )

def put_axisheader(f,al):
    """ 
    Write axisheader list to file (128-bytes)
    """

    f.write( struct.pack('>6s h 3I 6f 84s',*al) )
    return


def axisheader2dic(header):
    """ 
    Convert axisheader list into dictionary
    """

    dic = dict()

    dic["nucleus"]        = str(header[0]).strip('\x00') 
    dic["spectral_shift"] = header[1]
    dic["npoints"] = header[2]
    dic["size"]    = header[3]
    dic["bsize"]   = header[4]
    dic["spectrometer_freq"] = header[5]
    dic["spectral_width"]    = header[6]
    dic["xmtr_freq"]         = header[7]
    dic["zero_order"]     = header[8]
    dic["first_order"]    = header[9]
    dic["first_pt_scale"] = header[10]
    dic["extended"] = header[11]

    return dic


def dic2axisheader(dic):
    """ 
    Convert axisheader dictionary to list
    """

    al = [0] * 12
    
    al[0]  = dic["nucleus"] 
    al[1]  = dic["spectral_shift"]
    al[2]  = dic["npoints"]
    al[3]  = dic["size"]
    al[4]  = dic["bsize"]
    al[5]  = dic["spectrometer_freq"]
    al[6]  = dic["spectral_width"]
    al[7]  = dic["xmtr_freq"]
    al[8]  = dic["zero_order"]
    al[9]  = dic["first_order"]
    al[10] = dic["first_pt_scale"]
    al[11] = dic["extended"]

    return al
