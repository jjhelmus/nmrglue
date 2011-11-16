"""
Functions for reading and writing NMRPipe binary files and table (.tab) files

NMRPipe file structure is described in the NMRPipe man pages and fdatap.h

""" 

# standard library modules
import struct 
import datetime
import os
from StringIO import StringIO

# external modules
import numpy as np

# nmrglue modules
import fileiobase

# these conversion functions are also accessable in the pipe module
from table import pipe2glue,glue2pipe,guess_pformat


#########################
# table reading/writing #
#########################

def read_table(filename):
    """
    Read a NMRPipe database table (.tab) file.

    Parameters:

    * filename  Name of file to read

    Returns: pcomments,pformat,rec

    * pcomments List of NMRPipe comment lines.
    * pformats  List of NMRPipe table column formats strings.
    * rec       Records array with named fields.

    """

    # divide up into comment lines and data lines
    specials = ["VARS","FORMAT","NULLSTRING","NULLVALUE","REMARK","DATA"]
    f = open(filename,'r')
    cl = []
    dl = []
    for line in f:
        for k in specials:
            if line[:len(k)]==k:
                cl.append(line)
                break
        else:
            dl.append(line)
    f.close()
    
    # pull out and parse the VARS line
    vl = [i for i,l in enumerate(cl) if l[:4]=="VARS"]
    if len(vl)!=1:
        raise IOError("%s has no/more than one VARS line"%(filename))
    dtd = {'names':cl.pop(vl[0]).split()[1:]}

    # pull out and parse the FORMAT line
    fl = [i for i,l in enumerate(cl) if l[:6]=="FORMAT"]
    if len(fl)!=1:
        raise IOError("%s has no/more than one FORMAT line"%(filename))
    pformat = cl.pop(fl[0]).split()[1:]
    p2f = {'d':'i4','f':'f8','e':'f8','s':'S256'}   # pipe -> format dict.
    dtd['formats'] = [p2f[i[-1]] for i in pformat]

    # DEBUG
    #print  dtd['names'],dtd['formats']

    s = StringIO("".join(dl))

    rec = np.recfromtxt(s,dtype=dtd,comments='XXXXXXXXXXX')
    return cl,pformat,np.atleast_1d(rec)


def write_table(filename,pcomments,pformats,rec,overwrite=False):
    """
    Write a NMRPipe database table (.tab) file.

    Parameters:

    * filename  Name of file to write to.
    * pcomments List of NMRPipe comment lines.
    * pformats  List of NMRPipe table column formats strings.
    * rec       Records array of table.
    * overwrite Set True to overwrite file if it exists.

    """
    if len(rec[0])!=len(pformats):
        s = "number of rec columns %i and pformat elements %i do not match"
        raise ValueError( s%(len(rec[0]),len(pformats) ) )
    
    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite)

    # write out the VARS line
    names = rec.dtype.names
    s = "VARS   "+" ".join(names)+"\n"
    f.write(s)

    # write out the FORMAT line
    s = "FORMAT "+" ".join(pformats)+"\n"
    f.write(s)

    # write out any comment lines
    for c in pcomments:
        f.write(c)

    # write out each line of the records array
    s = " ".join(pformats)+"\n"
    for row in rec:
        f.write(s%tuple(row))

    f.close()
    return

###################
# unit conversion #
###################

def make_uc(dic,data,dim=-1):
    """ 
    Make a unit conversion object
    """

    if dim == -1:
        dim = data.ndim - 1 # last dimention 

    fn = "FDF" + str(int(dic["FDDIMORDER"][data.ndim-1-dim]))
    size = float(data.shape[dim])
    
    # check for quadrature in indirect dimentions
    if (dic[fn+"QUADFLAG"] != 1) and (dim !=data.ndim-1):
        size = size/2.
        cplx = True
    else:
        cplx = False

    sw = dic[fn+"SW"]
    if sw == 0.0:   
        sw = 1.0
    obs = dic[fn+"OBS"]
    if obs == 0.0:
        obs = 1.0

    car = dic[fn+"CAR"]*obs

    # NMRPipe keeps the carrier constant during extractions storing the 
    # location of this point as CENTER.  This may not be the actual "center" of
    # the spectrum and may not even be a valid index in that dimension. We need
    # to re-center the carrier value so that actually represents the 
    # frequency of the central point in the dimension.
    car = car + sw/size * ( dic[fn+"CENTER"]-1.-size/2. )
    

    return fileiobase.unit_conversion(size,cplx,sw,obs,car)


############################
# dictionary/data creation #
############################

fd2dphase_dic = {"magnitude":0,"tppi":1,"states":2,"image":3}

def create_data(data):
    """ 
    Create a NMRPipe data array (recast into float32 or complex64)
    """

    if np.iscomplexobj(data):   # check quadrature
        return np.array(data,dtype="complex64")
    else:
        return np.array(data,dtype="float32")


########################
# universal dictionary #
########################

def guess_udic(dic,data):
    """ 
    Guess parameter of universal dictionary from dic,data pair
    """

    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in xrange(data.ndim):
         
        udic[i]["size"] = data.shape[i] # size from data shape
        
        # determind NMRPipe axis name
        fn = ["FDF2","FDF1","FDF3","FDF4"][(data.ndim-1)-i]
        
        # directly corresponding
        udic[i]["sw"] = dic[fn+"SW"]
        udic[i]["obs"] = dic[fn+"OBS"] 
        udic[i]["car"] = dic[fn+"CAR"]*dic[fn+"OBS"] # ppm->hz 
        udic[i]["label"] = dic[fn+"LABEL"]

        if dic[fn+"QUADFLAG"] == 1: # real data
            udic[i]["complex"] = False
        else:
            udic[i]["complex"] = True

        if dic[fn+"FTFLAG"] == 0: # time domain
            udic[i]["time"] = True
            udic[i]["freq"] = False
        else:
            udic[i]["time"] = False
            udic[i]["freq"] = True

        if i != 0:
            if  dic["FD2DPHASE"] == 0:
                udic[i]["encoding"] = "unknown" # XXX magnitude
            elif dic["FD2DPHASE"] == 1:
                udic[i]["encoding"] = "tppi"
            elif dic["FD2DPHASE"] == 2:
                udic[i]["encoding"] = "states"
            elif dic["FD2DPHASE"] == 2:
                udic[i]["encoding"] = "unknown" # XXX image
            else:
                udic[i]["encoding"] = "unknown"

    return udic


def create_dic(udic,datetimeobj=datetime.datetime.now()):
    """ 
    Crate a NMRPipe dictiony from universal dictionary
    
    Parameters:

    * udic        Universal dictionary
    * datetimeobj datetime object
    * user        Name of user
   
    Does not update dictionary keys that are unknown such as MIN/MAX,
    apodization and processing parameters, sizes in none-current domain. 
    Also rounding of parameter is different than NMRPipe.

    """

    # create the blank dictionary
    dic = create_empty_dic()    # create the empty dictionary
    dic = datetime2dic(datetimeobj,dic) # add the datetime to the dictionary

    # fill global dictionary parameters
    dic["FDDIMCOUNT"] = float(udic["ndim"])

    # FD2DPHASE
    if udic[0]["encoding"] == "tppi":
        dic["FD2DPHASE"] = 1.0
    elif udic[0]["encoding"] == "states":
        dic["FD2DPHASE"] = 2.0
    else:
        dic["FD2DPHASE"] = 0.0    
 
    # fill in parameters for each dimension
    for i,adic in enumerate([udic[k] for k in xrange(udic["ndim"])]):
        n = int((dic["FDDIMCOUNT"]-1)-i)
        dic = add_axis_to_dic(dic,adic,n)

   
    if dic["FDDIMCOUNT"] >= 3: # at least 3D    
        dic["FDFILECOUNT"] = dic["FDF3SIZE"] * dic["FDF4SIZE"]

    if (dic["FDF1QUADFLAG"]==dic["FDF2QUADFLAG"]==dic["FDF3QUADFLAG"]) and (
       dic["FDF1QUADFLAG"]==dic["FDF4QUADFLAG"]==1):
        dic["FDQUADFLAG"] = 1.0


    return dic


def add_axis_to_dic(dic,adic,n):
    """ 
    Add an axis to NMRPipe dictionary

    n is 0,1,2,... (0 is direct dim, 1 first indirect...)

    """

    # determind F1,F2,F3,...
    fn = ["FDF2","FDF1","FDF3","FDF4"][n]
    
    # parameter directly in dictionary
    dic[fn+"SW"]   = float(adic["sw"])
    dic[fn+"OBS"]  = float(adic["obs"])
    dic[fn+"CAR"]  = float(adic["car"]/adic["obs"])
    dic[fn+"LABEL"] = adic["label"]

    if adic["complex"]:
        dic[fn+"QUADFLAG"]  = 0.0
    else:
        dic[fn+"QUADFLAG"]  = 1.0

    
    # determine R|I size
    if adic["complex"] and n!=0:
        psize = adic["size"]/2.
    else:
        psize = adic["size"]/1.

    # origin calculation size
    osize = psize

    # set FT/TD SIZE and FTFLAG depending on domain
    if adic["time"]:    
        dic[fn+"TDSIZE"] = psize
        dic[fn+"FTFLAG"] = 0.0
    else:
        dic[fn+"FTSIZE"] = psize
        dic[fn+"FTFLAG"] = 1.0
    
    # apodization and center
    dic[fn+"APOD"]   = dic[fn+"TDSIZE"]
   
    if n==0 or dic["FD2DPHASE"]!=1:
        dic[fn+"CENTER"] = int(psize / 2.)+1.
    else:   # TPPI requires division by 4
        dic[fn+"CENTER"] = int(psize/4.)+1
        osize = psize/2.

    # origin (last point) is CAR*OBS-SW*(N/2-1)/N
    # see Fig 3.1 on p.36 of Hoch and Stern
    #print "fn:",n
    #print  "CAR:",dic[fn+"CAR"]
    #print  "OBS:",dic[fn+"OBS"]
    #print  "SW:",dic[fn+"SW"]
    #print  "osize:",osize
    #print  "CENTER:",dic[fn+"CENTER"]
    dic[fn+"ORIG"] = dic[fn+"CAR"]*dic[fn+"OBS"] - dic[fn+"SW"] * \
        (osize-dic[fn+"CENTER"])/osize
 
    if n==0:  # direct dim
        dic["FDSIZE"]     = psize
        dic["FDREALSIZE"] = psize
    
    if n==1:  # first indirect
        dic["FDSPECNUM"] = float(adic["size"]) # R+I
    
    if n==2:  # second indirect
        if adic["complex"]:
            dic["FDF3SIZE"] = psize*2
        else:
            dic["FDF3SIZE"] = psize
        
    if n==3:  # third indirect
        if adic["complex"]:
            dic["FDF4SIZE"] = psize*2
        else:
            dic["FDF3SIZE"] = psize
    return dic


def create_empty_dic():
    """ 
    Creates a nmrpipe dictionary with default values
    """

    dic =  fdata2dic(np.zeros((512),dtype="float32"))

    # parameters which are 1
    dic["FDF1CENTER"] = 1.
    dic["FDF2CENTER"] = 1.
    dic["FDF3CENTER"] = 1.
    dic["FDF4CENTER"] = 1.

    dic["FDF3SIZE"] = 1.
    dic["FDF4SIZE"] = 1.

    dic["FDF1QUADFLAG"] = 1.
    dic["FDF2QUADFLAG"] = 1.
    dic["FDF3QUADFLAG"] = 1.
    dic["FDF4QUADFLAG"] = 1.

    dic["FDSPECNUM"] = 1.
    dic["FDFILECOUNT"] = 1.
    dic["FD2DVIRGIN"] = 1.
    # dimention ordering

    dic["FDDIMORDER1"] = 2.0
    dic["FDDIMORDER2"] = 1.0
    dic["FDDIMORDER3"] = 3.0
    dic["FDDIMORDER4"] = 4.0
    dic["FDDIMORDER"] = [2.0,1.0,3.0,4.0]

    # string and such
    dic["FDF1LABEL"] = "Y"
    dic["FDF2LABEL"] = "X"
    dic["FDF3LABEL"] = "Z"
    dic["FDF4LABEL"] = "A"

    # misc values
    dic["FDFLTFORMAT"] = struct.unpack('f','\xef\xeenO')[0]
    dic["FDFLTORDER"] = float(2.3450000286102295)

    return dic


def datetime2dic(dt,dic):
    """ 
    Add datatime object to dictionary
    """

    dic["FDYEAR"]  = float(dt.year)
    dic["FDMONTH"] = float(dt.month)
    dic["FDDAY"]   = float(dt.day)

    dic["FDHOURS"] = float(dt.hour)
    dic["FDMINS"]  = float(dt.minute)
    dic["FDSECS"]  = float(dt.second)

    return dic


def dic2datetime(dic):
    """ 
    Create a datetime object from dictionary
    """
    year   = int(dic["FDYEAR"])
    month  = int(dic["FDMONTH"])
    day    = int(dic["FDDAY"])
    hour   = int(dic["FDHOURS"])
    minute = int(dic["FDMINS"])
    second = int(dic["FDSECS"])

    return datetime.datetime(year,month,day,hour,minute,second)

################
# file reading #
################

def read(filename):
    """
    Read a NMRPipe binary file returning a dic,data pair.

    For standard multi-file 3D/4D NMRPipe data sets, filename should be a 
    filemask (for example "/ft/test%03d.ft3") with a "%" formatter.  If only
    one file of a 3D/4D data set is provided only that 2D slice of the data is
    read (for example "/ft/test001.ft3" results in a 2D data set being read).

    NMRPipe data streams stored as files (one file 3D/4D data sets made using
    xyz2pipe) can be read by providing the file name of the stream.  The entire
    data set is read into memory.

    """
    if filename.count("%")==1:
        filemask = filename
        filename = filename%1
    elif filename.count("%")==2:
        filemask = filename
        filename = filename%(1,1)
    else:
        filemask = None

    fdata = get_fdata(filename)
    dic = fdata2dic(fdata)  
    order = dic["FDDIMCOUNT"]

    if order == 1:
        return read_1D(filename)
    if order == 2:
        return read_2D(filename)
    if dic["FDPIPEFLAG"] != 0:  # open streams
        return read_stream(filename)
    if filemask == None:     # if no filemask open as 2D
        return read_2D(filename)
    if order == 3:
        return read_3D(filemask)
    if order == 4:
        return read_4D(filemask)
    
    raise ValueError,'unknown dimensionality: %s'%order


def read_lowmem(filename):
    """
    Read a NMRPipe binary file with minimal memory usage.

    For standard multi-file 3D/4D NMRPipe data sets, filename should be a
    filemask (for example "/ft/test%03d.ft3") with a "%" formatter.  If only
    one file of a 3D/4D data set is provided only that 2D slice of the data is
    read (for example "/ft/test001.ft3" results in a 2D data set being read).

    NMRPipe data streams stored as files (one file 3D/4D data sets made using
    xyz2pipe) can be read by providing the file name of the stream.  

    """
    if filename.count("%")==1:
        filemask = filename
        filename = filename%1
    elif filename.count("%")==2:
        filemask = filename
        filename = filename%(1,1)
    else:
        filemask = None

    fdata = get_fdata(filename)
    dic = fdata2dic(fdata)  
    order = dic["FDDIMCOUNT"]
    
    if order == 1:
        return read_1D(filename)    # there is no 1D low memory option
    if order == 2:
        return read_lowmem_2D(filename)
    if dic["FDPIPEFLAG"] != 0:  # open streams
        return read_lowmem_stream(filename)
    if filemask == None:    # if no filemask open as 2D
        return read_lowmem_2D(filename)
    if order == 3:
        return read_lowmem_3D(filemask)
    if order == 4:
        return read_lowmem_4D(filemask)
    
    raise ValueError,'unknown dimentionality: %s'%order


# dimension specific reading

def read_1D(filename):
    """ 
    Read a 1D NMRPipe binary file returning a dic,data pair
    """
    fdata,data = get_fdata_data(filename)   # get the fdata and data arrays
    dic = fdata2dic(fdata)  # convert the fdata block to a python dictionary    
    data = reshape_data(data,find_shape(dic))    # reshape data
    
    # unappend imaginary data if needed
    if dic["FDF2QUADFLAG"] != 1:
        data = unappend_data(data)

    return (dic,data)

def read_2D(filename):
    """
    Read a 2D NMRPipe file or NMRPipe data stream returning a dic,data pair
    """
    fdata,data = get_fdata_data(filename)   # get the fdata and data arrays
    dic = fdata2dic(fdata)  # convert the fdata block to a python dictionary
    data = reshape_data(data,find_shape(dic))    # reshape data

    # unappend imaginary data if needed
    if dic["FDTRANSPOSED"] == 1 and dic["FDF1QUADFLAG"] != 1:
        data = unappend_data(data)
    elif dic["FDTRANSPOSED"] == 0 and dic["FDF2QUADFLAG"] != 1:
        data = unappend_data(data)

    return (dic,data)
        
def read_lowmem_2D(filename):
    """
    Read a 2D NMRPipe file or NMRPipe data stream with minimal memory usage
    """
    dic = fdata2dic(get_fdata(filename))
    order = dic["FDDIMCOUNT"]
    if order == 2:
        data = pipe_2d(filename)
    if order == 3:
        data = pipestream_3d(filename)
    if order == 4:
        data = pipestream_4d(filename)
    return dic,data

def read_stream(filename):
    """
    Read a NMRPipe data stream (one file 3D or 4D files)
    """
    return read_2D(filename)

def read_lowmem_stream(filename):
    """
    Read a NMRPipe data stream with minimal memory usage
    """
    return read_lowmem_2D(filename)

def read_3D(filemask):
    """
    Read a 3D NMRPipe binary file returning a dic,data pair.
    """
    dic,data = read_lowmem_3D(filemask)
    data = data[:,:,:]  # read all the data
    return dic,data

def read_lowmem_3D(filemask):
    """ 
    Read a 3D NMRPipe binary file with minimal memory usage.
    """
    if '%' not in filemask: # data streams should be read with read_stream
        return read_lowmem_stream(filemask)
    data = pipe_3d(filemask)    # create a new pipe_3d object
    dic = fdata2dic(get_fdata(filemask%(1)))
    return (dic,data)

def read_4D(filemask):
    """
    Read a 3D NMRPipe binary file returning a dic,data pair.

    This function should not be used to read NMRPipe data streams stored in a
    single file (one file 3D/4D data sets made using xyz2pipe), read_2D
    should be used.
    """
    dic,data = read_lowmem_4D(filemask)
    data = data[:,:,:,:]  # read all the data
    return dic,data

def read_lowmem_4D(filemask):
    """ 
    Read a NMRPipe binary file with minimal memory usage.

    This function should not be used to read NMRPipe data streams stored in a
    single file (one file 3D/4D data sets made using xyz2pipe), read_lowmem_2D
    should be used.
    """
    if '%' not in filemask: # data streams should be read with read_stream
        return read_lowmem_stream(filemask)

    data = pipe_4d(filemask)    # create a new pipe_3d object
    if data.singleindex:
        dic = fdata2dic(get_fdata(filemask%(1)))
    else:
        dic = fdata2dic(get_fdata(filemask%(1,1)))
    return (dic,data)

#####################
# writing functions #
#####################

def write(filename,dic,data,overwrite=False):
    """
    Write a NMRPipe file 

    For 3D data if filename has no '%' formatter then the data is written as a 
    3D NMRPipe data stream.  When the '%' formatter is provided the data is 
    written out as a standard NMRPipe 3D multi-file 3D.

    For 4D data, filename can have one, two or no '%' formatters resulting in
    a single index file (test%03d.ft), two index file(test%02d%03d.ft), or 
    one file data stream (test.ft4).

    dic["FDPIPEFLAG"] is not changed or checked when writing, please check
    that this value is 0.0 for standard non-data stream files, and 1.0 for data
    stream files or an file may be written with an incorrect header.

    Set overwrite to True to overwrite files that exist.

    """
    # load all data if the data is not a numpy ndarray
    if type(data) != np.ndarray:
        data = data[:]

    if filename.count("%")==0:
        return write_single(filename,dic,data,overwrite)
    elif data.ndim==3:
        return write_3D(filename,dic,data,overwrite)
    elif data.ndim==4:
        return write_4D(filename,dic,data,overwrite)
    
    raise ValueError,'unknown filename/dimension'

def write_single(filename,dic,data,overwrite=False):
    """
    Write data to a single NMRPipe file from memory.  This write 1D and 2D 
    files completely as well as NMRPipe data streams. 2D planes of 3D and
    4D files should be written with this function.

    """
    # append imaginary and flatten
    if data.dtype=="complex64":
        data = append_data(data)
    data = unshape_data(data)

    # create the fdata array
    fdata = dic2fdata(dic)

    # write the file
    put_data(filename,fdata,data,overwrite)
    return

def write_3D(filemask,dic,data,overwrite=False):
    """
    Write a standard multi-file 3D NMRPipe file
    """
    lenZ,lenY,lenX = data.shape
    for zi in range(lenZ):
        fn = filemask%(zi+1)
        plane  = data[zi]
        write_single(fn,dic,plane,overwrite)
    return

def write_4D(filemask,dic,data,overwrite=False):
    """
    Write a one or two index 4D NMRPipe file
    """
    lenA,lenZ,lenY,lenX = data.shape
    for ai in range(lenA):
        for zi in range(lenZ):

            if filemask.count("%")==2:
                fn = filemask%(ai+1,zi+1)
            else:
                fn = filemask%(ai*lenZ+zi+1)
            
            plane = data[ai,zi]
        
            # update dictionary if needed
            if dic["FDSCALEFLAG"] == 1:
                dic["FDMAX"]     = plane.max()
                dic["FDDISPMAX"] = dic["FDMAX"]
                dic["FDMIN"]     = plane.min()
                dic["FDDISPMIN"] = dic["FDMIN"]
            write_single(fn,dic,plane,overwrite)
    return

def write_lowmem(filename,dic,data,overwrite=False):
    """ 
    Write a NMRPipe file using minimal memory
    """

    if data.ndim==1:
        return write_single(filename,dic,data,overwrite)
    if data.ndim==2:
        return write_lowmem_2D(filename,dic,data,overwrite)
    if data.ndim==3:
        if "%" in filename:
            return write_lowmem_3D(filename,dic,data,overwrite)
        else:
            return write_lowmem_3Ds(filename,dic,data,overwrite)
    if data.ndim==4:
        if "%" in filename:
            return write_lowmem_4D(filename,dic,data,overwrite)
        else:
            return write_lowmem_4Ds(filename,dic,data,overwrite)
    
    raise ValueError,'unknown dimensionality: %s'%data.ndim
    

def write_lowmem_2D(filename,dic,data,overwrite=False):
    """
    Write a 2D NMRPipe file using minimal memory (trace by trace)
    """
    fh = fileiobase.open_towrite(filename,overwrite=overwrite)
    
    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh,fdata)

    # put data trace by trace
    lenY,lenX = data.shape
    for y in xrange(lenY):
        put_trace(fh,data[y])
    fh.close()
    return 

def write_lowmem_3D(filename,dic,data,overwrite=False):
    """
    Write a standard multi-file 3D NMRPipe file using minimal memory 
    (trace by trace)

    MIN/MAX parameters are not updated in the NMRPipe headers.
    """
    # create the fdata array
    fdata = dic2fdata(dic)

    # put data trace by trace
    lenZ,lenY,lenX = data.shape
    for z in xrange(lenZ):
        # open the file to store the 2D plane
        fh = fileiobase.open_towrite(filename%(z+1),overwrite=overwrite)
        put_fdata(fh,fdata)
        for y in xrange(lenY):
            put_trace(fh,data[z,y])
        fh.close()
    return 


def write_lowmem_3Ds(filename,dic,data,overwrite=False):
    """
    Write 3D NMRPipe data stream file using minimal memory (trace by trace)
    """
    fh = fileiobase.open_towrite(filename,overwrite=overwrite)

    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh,fdata)

    # put data trace by trace
    lenZ,lenY,lenX = data.shape
    for z in xrange(lenZ):
        for y in xrange(lenY):
            put_trace(fh,data[z,y])
    fh.close()
    return 


def write_lowmem_4D(filename,dic,data,overwrite=False):
    """
    Write a standard multi-file (single or double index 4D NMRPipe file using 
    minimal memory (trace by trace).  

    MIN/MAX parameters are not updated in the NMRPipe headers.
    """
    # create the fdata array
    fdata = dic2fdata(dic)

    # put data trace by trace
    lenA,lenZ,lenY,lenX = data.shape
    for a in xrange(lenA):
        for z in xrange(lenZ):
            # open the file to store the 2D plane
            if filename.count("%")==1:
                fname = filename%(a*lenZ+z+1)
            else:
                fname = filename%(a+1,z+1)
            fh = fileiobase.open_towrite(fname,overwrite=overwrite)
            put_fdata(fh,fdata)
            for y in xrange(lenY):
                put_trace(fh,data[a,z,y])
            fh.close()
    return 


def write_lowmem_4Ds(filename,dic,data,overwrite=False):
    """
    Write 4D NMRPipe data stream file using minimal memory (trace by trace)
    """
    fh = fileiobase.open_towrite(filename,overwrite=overwrite)

    # create the fdata array and put to disk
    fdata = dic2fdata(dic)
    put_fdata(fh,fdata)

    # put data trace by trace
    lenA,lenZ,lenY,lenX = data.shape
    for a in xrange(lenA):
        for z in xrange(lenZ):
            for y in xrange(lenY):
                put_trace(fh,data[a,z,y])
    fh.close()
    return 

###############
# put to disk #
###############


def put_fdata(fh,fdata):
    """
    Put fdata to NMRPipe file described by file object fh
    """
    if fdata.dtype != 'float32':
        raise TypeError,'fdata.dtype is not float32'
    fh.write(fdata.tostring())
    return


def put_trace(fh,trace):
    """
    Put 1D trace (real or complex) to NMRPipe file described by file object fh
    """
    if trace.dtype == 'complex64':
        trace = append_data(trace)
    if trace.dtype != 'float32':
        raise TypeError,'trace.dtype is not float32'

    fh.write(trace.tostring())
    return


def put_data(filename,fdata,data,overwrite=False):
    """
    Put fdata and data to 2D NMRPipe.
    """

    if data.dtype != 'float32':
        print data.dtype
        raise TypeError,'data.dtype is not float32'
    if fdata.dtype != 'float32':
        raise TypeError,'fdata.dtype is not float32'

    # write the file
    f = fileiobase.open_towrite(filename,overwrite=overwrite)
    f.write(fdata.tostring())
    f.write(data.tostring())
    f.close()
    return


def write_slice_3D(filemask,dic,data,shape,(sz,sy,sx) ):
    """ 
    Write a slice of a 3D data array to file

    Opens (or if necessary creates) 2D NMRPipe file(s) to write 
    data, where total 3D file size is given by shape.

    Parameters:
    * filemask      String with single formatting operator (%)
    * data          3D array of data
    * dic           Dictionary to write when/if files are created
    * shape         3-tuple of integers indicating the overall matrix shape
    * (sz,sy,sx)    3-tuple of slice object which specify location of data

    This function memmaps 2D NMRPipe files for speed. It only writes 
    dictionaries to file when created, leaving them unmodified if the file
    exists.  
    
    Only error checking is that data is 3D. 

    Users are not expected to use this function, rather use the iter3D object

    """
    
    if data.ndim != 3:
        raise ValueError,"passed array must be 3D"

    # unpack the shape
    dz,dy,dx = shape
    
    # create list of file names
    fnames = [filemask % i for i in range(1,dz+1)]

    # loop over the requested z-slice
    for i,f in enumerate(fnames[sz]):
    
        #print "i:",i,"f:",f

        if os.path.isfile(f) == False:
            # file doesn't exist, create a empty one        
            ndata = np.zeros( (dy,dx),dtype=data.dtype)
            write_2D(f,dic,data,False)
            del(ndata)
        
        # mmap the [new] file
        mdata = np.memmap(f,dtype='float32',offset=512*4,mode='r+')
        # reshape 
        mdata = mdata.reshape((dy,dx))
            
        # unpack into rdata,[idata] depending on quadrature
        if data.dtype == 'complex64':
            h = mdata.shape[-1]/2.0
            rdata = mdata[...,:h]
            idata = mdata[...,h:]
        else:
            rdata = mdata

        # write the data out, flush and close
        rdata[sy,sx] = data.real[i]
        rdata.flush()
        if data.dtype == 'complex64':
            idata[sy,sx] = data.imag[i]
            idata.flush()
            del(idata)

        # clean up
        del(rdata)
        del(mdata)


# iter3D tools (xyz2pipe and pipe2xyz replacements)


#Notes for iter3D implementation
#
#'x'/'y' in_lead
#==============
#Reading
#-------
#- When passed x must transposed 1,2 if dic["FDTRANSPOSED"] == 1
# (might need to call pipe_proc.tp)
#- if 'y' passed then cann pipe_proc.tp unless dic["FDTRANSPOED"]
#- save 'good' dictionary and return each loop
#
#Looping
#-------
#- will loop until data.shape[0] reached
#- returns dic, XY or YX plane
#
#Writing
#-------
#- if 'y' out then need final pipe_proc.tp of data, if 'x' do nothing
#- reshape data to 1,plane.shape[0],plane.shape[1]
#- size becomes data.shape[0],plane.shape[0],plane.shape[1]
#- sz = slice(i,i+1,1) sy=sx=slice(None)
#
#
#'z' in_lead
#===========
#Reading
#-------
#- Untranspose if dic["TRANSPOSED"] == 1 (call pipe_proc.tp)
#- transpose (1,2,0)
#- ORDER 1,2,3 = 3,1,2 and array
#- update "FDSLICECOUNT" and "FDSIZE" taking into accound complex packing
#- also update "FDSPECNUM"
#- call write_slice3D
#- store shape as self.max_iter
#
#Looping
#-------
#- grab the slice and pack_complex if needed
#- returns dic,ZX-plane
#
#Writing
#-------
#- if out_lead = 'x' needs final pipe_proc.tp of data, if 'z' do nothing
#- reshape data to 1,plane.shape[0],plane.shape[1]
#- transposed data to 2,0,1 (or combine with above step
#- update "FDSIZE" and "FDSPECNUM"
#- remove min/max
#- update FDDIMORDER and ORDER1,2,3
#- size plane.shape[0],self.max_iter,plane.shape[2]
#- sz = slice(None)=sx
#- sy = slice(i,i+1,1)


def pack_complex(data):
    """
    Pack inteleaved real,imag array into complex array
    """
    return np.array(data[...,::2]+data[...,1::2]*1.j,dtype="complex64")


def transpose_3D(dic,data,(a1,a2,a3)=(2,1,0) ):
    """ 
    Transpose pipe_3d object and dictionary
    """

    rdic = dict(dic)    # create a copy of the dictionary
    
    # transpose the data
    data = data.transpose( (a1,a2,a3) )

    # transpose the dictionary
    s3 = "FDDIMORDER"+str(int(3-a1))    # 3rd axis is 0th axis in data_nd
    s2 = "FDDIMORDER"+str(int(3-a2))    # 2nd axis is 1st axis in data_nd
    s1 = "FDDIMORDER"+str(int(3-a3))    # 1st axis is 3nd axis in data_nd

    rdic["FDDIMORDER1"] = dic[s1]
    rdic["FDDIMORDER2"] = dic[s2]
    rdic["FDDIMORDER3"] = dic[s3]

    rdic['FDDIMORDER'] = [ rdic["FDDIMORDER1"], rdic["FDDIMORDER2"],
                           rdic["FDDIMORDER3"], rdic["FDDIMORDER4"] ]

    # set the shape dictionary parameters
    fn = "FDF"+str(int(rdic["FDDIMORDER1"]))
    if rdic[fn+"QUADFLAG"] != 1.0:   # last axis is complex
        rdic["FDSIZE"] = data.shape[2]/2.
    else:   # last axis is singular
        rdic["FDSIZE"] = data.shape[2]

    rdic["FDSLICECOUNT"] = data.shape[1]
    rdic["FDSPECNUM"] = rdic["FDSLICECOUNT"]

    return rdic,data


class iter3D(object):
    """ 
    Object which allows for graceful iteration over 3D NMRPipe files
    
    iter3D.iter() returns a (dic,plane) tuple which can be written using
    the x.writeplane function.
    
    When processing 3D files with iter3D object(s) the following dictionary 
    parameters may not have the same values as NMRPipe processing scripts 
    return:
    
    FDSLICECOUNT and
    
    FDMAX,FDDISMAX,FDMIN,FDDISPMIN when FDSCALEFLAG == 0
    
    Example::

        #3D data processing
        xiter = iter3D("data/test%03d.fid","x","x")
        for dic,YXplane in xiter:
            # process X and Y axis
            xiter.write("ft/test%03d.ft2",YXplane,dic)
        ziter = iter3D("ft/test%03d.ft2","z","z")
        for dic,XZplane in ziter:
            # process Z axis
            ziter.write("ft/test%03d.ft3",XZplane,dic)
    
    """

    def __init__(self,filemask,in_lead="x",out_lead="DEFAULT"):
        """
        Create a iter3D object

        Parameters:
        * filemask  string with single formatter (%) of NMRPipe files to read
        * in_lead   Axis name ('x','y','z') of last (1st) axis in outputed 2D 
        * out_lead  Axis name ('x','y','z') of axis to be written typically
                    this is the same as in_lead

        =======     ===============
        In-lead     Iterated Planes
        =======     ===============
        "x"         ('y','x')
        "y"         ('x','y')
        "z"         ('x','z') 
        =======     ===============

        """

        # check for invalid in_lead, out_lead
        if in_lead not in ["x","y","z"]:
            raise ValueError,"in_lead must be 'x','y' or 'z'"

        if out_lead not in ["x","y","z","DEFAULT"]:
            raise ValueError,"out_lead must be 'x','y','z' or 'DEFAULT'"

        if out_lead == "DEFAULT":
            out_lead = in_lead

        if in_lead in ["x","y"] and out_lead not in ["x","y"]:
            raise ValueError,"Invalid in_lead, out_lead pair"

        if in_lead == "z" and out_lead not in ["x","z"]:
            raise ValueError,"Invalid in_lead, out_lead pair"

        self.in_lead  = in_lead
        self.out_lead = out_lead

        self.dic,self.pipe_3d = read_3D(filemask)
  
        # uptranspose data if needed
        if self.dic["FDTRANSPOSED"] == 1.0:
            # need to switch X and Y (0,2,1)
            self.dic,self.pipe_3d = transpose_3D(self.dic,self.pipe_3d,(0,2,1))

        # self.pipe_3d and self.dic are now REALLY ZYX order
       
        # now prep pipe_3d for slicing and be make
        # idic the iterator dictionary

        self.i = -1  # counter

        if self.in_lead == "x":
            # leave as is Z(YX)
            self.needs_pack_complex = False
            self.idic = dict(self.dic)
            self.i_max = int(self.pipe_3d.shape[0])

        elif self.in_lead == "y":
            # transpose to Z(XY)
            self.idic,self.pipe_3d = transpose_3D(self.dic,self.pipe_3d,(0,2,1))
            self.needs_pack_complex = False
            self.i_max = int(self.pipe_3d.shape[0])

        elif self.in_lead == "z":
            # transpose to Y(XZ)
            self.idic,self.pipe_3d = transpose_3D(self.dic,self.pipe_3d,(1,2,0))
            fn = "FDF"+str(int(self.idic["FDDIMORDER1"]))
            if self.idic[fn+"QUADFLAG"] != 1.0:   # z axis is complex
                self.needs_pack_complex = True
            else:
                self.needs_pack_complex = False
            self.i_max = int(self.pipe_3d.shape[0])
        else:
            raise Error,"You should NEVER get here"


    def __iter__(self):
        """ 
        x.__iter__() <==> iter(x)
        """
        return self

    def next(self):
        """ 
        Return the next dic,plane or raise StopIteration
        """
        self.i = self.i + 1
        if self.i >= self.i_max:
            raise StopIteration
        else:
            plane = self.pipe_3d[self.i]
            if self.needs_pack_complex:
                plane = pack_complex(plane)
            return (dict(self.idic),plane)

    def reinitialize(self):
        """ 
        Restart iterator at first dic,plane
        """
        self.i = -1


    def write(self,filemask,plane,dic):
        """ 
        Write out current plane
        """

        # make the plane a 3D array
        plane = plane.reshape(1,plane.shape[0],plane.shape[1])

        if self.in_lead != self.out_lead:
            # transpose the last two axes
            dic,plane = transpose_3D(dic,plane,(0,2,1))
            

        if self.in_lead == "x" or self.in_lead=="y":
            shape = ( self.i_max,plane.shape[1],plane.shape[2] )
            sz = slice(self.i,self.i+1,1)
            sx = slice(None)
            sy = slice(None)

        elif self.in_lead == "z":
            # reorder from YXZ -> ZYX
            dic,plane = transpose_3D(dic,plane,(2,0,1))
            
            # turn scale flag off
            dic["FDSCALEFLAG"] = 0.0
            # the Y size is incorrect
            dic["FDSPECNUM"] = self.i_max

            # update the file count XXX these should be done bettwe
            dic["FDFILECOUNT"] = plane.shape[0]
            dic["FDF3SIZE"] = plane.shape[0]

            shape = ( plane.shape[0],self.i_max,plane.shape[2] )
            sx = slice(None)
            sy = slice(self.i,self.i+1,1)
            sz = slice(None)

        else:
            raise Error,"You should NEVER get here"

        #print "Writing out slice :",self.i
        #print "shape:",shape
        #print "plane.shape",plane.shape
        #print "sx,sy,sz",sx,sy,sz
        #print dic["FDFILECOUNT"]
        write_slice_3D(filemask,dic,plane,shape,(sz,sy,sx) )

#####################
# Shaping functions #
#####################

def find_shape(dic):
    """ 
    Find the shape (tuple) of data in a NMRPipe file from dictionary
    
    1-tuple is returned for 1D data, 2-tuple for 2D and non-stream 3D/4D data,
    3-tuple or 4-tuple for stream 3D/4D data.

    The last dimension of the tuple is length of the data in the file, the
    actual length of the data matrix may be half of this if the data is
    complex.
    
    """

    if dic["FDDIMCOUNT"] == 1: # 1D Data

        if dic["FDF2QUADFLAG"] == 1:
            multi = 1.0
        else:
            multi = 2.0

        dim1 = int(dic["FDSIZE"]*multi)
        return (dim1)

    else: # 2D+ Data
        
        if dic["FDF1QUADFLAG"] == 1 and dic["FDTRANSPOSED"] == 1:
            multi = 1.0
        elif dic["FDF2QUADFLAG"] == 1 and dic["FDTRANSPOSED"] == 0:
            multi = 1.0
        else:
            multi = 2.0

        dim1 = int(dic["FDSIZE"]*multi)
        dim2 = int(dic["FDSPECNUM"])

        # when the direct dim is singular and the indirect 
        # dim is complex FDSPECNUM is half of the correct value
        if dic["FDQUADFLAG"] == 0 and multi == 1.0:
            dim2 = dim2*2

        
        # check for 3/4D data stream format files (made using xyz2pipe)
        if dic["FDDIMCOUNT"] == 3 and dic["FDPIPEFLAG"] != 0:
            dim3 = int(dic["FDF3SIZE"])
            return (dim3,dim2,dim1)
        if dic["FDDIMCOUNT"] == 4 and dic["FDPIPEFLAG"] != 0:
            dim3 = int(dic["FDF3SIZE"])
            dim4 = int(dic["FDF4SIZE"])
            return (dim4,dim3,dim2,dim1)

        return (dim2,dim1)


def reshape_data(data,shape):
    """ 
    Reshape data or return 1D data after warning
    """
    try:
        return data.reshape(shape)
    except ValueError:
            print "Warning:",data.shape,"cannot be shaped into",shape
            return data


def unshape_data(data):
    """ 
    Returns 1D version of data
    """
    return data.flatten()


def unappend_data(data):
    """ 
    Returns complex data with last axis (-1) unappended

    Data should have imaginary data vector appended to real data vector

    """
    h = data.shape[-1]/2.0
    return np.array(data[...,:h]+data[...,h:]*1.j,dtype="complex64")


def append_data(data):
    """ Return data with last axis (-1) appeneded

    Data should be complex

    """
    return np.concatenate( (data.real,data.imag) , axis=-1)


###################
# fdata functions #
###################

def fdata2dic(fdata):
    """ 
    Convert a fdata array to fdata dictionary
    
    Converts the raw 512x4-byte NMRPipe header into a python dictionary
    with keys as given in fdatap.h
    """

    dic = dict()

    # Populate the dictionary with FDATA which contains numbers
    for key in fdata_dic.keys():
        dic[key] = float(fdata[ int( fdata_dic[key] ) ])
    
    # make the FDDIMORDER
    dic["FDDIMORDER"] = [dic["FDDIMORDER1"],dic["FDDIMORDER2"],  \
                         dic["FDDIMORDER3"],dic["FDDIMORDER4"]]
    
    # Populate the dictionary with FDATA which contains strings
    dic["FDF2LABEL"] = struct.unpack('8s',fdata[16:18])[0].rstrip('\x00')
    dic["FDF1LABEL"] = struct.unpack('8s',fdata[18:20])[0].rstrip('\x00')
    dic["FDF3LABEL"] = struct.unpack('8s',fdata[20:22])[0].rstrip('\x00')
    dic["FDF4LABEL"] = struct.unpack('8s',fdata[22:24])[0].rstrip('\x00')
    dic["FDSRCNAME"]  = struct.unpack('16s' ,fdata[286:290])[0].rstrip('\x00')
    dic["FDUSERNAME"] = struct.unpack('16s' ,fdata[290:294])[0].rstrip('\x00')
    dic["FDTITLE"]    = struct.unpack('60s' ,fdata[297:312])[0].rstrip('\x00')
    dic["FDCOMMENT"]  = struct.unpack('160s',fdata[312:352])[0].rstrip('\x00')
    dic["FDOPERNAME"] = struct.unpack('32s' ,fdata[464:472])[0].rstrip('\x00')

    return dic


def dic2fdata(dic):
    """ 
    Converts a NMRPipe dictionary into an array
    """

    # A 512 4-byte array to hold the nmrPipe header data
    fdata = np.zeros(512,'float32')

    # Populate the array with the simple numbers
    for key in fdata_nums.keys():
        fdata[ int( fdata_dic[key])] = float(dic[key])

    # Check that FDDIMORDER didn't overwrite FDDIMORDER1
    fdata[ int(fdata_dic["FDDIMORDER1"]) ] = dic["FDDIMORDER1"] 

    # Pack the various strings into terminated strings of the correct length
    # then into floats in the fdata array
    fdata[16:18] = struct.unpack('2f', struct.pack('8s',dic["FDF2LABEL"]) )
    fdata[18:20] = struct.unpack('2f', struct.pack('8s',dic["FDF1LABEL"]) )
    fdata[20:22] = struct.unpack('2f', struct.pack('8s',dic["FDF3LABEL"]) )
    fdata[22:24] = struct.unpack('2f', struct.pack('8s',dic["FDF4LABEL"]) )

    # and the longer strings (typically blank)
    fdata[286:290]=struct.unpack( '4f', struct.pack( '16s',dic["FDSRCNAME"]))
    fdata[290:294]=struct.unpack( '4f', struct.pack( '16s',dic["FDUSERNAME"]))
    fdata[297:312]=struct.unpack('15f', struct.pack( '60s',dic["FDTITLE"]))
    fdata[312:352]=struct.unpack('40f', struct.pack('160s',dic["FDCOMMENT"]))
    fdata[464:472]=struct.unpack( '8f', struct.pack( '32s',dic["FDOPERNAME"]))

    return fdata


#################################
# raw reading of data from file #
#################################

def get_fdata(filename):
    """
    Get an array of length 512 holding NMRPipe header
    """
    fdata =  np.fromfile(filename,'float32',512)

    if fdata[2] - 2.345 > 1e-6:    # fdata[2] should be 2.345
        fdata = fdata.byteswap()
    
    return fdata


def get_data(filename):
    """
    Get array of data
    """

    data = np.fromfile(filename,'float32')

    if data[2] - 2.345 > 1e-6:  # check for byteswap
        data = data.byteswap()

    return data[512:]


def get_fdata_data(filename):
    """ 
    Get fdata and data array, returns (fdata,data)
    """
    data = np.fromfile(filename,'float32')
    if data[2] - 2.345 > 1e-6:  # check for byteswap
        data = data.byteswap()

    return data[:512],data[512:]


##############################################
# low memory numpy.ndarray emulating objects #
##############################################

def get_trace(fhandle,ntrace,pts,bswap,cplex):
    """
    Get a single trace from a NMRPipe file 

    Parameters:

    * fhandle   File object of open NMRPipe file.
    * ntrace    Trace numbers (starting from 0).
    * pts       Number of points in trace, R|I.
    * bswap     True/False to perform byteswap on trace.
    * cplex     True/False to unappend imaginary data.

    """

    if cplex:
        tpts = pts*2    # read twice as many points if data is complex
    else:
        tpts = pts
    
    fhandle.seek( 4*(512+ntrace*tpts) ) # seek to the start of the trace 
    trace = np.fromfile(fhandle,'float32',tpts)

    if bswap:
        trace = trace.byteswap()
    if cplex:
      return unappend_data(trace) 
    else:
        return trace


class pipe_2d(fileiobase.data_nd):
    """
    Emulate a numpy.ndarray objects without loading data into memory for low 
    memory reading of 2D NMRPipe files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes 
      ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order=(0,1)):
        """
        Create and set up object
        """

        # read and parse the NMRPipe header
        fdata = get_fdata(filename) # get the header data
        if fdata[2] - 2.345 > 1e-6: # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False
        
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # set object attributes 
        self.filename = filename
        self.order = order

        # check last axis quadrature
        fn = "FDF"+str(int(dic["FDDIMORDER1"]))
        if dic[fn+"QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype=np.dtype('float32')
        else:
            self.cplex = True
            self.dtype=np.dtype('complex64')
            fshape[1] = fshape[1]/2

        # finalize
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes
        
    def __fcopy__(self,order):
        """
        Create a copy
        """
        n = pipe_2d(self.filename,order)
        return n

    def __fgetitem__(self,(sY,sX)):
        """
        Return ndarray of selected values

        (sY,sX) is a well formated tuple of slices
        """
        f = open(self.filename,'r') # open the file for reading
        
        # determine which objects should be selected
        lenY,lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]

        # create an empty array to store the selected slice
        out=np.empty( (len(ych),len(xch)) ,dtype=self.dtype)

        # read in the data trace by trace
        for yi,y in enumerate(ych):
            ntrace = y
            trace = get_trace(f,ntrace,lenX,self.bswap,self.cplex)
            out[yi] = trace[sX] 
        f.close()
        return out


# There are two types of NMRPipe 3D files: 
# 1) streams which are single file data sets made with xyz2pipe.
# 2) multiple file data test, names test%03d.ft3, etc.
# Low memory objects exist for both, choose the correct one, or let read
# do it for you.

class pipe_3d(fileiobase.data_nd):
    """
    Emulate a numpy.ndarray objects without loading data into memory for low
    memory reading of 3D NMRPipe files (multiple file data sets)

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.
    
    """

    def __init__(self,filemask,order=(0,1,2),fcheck=False):    
        """
        Create and set up object, check that files exist if fcheck is True
        """
        filename = filemask%1

        # read and parse the NMRPipe header in the first file of the 3D
        fdata = get_fdata(filename) # get the header data
        if fdata[2] - 2.345 > 1e-6: # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        # find the shape of the first two dimensions
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))[-2:]
        
        # find the length of the third dimension
        f3 = "FDF"+str(int(dic["FDDIMORDER3"]))
        lenZ = int(dic[f3+"SIZE"])
        fshape.insert(0,lenZ)   # insert as leading size of fshape
        
        # check that all files exist if fcheck is set
        if fcheck:
            for i in range(1,lenZ+1):
                if os.path.exists(filemask%i)==False:
                    raise IOError("File not found: "+str(filemask%i))

        # check last axis quadrature
        fn = "FDF"+str(int(dic["FDDIMORDER1"]))
        if dic[fn+"QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype=np.dtype('float32')
        else:
            self.cplex = True
            self.dtype=np.dtype('complex64')
            fshape[2] = fshape[2]/2
        
        # finalize
        self.filemask = filemask
        self.order = order
        self.fshape = fshape
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self,order):
        """ 
        Create a copy 
        """
        n = pipe_3d(self.filemask,order)
        return n

    def __fgetitem__(self, (sZ,sY,sX) ):
        """ 
        Return ndarray of selected values

        (sZ,sY,sX) is a well formated tuple of slices
        """
        # determine which objects should be selected
        lenZ,lenY,lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]

        # create an empty array to store the selected slice
        out=np.empty( (len(zch),len(ych),len(xch)) ,dtype=self.dtype)

        # read in the data file by file and trace by trace
        for zi,z in enumerate(zch):
            f = open(self.filemask%(z+1),'r')   # open the current Z axis file
            for yi,y in enumerate(ych):
                ntrace = y
                trace = get_trace(f,ntrace,lenX,self.bswap,self.cplex)
                out[zi,yi] = trace[sX]
            f.close()
        return out       


class pipestream_3d(fileiobase.data_nd):
    """
    Emulate a numpy.ndarray objects without loading data into memory for low 
    memory reading of 3D NMRPipe data stream (one file 3D):

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes 
      ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order=(0,1,2)):
        """
        Create and set up object
        """

        # read and parse the NMRPipe header
        fdata = get_fdata(filename) # get the header data
        if fdata[2] - 2.345 > 1e-6: # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False
        
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # check last axis quadrature
        fn = "FDF"+str(int(dic["FDDIMORDER1"]))
        if dic[fn+"QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype=np.dtype('float32')
        else:
            self.cplex = True
            self.dtype=np.dtype('complex64')
            fshape[2] = fshape[2]/2

        # finalize
        self.filename = filename
        self.order = order
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes
        
    def __fcopy__(self,order):
        """
        Create a copy
        """
        n = pipestream_3d(self.filename,order)
        return n

    def __fgetitem__(self,(sZ,sY,sX)):
        """
        Return ndarray of selected values

        (sZ,sY,sX) is a well formated tuple of slices
        """
        f = open(self.filename,'r') # open the file for reading
        
        # determine which objects should be selected
        lenZ,lenY,lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]

        # create an empty array to store the selected slice
        out=np.empty( (len(zch),len(ych),len(xch)) ,dtype=self.dtype)

        # read in the data trace by trace
        for zi,z in enumerate(zch):
            for yi,y in enumerate(ych):
                ntrace = y+z*lenY
                trace = get_trace(f,ntrace,lenX,self.bswap,self.cplex)
                out[zi,yi] = trace[sX]
        
        f.close()
        return out

# There are three types of NMRPipe 4D files:
# 1) streams which are single file data sets made with xyz2pipe.
# 2) single index multiple file data sets, named test%03d.ft4, etc.
# 3) two index muttuple file data sets, named test%02d%03d.ft2, made with
# pipe2xyz and conversion binary.
# Low memory objects exist for all three, choose the correct one, or let read
# do it for you.


class pipe_4d(fileiobase.data_nd):
    """
    Emulate a numpy.ndarray objects without loading data into memory for low
    memory reading of single/two index 4D NMRPipe files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.
    
    """

    def __init__(self,filemask,order=(0,1,2,3),fcheck=False):    
        """
        Create and set up object, check that files exist if fcheck is True
        """

        if filemask.count("%")==1:
            self.singleindex=True
            filename = filemask%(1)
        elif filemask.count("%")==2:
            self.singleindex=False
            filename = filemask%(1,1)
        else:
            raise ValueError("bad filemask")

        # read and parse the NMRPipe header in the first file of the 3D
        fdata = get_fdata(filename) # get the header data
        if fdata[2] - 2.345 > 1e-6: # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False

        # find the shape of the first two dimensions
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))[-2:]
        
        # find the length of the third dimension
        f3 = "FDF"+str(int(dic["FDDIMORDER3"]))
        lenZ = int(dic[f3+"SIZE"])
        fshape.insert(0,lenZ)   # insert as leading size of fshape
        
        # find the length of the fourth dimension
        f4 = "FDF"+str(int(dic["FDDIMORDER4"]))
        lenA = int(dic[f4+"SIZE"])
        fshape.insert(0,lenA)   # insert as leading size of fshape


        # check that all files exist if fcheck is set
        if fcheck:
            for ai in range(1,lenA+1):
                for zi in range(1,lenZ+1):
                    if self.singleindex:
                        fname = filemask%(a*lenZ+z+1)
                    else:
                        fname = filemask%(a+1,z+1)

                    if os.path.exists(fname)==False:
                        raise IOError("File not found: "+str(fname))

        # check last axis quadrature
        fn = "FDF"+str(int(dic["FDDIMORDER1"]))
        if dic[fn+"QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype=np.dtype('float32')
        else:
            self.cplex = True
            self.dtype=np.dtype('complex64')
            fshape[3] = fshape[3]/2
        
        # finalize
        self.filemask = filemask
        self.order = order
        self.fshape = fshape
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self,order):
        """ 
        Create a copy 
        """
        n = pipe_4d(self.filemask,order)
        return n

    def __fgetitem__(self, (sA,sZ,sY,sX) ):
        """ 
        Return ndarray of selected values

        (sZ,sY,sX) is a well formated tuple of slices
        """
        # determine which objects should be selected
        lenA,lenZ,lenY,lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]
        ach = range(lenA)[sA]

        # create an empty array to store the selected slice
        out=np.empty( (len(ach),len(zch),len(ych),len(xch)) ,dtype=self.dtype)

        # read in the data file by file, trace by trace
        for ai,a in enumerate(ach):
            for zi,z in enumerate(zch):                
                if self.singleindex:   # single index
                    f = open(self.filemask%(a*lenZ+z+1),'r')
                else:   # two index
                    f = open(self.filemask%(a+1,z+1),'r')
                for yi,y in enumerate(ych):
                    ntrace = y
                    trace = get_trace(f,ntrace,lenX,self.bswap,self.cplex)
                    out[ai,zi,yi] = trace[sX]
                f.close()
        return out


class pipestream_4d(fileiobase.data_nd):
    """
    Emulate a numpy.ndarray objects without loading data into memory for low 
    memory reading of 4D NMRPipe data stream (one file 4D):

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes 
      ordering.
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,filename,order=(0,1,2,3)):
        """
        Create and set up object
        """

        # read and parse the NMRPipe header
        fdata = get_fdata(filename) # get the header data
        if fdata[2] - 2.345 > 1e-6: # check if byteswapping will be necessary
            self.bswap = True
        else:
            self.bswap = False
        
        dic = fdata2dic(fdata)  # create the dictionary
        fshape = list(find_shape(dic))

        # set object attributes 
        self.filename = filename
        self.order = order

        # check last axis quadrature
        fn = "FDF"+str(int(dic["FDDIMORDER1"]))
        if dic[fn+"QUADFLAG"] == 1.0:
            self.cplex = False
            self.dtype=np.dtype('float32')
        else:
            self.cplex = True
            self.dtype=np.dtype('complex64')
            fshape[3] = fshape[3]/2

        # finalize
        self.fshape = tuple(fshape)
        self.__setdimandshape__()   # set ndim and shape attributes
        
    def __fcopy__(self,order):
        """
        Create a copy
        """
        n = pipestream_4d(self.filename,order)
        return n

    def __fgetitem__(self,(sA,sZ,sY,sX)):
        """
        Return ndarray of selected values

        (sA,sZ,sY,sX) is a well formated tuple of slices
        """
        f = open(self.filename,'r') # open the file for reading
        
        # determine which objects should be selected
        lenA,lenZ,lenY,lenX = self.fshape
        xch = range(lenX)[sX]
        ych = range(lenY)[sY]
        zch = range(lenZ)[sZ]
        ach = range(lenA)[sA]

        # create an empty array to store the selected slice
        out=np.empty( (len(ach),len(zch),len(ych),len(xch)) ,dtype=self.dtype)

        # read in the data trace by trace
        for ai,a in enumerate(ach):
            for zi,z in enumerate(zch):
                for yi,y in enumerate(ych):
                    ntrace = y+z*lenY+a*lenY*lenZ
                    trace = get_trace(f,ntrace,lenX,self.bswap,self.cplex)
                    out[ai,zi,yi] = trace[sX]
        
        f.close()
        return out




# data, see fdata.h

fdata_nums = {
'FDF4CENTER': '82', 'FDF2P0': '109', 'FDF2P1': '110', 'FDF1P1': '246',
'FDF2X1': '257', 'FDF1P0': '245', 'FDF3AQSIGN': '476', 'FDDISPMAX': '251',
'FDF4FTFLAG': '31', 'FDF3X1': '261', 'FDRANK': '180', 'FDF2C1': '418',
'FDF2QUADFLAG': '56', 'FDSLICECOUNT': '443', 'FDFILECOUNT': '442',
'FDMIN': '248', 'FDF3OBS': '10', 'FDF4APODQ2': '407', 'FDF4APODQ1': '406', 
'FDF3FTSIZE': '200', 'FDF1LB': '243', 'FDF4C1': '409', 'FDF4QUADFLAG': '54',
'FDF1SW': '229', 'FDTRANSPOSED': '221', 'FDSECS': '285', 'FDF1APOD': '428',
'FDF2APODCODE': '413', 'FDPIPECOUNT': '75',
'FDPEAKBLOCK': '362', 'FDREALSIZE': '97', 'FDF4SIZE': '32',
'FDF4SW': '29', 'FDF4ORIG': '30', 'FDF3XN': '262', 'FDF1OBS': '218',
'FDDISPMIN': '252', 'FDF2XN': '258', 'FDF3P1': '61', 'FDF3P0': '60',
'FDF1ORIG': '249', 'FDF2FTFLAG': '220', 'FDF1TDSIZE': '387', 'FDLASTPLANE': '78',
'FDF1ZF': '437', 'FDF4FTSIZE': '201', 'FDF3C1': '404', 'FDFLTFORMAT': '1',
'FDF4CAR': '69', 'FDF1FTFLAG': '222', 'FDF2OFFPPM': '480',
'FDSIZE': '99', 'FDYEAR': '296', 'FDF1C1': '423', 'FDUSER3': '72',
'FDF1FTSIZE': '98', 'FDMINS': '284', 'FDSCALEFLAG': '250', 'FDF3TDSIZE': '388',
'FDPARTITION': '65', 'FDF3FTFLAG': '13', 'FDF2APODQ1': '415',
'FD2DVIRGIN': '399', 'FDF2APODQ3': '417', 'FDF2APODQ2': '416',
'FD2DPHASE': '256', 'FDMAX': '247', 'FDF3SW': '11', 'FDF4TDSIZE': '389',
'FDPIPEFLAG': '57', 'FDDAY': '295', 'FDF2UNITS': '152', 'FDF4APODQ3': '408',
'FDFIRSTPLANE': '77', 'FDF3SIZE': '15', 'FDF3ZF': '438',
'FDF3ORIG': '12', 'FD1DBLOCK': '365', 'FDF1AQSIGN': '475', 'FDF2OBS': '119',
'FDF1XN': '260', 'FDF4UNITS': '59', 'FDDIMCOUNT': '9', 'FDF4XN': '264',
'FDUSER2': '71', 'FDF4APODCODE': '405', 'FDUSER1': '70', 'FDMCFLAG': '135',
'FDFLTORDER': '2', 'FDUSER5': '74', 'FDF3QUADFLAG': '51',
'FDUSER4': '73', 'FDTEMPERATURE': '157', 'FDF2APOD': '95', 'FDMONTH': '294',
'FDF4OFFPPM': '483', 'FDF3OFFPPM': '482', 'FDF3CAR': '68', 'FDF4P0': '62', 
'FDF4P1': '63', 'FDF1OFFPPM': '481', 'FDF4APOD': '53', 'FDF4X1': '263',
'FDLASTBLOCK': '359', 'FDPLANELOC': '14', 'FDF2FTSIZE': '96',
'FDF1X1': '259', 'FDF3CENTER': '81', 'FDF1CAR': '67', 'FDMAGIC': '0', 
'FDF2ORIG': '101', 'FDSPECNUM': '219', 'FDF2AQSIGN': '64',
'FDF1UNITS': '234', 'FDF2LB': '111', 'FDF4AQSIGN': '477', 'FDF4ZF': '439',
'FDTAU': '199', 'FDNOISE': '153', 'FDF3APOD': '50',
'FDF1APODCODE': '414', 'FDF2SW': '100', 'FDF4OBS': '28', 'FDQUADFLAG': '106',
'FDF2TDSIZE': '386', 'FDHISTBLOCK': '364', 
'FDBASEBLOCK': '361', 'FDF1APODQ2': '421', 'FDF1APODQ3': '422',
'FDF1APODQ1': '420', 'FDF1QUADFLAG': '55', 'FDF3UNITS': '58', 'FDF2ZF': '108',
'FDCONTBLOCK': '360', 'FDDIMORDER4': '27', 'FDDIMORDER3': '26', 
'FDDIMORDER2': '25', 'FDDIMORDER1': '24', 'FDF2CAR': '66', 'FDF3APODCODE': '400',
'FDHOURS': '283', 'FDF1CENTER': '80', 'FDF3APODQ1': '401', 'FDF3APODQ2': '402',
'FDF3APODQ3': '403', 'FDBMAPBLOCK': '363', 'FDF2CENTER': '79'}

fdata_dic = {
'FDF4CENTER': '82', 'FDF2P0': '109', 'FDF2P1': '110', 'FDF1P1': '246',
'FDF2X1': '257', 'FDF1P0': '245', 'FDF3AQSIGN': '476', 'FDDISPMAX': '251',
'FDF4FTFLAG': '31', 'FDF3X1': '261', 'FDRANK': '180', 'FDF2C1': '418',
'FDF2QUADFLAG': '56', 'FDSLICECOUNT': '443', 'FDFILECOUNT': '442',
'FDMIN': '248', 'FDF3OBS': '10', 'FDF4APODQ2': '407', 'FDF4APODQ1': '406', 
'FDF3FTSIZE': '200', 'FDF1LB': '243', 'FDF4C1': '409', 'FDF4QUADFLAG': '54',
'FDF1SW': '229', 'FDTRANSPOSED': '221', 'FDSECS': '285', 'FDF1APOD': '428',
'FDF2APODCODE': '413', 'FDPIPECOUNT': '75', 'FDOPERNAME': '464',
'FDF3LABEL': '20', 'FDPEAKBLOCK': '362', 'FDREALSIZE': '97', 'FDF4SIZE': '32',
'FDF4SW': '29', 'FDF4ORIG': '30', 'FDF3XN': '262', 'FDF1OBS': '218',
'FDDISPMIN': '252', 'FDF2XN': '258', 'FDF3P1': '61', 'FDF3P0': '60',
'FDF1ORIG': '249', 'FDF2FTFLAG': '220', 'FDF1TDSIZE': '387', 'FDLASTPLANE': '78',
'FDF1ZF': '437', 'FDF4FTSIZE': '201', 'FDF3C1': '404', 'FDFLTFORMAT': '1',
'FDF4CAR': '69', 'FDF1FTFLAG': '222', 'FDF2OFFPPM': '480', 'FDF1LABEL': '18',
'FDSIZE': '99', 'FDYEAR': '296', 'FDF1C1': '423', 'FDUSER3': '72',
'FDF1FTSIZE': '98', 'FDMINS': '284', 'FDSCALEFLAG': '250', 'FDF3TDSIZE': '388',
'FDTITLE': '297', 'FDPARTITION': '65', 'FDF3FTFLAG': '13', 'FDF2APODQ1': '415',
'FD2DVIRGIN': '399', 'FDF2APODQ3': '417', 'FDF2APODQ2': '416',
'FD2DPHASE': '256', 'FDMAX': '247', 'FDF3SW': '11', 'FDF4TDSIZE': '389',
'FDPIPEFLAG': '57', 'FDDAY': '295', 'FDF2UNITS': '152', 'FDF4APODQ3': '408',
'FDFIRSTPLANE': '77', 'FDF3SIZE': '15', 'FDF3ZF': '438', 'FDDIMORDER': '24',
'FDF3ORIG': '12', 'FD1DBLOCK': '365', 'FDF1AQSIGN': '475', 'FDF2OBS': '119',
'FDF1XN': '260', 'FDF4UNITS': '59', 'FDDIMCOUNT': '9', 'FDF4XN': '264',
'FDUSER2': '71', 'FDF4APODCODE': '405', 'FDUSER1': '70', 'FDMCFLAG': '135',
'FDFLTORDER': '2', 'FDUSER5': '74', 'FDCOMMENT': '312', 'FDF3QUADFLAG': '51',
'FDUSER4': '73', 'FDTEMPERATURE': '157', 'FDF2APOD': '95', 'FDMONTH': '294',
'FDF4OFFPPM': '483', 'FDF3OFFPPM': '482', 'FDF3CAR': '68', 'FDF4P0': '62', 
'FDF4P1': '63', 'FDF1OFFPPM': '481', 'FDF4APOD': '53', 'FDF4X1': '263',
'FDLASTBLOCK': '359', 'FDPLANELOC': '14', 'FDF2FTSIZE': '96', 'FDUSERNAME': '290',
'FDF1X1': '259', 'FDF3CENTER': '81', 'FDF1CAR': '67', 'FDMAGIC': '0', 
'FDF2ORIG': '101', 'FDSPECNUM': '219', 'FDF2LABEL': '16', 'FDF2AQSIGN': '64',
'FDF1UNITS': '234', 'FDF2LB': '111', 'FDF4AQSIGN': '477', 'FDF4ZF': '439',
'FDTAU': '199', 'FDF4LABEL': '22', 'FDNOISE': '153', 'FDF3APOD': '50',
'FDF1APODCODE': '414', 'FDF2SW': '100', 'FDF4OBS': '28', 'FDQUADFLAG': '106',
'FDF2TDSIZE': '386', 'FDHISTBLOCK': '364', 'FDSRCNAME': '286', 
'FDBASEBLOCK': '361', 'FDF1APODQ2': '421', 'FDF1APODQ3': '422',
'FDF1APODQ1': '420', 'FDF1QUADFLAG': '55', 'FDF3UNITS': '58', 'FDF2ZF': '108',
'FDCONTBLOCK': '360', 'FDDIMORDER4': '27', 'FDDIMORDER3': '26', 
'FDDIMORDER2': '25', 'FDDIMORDER1': '24', 'FDF2CAR': '66', 'FDF3APODCODE': '400',
'FDHOURS': '283', 'FDF1CENTER': '80', 'FDF3APODQ1': '401', 'FDF3APODQ2': '402',
'FDF3APODQ3': '403', 'FDBMAPBLOCK': '363', 'FDF2CENTER': '79'}
