"""
Functions to convert between NMR file formats

"""

import numpy as np
import datetime

# modules from nmrglue
import pipe
import varian
import bruker
import sparky
import fileiobase

class converter(object):
    """
    The converter object allowes conversion between NMR file formats

    Example::

        vdic,vdata = ng.varian.read("fid")
        C = ng.convert.converter()
        C.from_varian(vdic,vdata)
        pdic,pdata = C.to_pipe()
        ng.pipe.write("test.fid",pdic,pdata)

    """


    def __init__(self):
        """ 
        Create and set up the object 
        """

        pass

    # utility functions

    def __returndata(self):
        """ 
        Return data or emulated data after error checking
        """

        # Error checking
        if "_data" not in self.__dict__:
            raise IOError("converter not loaded with data")
        if "_udic" not in self.__dict__:
            raise IOError("converter not loaded with dictionary")
        if "_iproc" not in self.__dict__:
            raise IOError("converter not loaded with processing parameters")
        if "_oproc" not in self.__dict__:
            raise IOError("converter not loaded with processing parameters")
        if "_odtype" not in self.__dict__:
            raise IOError("converter not loaded output dtype")

        # Warnings
        if self._data.dtype.kind != np.dtype(self._odtype).kind:
            print "Warning: Incompatiable dtypes, conversion not recommended"


        # Return data
        if isinstance(self._data,np.ndarray):   # 'real' data
            return self.__procdata()

        else:   # return emulated data
            # create the needed emulated data class
            if self._data.ndim == 2:
                ptup = (self._iproc,self._oproc,self._odtype)
                return udata_2d(self._data,ptup,self._data.order)
            if self._data.ndim == 3:
                ptup = (self._iproc,self._oproc,self._odtype)
                return udata_3d(self._data,ptup,self._data.order)
            else:
                raise NotImplementedError("dimension of data not supported")


    def __procdata(self):
        """ 
        Process data as indicated by flags
        """

        # copy the data
        data = np.copy(self._data)

        # processing for input type
        # sign alt. indirect dimension
        if data.ndim >= 2 and "alt_id_sign" in self._iproc:
            data[1::2] = -data[1::2]

        if "realfactor" in self._iproc:
            data.real = data.real*self._iproc['realfactor']
            
        if "imagfactor" in self._iproc:
            data.imag = data.imag*self._iproc['imagfactor']

        # processing for output 
        # sign alt. indirect dimension
        if data.ndim >= 2 and "alt_id_sign" in self._oproc:
            data[1::2] = -data[1::2]

        if "realfactor" in self._oproc:
            data.real = data.real*self._oproc['realfactor']
            
        if "imagfactor" in self._oproc:
            data.imag = data.imag*self._oproc['imagfactor']

        return data.astype(self._odtype)


    # IMPORTERS (from_*)

    def from_universal(self,dic,data):
        """ 
        load data from universal dic,data pair
        """

        # set data
        self._data = data
        self._iproc = {}

        # set the dictionary
        self._udic = dic

    def from_varian(self,dic,data,udic=None):
        """ 
        load data and dictionary from varian fid pair

        Parameter:

        * dic     Varian dictionary
        * data    Varian data
        * udic    Universal dictionary. If not provided is guessed.

        """

        # set data
        self._data = data
        if udic != None and udic[0]['encoding'].lower() == "tppi":
            self._iproc = {"imagfactor":-1.0}
        else:   # states, etc needs sign alt. of indirect dim.
            self._iproc = {"alt_id_sign":True,"imagfactor":-1.0}

        # set the universal dictionary
        if udic != None:
            self._udic = udic
        else:
            self._udic = varian.guess_udic(dic,data) 


    def from_pipe(self,dic,data,udic=None):
        """ 
        Load data and dictionary from NMRPipe pair
        """

        # set data
        self._data = data
        self._iproc = {}
       
        # set the universal dictionary
        if udic != None:
            self._udic = udic
        else:
            self._udic = pipe.guess_udic(dic,data)


    def from_sparky(self,dic,data,udic=None):
        """ 
        Load data and dictionary from Sparky pair
        """

        # set data
        self._data = data
        self._iproc = {}


        # set the universal dictionary
        if udic != None:
            self._udic = udic
        else:
            self._udic = sparky.guess_udic(dic,data)


    def from_bruker(self,dic,data,udic=None):
        """ 
        Load data and dictionary from Bruker ser/fid file
        """

        # set data
        self._data = data
        self._iproc = {}

        # set the universal dictionary
        if udic != None:
            self._udic = udic
        else:
            self._udic = bruker.guess_udic(dic,data)

    

    # EXPORTERS (to_*)

    def to_universal(self):
        """ 
        Return universal dictionary and original data
        """

        # create dictionary
        dic = dict(self._udic)

        # add processing flags for output
        self._oproc = {}
        self._odtype = self._data.dtype

        return dic,self.__returndata()


    def to_pipe(self,datetimeobj=datetime.datetime.now()):
        """ 
        Return NMRPipe dic,data pair
        """
           
        # create dictionary
        dic = pipe.create_dic(self._udic,datetimeobj)

        # add processing flags for output
        self._oproc = {}
        if self._udic[self._udic["ndim"]-1]["complex"]:
            self._odtype = "complex64"
        else:
            self._odtype = "float32"

        return dic,self.__returndata()


    def to_varian(self):
        """ 
        Return Varian dic,data pair
        """

        # create dictionary
        dic = varian.create_dic(self._udic)

        # add processing flags for output
        self._oproc = {"alt_id_sign":True,"imagfactor":-1.0}
        self._odtype = "complex64"

        return dic,self.__returndata()


    def to_sparky(self,datetimeobj=datetime.datetime.now(),user='user'):
        """ 
        Return sparky dic,data pair
        """
        
        # create dictionary
        dic = sparky.create_dic(self._udic,datetimeobj,user)

        # add processing flags for output
        self._oproc = {}
        self._odtype = "float32"

        return dic,self.__returndata()


    def to_bruker(self):
        """ 
        Return Bruker dic,data pair
        """
    
        # create dictionary
        dic = bruker.create_dic(self._udic)

        # add processing flags for output
        self._oproc = {}
        self._odtype = "complex128"

        return dic,self.__returndata()


class udata_2d(fileiobase.data_2d):
    """
    udata_2d emulates a numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects
    * can be used in iteration with expected results
    * transpose and swapaxes functions create a new udata_2d object with
      new axes ordering.
    * has ndim,shape and dtype attributes

    This object processes other data_2d derived objects as they are read
    from disk.

    """

    def __init__(self,data,(iproc,oproc,odtype),order=["y","x"]):
        """ 
        Create and set up a udata_2d object
        """
        # copy the converter attributes
        self._iproc  = iproc
        self._oproc  = oproc
        self._odtype = odtype

        # create a copy of the data with correct order
        self.data = data.__fcopy__(order)

        # required data_2d attributes
        self.lenX = self.data.lenX
        self.lenY = self.data.lenY
        self.order = self.data.order
        self.shape = self.data.shape
        self.dtype = self._odtype
        self.ndim = 2

    def __fcopy__(self,order):
        """ 
        Create a copy with given order
        """

        ptup = (self._iproc,self._oproc,self._odtype)
        n = udata_2d(self.data,ptup,order)
        return n

    def __fgetitem__(self,(sY,sX)):
        """ 
        Returns ndarray of selected values

        sY,sX is a well formatted 2-tuple of slices

        """

        # get the raw data
        data = self.data.__fgetitem__( (sY,sX))

        # process the data for output

        # sign alt. indirect dimension
        if "alt_id_sign" in self._iproc:
            # XXX there is probably a better way to do this with tile, etc
            factor = np.ones( (self.lenY,1) )
            factor[1::2]=-1
            data = factor[sY]*data

        if "realfactor" in self._iproc:
            data.real = data.real*self._iproc['realfactor']
            
        if "imagfactor" in self._iproc:
            data.imag = data.imag*self._iproc['imagfactor']

        # processing for output 

        # sign alt. indirect dimension
        if "alt_id_sign" in self._oproc:
            factor = np.ones( (self.lenY,1) )
            factor[1::2]=-1
            data = factor[sY]*data

        if "realfactor" in self._oproc:
            data.real = data.real*self._oproc['realfactor']
            
        if "imagfactor" in self._oproc:
            data.imag = data.imag*self._oproc['imagfactor']

        return data.astype(self._odtype)


class udata_3d(fileiobase.data_3d):
    """
    udata_3d emulates a numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects
    * can be used in iteration with expected results
    * transpose and swapaxes functions create a new udata_3d object with
      new axes ordering.
    * has ndim,shape and dtype attributes

    This object processes other data_3d derived objects as they are read
    from disk.

    """

    def __init__(self,data,(iproc,oproc,odtype),order=["z","y","x"]):
        """ 
        Create and set up a udata_3d object
        """
        # copy the converter attributes
        self._iproc  = iproc
        self._oproc  = oproc
        self._odtype = odtype

        # create a copy of the data with correct order
        self.data = data.__fcopy__(order)

        # required data_2d attributes
        self.lenX = self.data.lenX
        self.lenY = self.data.lenY
        self.lenZ = self.data.lenZ
        self.order = self.data.order
        self.shape = self.data.shape
        self.dtype = self._odtype
        self.ndim = 3

    def __fcopy__(self,order):
        """ 
        Create a copy with given order
        """

        ptup = (self._iproc,self._oproc,self._odtype)
        n = udata_3d(self.data,ptup,order)
        return n

    def __fgetitem__(self,(sZ,sY,sX)):
        """ 
        Returns ndarray of selected values

        (sZ,sY,sX) is a  formatted tuple of slices

        """
        # get the raw data
        data = self.data.__fgetitem__( (sZ,sY,sX))

        # process the data for output

        # sign alt. indirect dimension 
        # XXX there is probably a better way to do this...
        if "alt_id_sign" in self._iproc:
           
            # alternate Y axis 
            factor1 = np.ones( (self.lenY,1) )
            factor1[1::2]=-1
            for i,j in enumerate(factor1[sY]):
                if j==-1:
                    np.atleast_3d(data)[:,i,:]*=-1

            # alternate Z axis
            factor2 = np.ones( (self.lenZ,1) )
            factor2[1::2]=-1
            for i,j in enumerate(factor2[sZ]):
                if j==-1:
                    np.atleast_3d(data)[i,:,:]*=-1

        if "realfactor" in self._iproc:
            data.real = data.real*self._iproc['realfactor']
            
        if "imagfactor" in self._iproc:
            data.imag = data.imag*self._iproc['imagfactor']

        # processing for output 

        # sign alt. indirect dimension
        if "alt_id_sign" in self._oproc:
            factor = np.ones( (self.lenY,1) )
            factor[1::2]=-1
            data = factor[sY]*data

        if "realfactor" in self._oproc:
            data.real = data.real*self._oproc['realfactor']
            
        if "imagfactor" in self._oproc:
            data.imag = data.imag*self._oproc['imagfactor']

        return data.astype(self._odtype)
