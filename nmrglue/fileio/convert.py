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

        vdic,vdata = ng.varian.read("varian_dir")
        C = ng.convert.converter()
        C.from_varian(vdic,vdata)
        pdic,pdata = C.to_pipe()
        ng.pipe.write("test.fid",pdic,pdata)

    """
    def __init__(self):
        """ 
        Create a converter object
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
        if isinstance(self._data,np.ndarray):   # in memory data
            return self.__procdata()

        else:   # return emulated data
            iproc = self._iproc
            oproc = self._oproc
            odtype = self._odtype
            order = self._data.order
            return udata_nd(self._data,iproc,oproc,odtype,order)

    def __procdata(self):
        """ 
        Process data as indicated by flags
        """

        # copy the data
        data = np.copy(self._data)

        # processing for input type
        # sign alt. indirect dimension
        if data.ndim >= 2 and "alt_id_sign" in self._iproc:
            #data[1::2] = -data[1::2]
            s = [slice(None,None,None)]*data.ndim
            for i in range(data.ndim-1):
                s[i] = slice(1,None,2)
                data[s] = -data[s]
                s[i] = slice(None,None,None)

        if "realfactor" in self._iproc:
            data.real = data.real*self._iproc['realfactor']
            
        if "imagfactor" in self._iproc:
            data.imag = data.imag*self._iproc['imagfactor']

        # processing for output 
        # sign alt. indirect dimension
        if data.ndim >= 2 and "alt_id_sign" in self._oproc:
            s = [slice(None,None,None)]*data.ndim
            for i in range(data.ndim-1):
                s[i] = slice(1,None,2)
                data[s] = -data[s]
                s[i] = slice(None,None,None)

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

class udata_nd(fileiobase.data_nd):
    """
    Wrap other fileiobase.data_nd derived objects with input/output conversion
    when slices are requested.  

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    """
    
    def __init__(self,edata,iproc,oproc,odtype,order=None):
        """
        create and set up
        """
        # set converter attributes
        self._iproc = iproc         # input processing dictionary
        self._oproc = oproc         # output processing dictionary
        self._odtype = odtype       # output dtype
        self.edata = edata          # file
        
        # required data_nd attributes
        self.order = order  
        self.fshape = edata.fshape
        self.dtype = odtype
        self.__setdimandshape__()   # set ndim and shape attributes

    def __fcopy__(self,order):
        """
        Create a copy
        """
        n = udata_nd(self.edata,self._iproc,self._oproc,self._odtype,order)
        return n

    def __fgetitem__(self,slices):
        """
        Return ndarray of selected values

        slices is a well formateed n-tuple of slices
        """
        data = self.edata.__fgetitem__(slices)
        
        # input processing 
        if "alt_id_sign" in self._iproc:    # sign alt. indirect dimension
            if "alt_id_sign" not in self._oproc:    # skip if in both
                fslice = slices[:-1]
                ffshape = self.fshape[:-1]
                nd_iter = fileiobase.ndtofrom_iter(ffshape,fslice)
                for out_index,in_index in nd_iter:
                    # negate the trace if there is an odd number of 
                    # odd number indices in the slice
                    if np.mod(in_index,2).sum()%2 == 1: 
                        data[out_index] = -data[out_index]


        if "realfactor" in self._iproc:
            data.real = data.real*self._iproc['realfactor']
        if "imagfactor" in self._iproc:
            data.imag = data.imag*self._iproc['imagfactor']

        # output processing
        if "alt_id_sign" in self._oproc:
            if "alt_id_sign" not in self._iproc:
                fslice = slices[:-1]
                ffshape = self.fshape[:-1]
                nd_iter = fileiobase.ndtofrom_iter(ffshape,fslice)
                for out_index,in_index in nd_iter:
                    # negate the trace if there is an odd number of 
                    # odd number indices in the slice
                    if np.mod(in_index,2).sum()%2 == 1: 
                        data[out_index] = -data[out_index]

           

        if "realfactor" in self._oproc:
            data.real = data.real*self._oproc['realfactor']
        if "imagfactor" in self._oproc:
            data.imag = data.imag*self._oproc['imagfactor']
        
        return data.astype(self._odtype)
