"""
Functions to convert between NMR file formats
"""

import datetime
from warnings import warn

import numpy as np

from . import pipe
from . import varian
from . import bruker
from . import sparky
from . import rnmrtk
from . import fileiobase


class converter(object):
    """
    Object which allows conversion between NMR file formats, including low
    memory data objects.

    Conversion between NMR file formats with this class involves three steps.
    First a new converter object must be created.  Then the converter must be
    loaded with data using a ``from_`` method.  Finally, the dictionary and
    data representation of a NMR data in the desired format is extracted using
    a ``to_`` method.  This can then be written to disk.

    Example conversion::

        vdic, vdata = ng.varian.read("varian_dir")
        C = ng.convert.converter()
        C.from_varian(vdic, vdata)
        pdic, pdata = C.to_pipe()
        ng.pipe.write("test.fid", pdic, pdata)

    Spectral parameters can be provided directly by passing a Universal
    dictionary to any of the ``from_`` methods.  If not provided the spectral
    parameters are guessed from the file format's dictionary of parameters.

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
            raise IOError("converter not loaded with output dtype")

        # Warnings
        if self._data.dtype.kind != np.dtype(self._odtype).kind:
            warn("Incompatible dtypes, conversion not recommended")

        # Return data
        if isinstance(self._data, np.ndarray):   # in memory data
            return self.__procdata()

        else:   # return emulated data
            iproc = self._iproc
            oproc = self._oproc
            odtype = self._odtype
            order = self._data.order
            return udata_nd(self._data, iproc, oproc, odtype, order)

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
            s = [slice(None, None, None)] * data.ndim
            for i in range(data.ndim - 1):
                s[i] = slice(1, None, 2)
                data[s] = -data[s]
                s[i] = slice(None, None, None)

        if "realfactor" in self._iproc:
            data.real = data.real * self._iproc['realfactor']

        if "imagfactor" in self._iproc and np.iscomplexobj(data):
            data.imag = data.imag * self._iproc['imagfactor']

        # processing for output
        # sign alt. indirect dimension
        if data.ndim >= 2 and "alt_id_sign" in self._oproc:
            s = [slice(None, None, None)] * data.ndim
            for i in range(data.ndim - 1):
                s[i] = slice(1, None, 2)
                data[s] = -data[s]
                s[i] = slice(None, None, None)

        if "realfactor" in self._oproc:
            data.real = data.real * self._oproc['realfactor']

        if "imagfactor" in self._oproc and np.iscomplexobj(data):
            data.imag = data.imag * self._oproc['imagfactor']

        return data.astype(self._odtype)

    # IMPORTERS (from_*)
    def from_universal(self, dic, data):
        """
        Load converter with Universal data.

        Parameters
        ----------
        dic : dict
            Dictionary of universal parameters.
        data : array_like
            NMR data.

        """
        # set data
        self._data = data
        self._iproc = {}

        # set the dictionary
        self._udic = dic

    def from_varian(self, dic, data, udic=None):
        """
        Load converter with Agilent/Varian data.

        Parameters
        ----------
        dic : dict
            Dictionary of Agilent/Varian parameters.
        data : array_like
            NMR data.
        udic : dict, optional
            Universal dictionary, if not provided will be guesses from dic.

        """
        # set data
        self._data = data
        if udic is not None and udic[0]['encoding'].lower() == "tppi":
            self._iproc = {"imagfactor": -1.0}
        else:   # states, etc needs sign alt. of indirect dim.
            self._iproc = {"alt_id_sign": True, "imagfactor": -1.0}

        # set the universal dictionary
        if udic is not None:
            self._udic = udic
        else:
            self._udic = varian.guess_udic(dic, data)

    def from_rnmrtk(self, dic, data, udic=None, agilent_compatible=False):
        """
        Load converter with RNMRTK data.

        Parameters
        ----------
        dic : dict
            Dictionary of RNMRTK parameters.
        data : array_like
            NMR data.
        udic : dict, optional
            Universal dictionary, if not provided will be guesses from dic.
        agilent_compatible : bool, optional
            True when RNMRTK data is being compared to Agilent/Varian data.

        """
        # set data
        self._data = data

        # set input processing filters.
        if agilent_compatible:
            self._iproc = {"alt_id_sign": True, "imagfactor": -1.0}
        else:
            self._iproc = {}

        # set the universal dictionary
        if udic is not None:
            self._udic = udic
        else:
            self._udic = rnmrtk.guess_udic(dic, data)

    def from_pipe(self, dic, data, udic=None):
        """
        Load converter with NMRPipe data.

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : array_like
            NMR data.
        udic : dict, optional
            Universal dictionary, if not provided will be guesses from dic.

        """
        # set data
        self._data = data
        self._iproc = {}

        # set the universal dictionary
        if udic is not None:
            self._udic = udic
        else:
            self._udic = pipe.guess_udic(dic, data)

    def from_sparky(self, dic, data, udic=None):
        """
        Load converter with Sparky data.

        Parameters
        ----------
        dic : dict
            Dictionary of Sparky parameters.
        data : array_like
            NMR data.
        udic : dict, optional
            Universal dictionary, if not provided will be guesses from dic.

        """
        # set data
        self._data = data
        self._iproc = {}

        # set the universal dictionary
        if udic is not None:
            self._udic = udic
        else:
            self._udic = sparky.guess_udic(dic, data)

    def from_bruker(self, dic, data, udic=None):
        """
        Load converter with Bruker data.

        Parameters
        ----------
        dic : dict
            Dictionary of Bruker parameters.
        data : array_like
            NMR data.
        udic : dict, optional
            Universal dictionary, if not provided will be guesses from dic.

        """
        # set data
        self._data = data
        self._iproc = {}

        # set the universal dictionary
        if udic is not None:
            self._udic = udic
        else:
            self._udic = bruker.guess_udic(dic, data)

    # EXPORTERS (to_*)
    def to_universal(self):
        """
        Return Universal format data.

        Returns
        -------
        dic : dict
            Dictionary of Universal parameters.
        data : array_like
            NMR data in format as provided.

        """
        # create dictionary
        dic = dict(self._udic)

        # add processing flags for output
        self._oproc = {}
        self._odtype = self._data.dtype

        return dic, self.__returndata()

    def to_pipe(self, datetimeobj=datetime.datetime.now()):
        """
        Return NMRPipe format data.

        Parameters
        ----------
        datetime : datetime object, optional
            Datetime object to include in the NMRPipe parameters.  The current
            date and time is used by default.

        Returns
        -------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : array_like
            NMR data in NMRPipe format.

        """
        # create dictionary
        dic = pipe.create_dic(self._udic, datetimeobj)

        # add processing flags for output
        self._oproc = {}
        if self._udic[self._udic["ndim"] - 1]["complex"]:
            self._odtype = "complex64"
        else:
            self._odtype = "float32"

        return dic, self.__returndata()

    def to_rnmrtk(self, agilent_compatible=False, dim_order=None):
        """
        Return RNMRTK format data.

        Parameters
        ----------
        agilent_compatible : bool, optional
            True when RNMRTK data is being compared to Agilent/Varian data.
        dim_order : list, optional
            List mapping axis numbers in the universal dictionary to the to the
            order in which they will appear in the RNMRTK dictionary.  If None,
            the default, [0, 1, 2, ...] will be used.


        Returns
        -------
        dic : dict
            Dictionary of RNMRTK parameters.
        data : array_like
            NMR data in RNMRTK format.

        """
        # create dictionary
        dic = rnmrtk.create_dic(self._udic)

        # add processing flags for output
        if agilent_compatible:
            self._oproc = {"alt_id_sign": True, "imagfactor": -1.0}
        else:
            self._oproc = {}

        if self._udic[self._udic["ndim"] - 1]["complex"]:
            self._odtype = "complex64"
        else:
            self._odtype = "float32"

        return dic, self.__returndata()

    def to_varian(self):
        """
        Return Agilent/Varian format data.

        Returns
        -------
        dic : dict
            Dictionary of Agilent/Varian parameters.
        data : array_like
            NMR data in Agilent/Varian format.

        """
        # create dictionary
        dic = varian.create_dic(self._udic)

        # add processing flags for output
        self._oproc = {"alt_id_sign": True, "imagfactor": -1.0}
        self._odtype = "complex64"

        return dic, self.__returndata()

    def to_sparky(self, datetimeobj=datetime.datetime.now(), user='user'):
        """
        Return Sparky format data.

        Parameters
        ----------
        datetime : datetime object, optional
            Datetime object to include in the Sparky parameters.  The current
            date and time is used by default.
        user : str, optional
            Username to include in the Sparky parameters. 'user' is the
            default.

        Returns
        -------
        dic : dict
            Dictionary of Sparky parameters.
        data : array_like
            NMR data in Sparky format.

        """
        # create dictionary
        dic = sparky.create_dic(self._udic, datetimeobj, user)

        # add processing flags for output
        self._oproc = {}
        self._odtype = "float32"

        return dic, self.__returndata()

    def to_bruker(self):
        """
        Return Bruker format data.

        Returns
        -------
        dic : dict
            Dictionary of Bruker parameters.
        data : array_like
            NMR data in Bruker format.

        """
        # create dictionary
        dic = bruker.create_dic(self._udic)

        # add processing flags for output
        self._oproc = {}
        self._odtype = "complex128"

        return dic, self.__returndata()


class udata_nd(fileiobase.data_nd):
    """
    Wrap other fileiobase.data_nd derived objects with input/output conversion
    when slices are requested.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    edata : fileiobase.data_nd derived object
        Data object to wrap.
    iproc : dict
        Dictionary of processing required by input format.
    oproc :
        Dictionary of processing required by output format.
    odtype : dtype
        Output dtype.
    order : tuple
        Axis ordering relative to input data.

    Notes
    -----
    The iproc and oproc dictionary can contains the following keys and values.

    ===========     ==========  ==========================================
    key             value       Description
    ===========     ==========  ==========================================
    alt_id_sign     True/False  True alternates signs along indirect dims.
    realfactor      float       Real channel scaling factor.
    imagfactor      float       Imaginary channel scaling factor.
    ===========     ==========  ==========================================

    """

    def __init__(self, edata, iproc, oproc, odtype, order=None):
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

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = udata_nd(self.edata, self._iproc, self._oproc, self._odtype, order)
        return n

    def __fgetitem__(self, slices):
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
                nd_iter = fileiobase.ndtofrom_iter(ffshape, fslice)
                for out_index, in_index in nd_iter:
                    # negate the trace if there is an odd number of
                    # odd number indices in the slice
                    if np.mod(in_index, 2).sum() % 2 == 1:
                        data[out_index] = -data[out_index]

        if "realfactor" in self._iproc:
            data.real = data.real * self._iproc['realfactor']
        if "imagfactor" in self._iproc:
            data.imag = data.imag * self._iproc['imagfactor']

        # output processing
        if "alt_id_sign" in self._oproc:
            if "alt_id_sign" not in self._iproc:
                fslice = slices[:-1]
                ffshape = self.fshape[:-1]
                nd_iter = fileiobase.ndtofrom_iter(ffshape, fslice)
                for out_index, in_index in nd_iter:
                    # negate the trace if there is an odd number of
                    # odd number indices in the slice
                    if np.mod(in_index, 2).sum() % 2 == 1:
                        data[out_index] = -data[out_index]

        if "realfactor" in self._oproc:
            data.real = data.real * self._oproc['realfactor']
        if "imagfactor" in self._oproc:
            data.imag = data.imag * self._oproc['imagfactor']

        return data.astype(self._odtype)
