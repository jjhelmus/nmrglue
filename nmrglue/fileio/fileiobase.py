"""
fileiobase provides general purpose NMR file IO functions and classes
used by multiple nmrglue.fileio modules.
"""

import os
import string
import itertools

import numpy as np


def create_blank_udic(ndim):
    """
    Create a blank universal dictionary for a spectrum of dimension ndim.
    """
    udic = dict()
    udic["ndim"] = ndim

    for i in xrange(ndim):
        d = dict()
        d["sw"] = 999.99        # spectral width in Hz
        d["complex"] = True     # Quadrature, True when dimension is complex
        d["obs"] = 999.99       # Observation frequency in MHz
        d["car"] = 999.99       # Carrier frequency in Hz
        d["size"] = 1           # Number of points in dimension based on the
                                # shape of the data array. As such the direct
                                # dimension (-1) size is R|I, all indirect
                                # dimensions are R+I

        d["label"] = ["X", "Y", "Z", "A"][i]    # name of dimension

        # encoding of dimension, ie states, tppi, etc.  The direct dimension
        # should be listed as direct.
        if i == ndim - 1:
            d["encoding"] = "direct"
        else:
            d["encoding"] = "states"

        # time and freq flags for domain of dimension
        d["time"] = True
        d["freq"] = False

        udic[i] = d

    return udic


class unit_conversion():
    """
    Provides methods to convert between common NMR units

    Parameters
    ----------
    size : int
        Number of points in dimension (R|I).
    cplex : bool
        True if dimension is complex, False is real.
    sw : float
        Spectral width in Hz.
    obs : float
        Observation frequency in MHz.
    car : float
        Carrier frequency in Hz.

    """
    def __init__(self, size, cplx, sw, obs, car):
        """
        create and set up a unit_conversion object
        """
        # fundamental units
        self._size = size
        self._cplx = cplx
        self._sw = sw
        self._obs = obs
        self._car = car

        # derived units (these are in ppm)
        self._delta = -self._sw / (self._size * self._obs)
        self._first = self._car / self._obs - self._delta * self._size / 2.

    # individual unit conversion functions
    def __percent2pts(self, percent):
        return percent * (self._size - 1) / 100.0

    def __pts2percent(self, pts):
        return pts * 100 / (self._size - 1.0)

    def __hz2pts(self, hz):
        return ((hz / self._obs) - self._first) / self._delta

    def __pts2hz(self, pts):
        return (pts * self._delta + self._first) * self._obs

    def __ppm2pts(self, ppm):
        return (ppm - self._first) / self._delta

    def __pts2ppm(self, pts):
        return (pts * self._delta) + self._first

    # times based units: seconds, ms, and us
    def __sec2pts(self, sec):
        return sec * self._sw

    def __pts2sec(self, pts):
        return pts * 1. / self._sw

    def __ms2pts(self, ms):
        return ms * self._sw / 1.e3

    def __pts2ms(self, pts):
        return pts * 1.e3 / self._sw

    def __us2pts(self, us):
        return us * self._sw / 1.e6

    def __pts2us(self, pts):
        return pts * 1.e6 / self._sw

    # routers
    def __unit2pnt(self, val, units):
        """
        Convert units to points
        """
        units = units.upper()
        if units == "PPM":
            pts = self.__ppm2pts(val)
        elif units == "HZ":
            pts = self.__hz2pts(val)
        elif units == "%" or units == "PERCENT":
            pts = self.__percent2pts(val)
        elif units == "SEC" or units == "SECOND" or units == "S":
            pts = self.__sec2pts(val)
        elif units == "MS":
            pts = self.__ms2pts(val)
        elif units == "US":
            pts = self.__us2pts(val)
        else:
            raise ValueError("invalid unit type")
        #if self._cplx:
        #    return pts+round(pts)
        #else:
        return pts

    def __pnt2unit(self, val, units):
        """
        Convert points to units
        """
        units = units.upper()
        #if self._cplx:
        #    val = val-round(val)
        if units == "PPM":
            k = self.__pts2ppm(val)
        elif units == "HZ":
            k = self.__pts2hz(val)
        elif units == "%" or units == "PERCENT":
            k = self.__pts2percent(val)
        elif units == "SEC" or units == "SECOND" or units == "S":
            k = self.__pts2sec(val)
        elif units == "MS":
            k = self.__pts2ms(val)
        elif units == "US":
            k = self.__pts2us(val)
        else:
            raise ValueError("invalid units")
        return k

    def __str2pnt(self, s):
        """
        Convert string with units to points
        """
        units = s.strip(string.digits + string.whitespace + "." + "-").upper()
        val = float(s.strip(string.ascii_letters + string.whitespace + "%"))
        return self.__unit2pnt(val, units)

    def __convert(self, val, unit=None):
        """
        Convert string or value/unit pair
        """
        if type(val) == str:
            return self.__str2pnt(val)
        else:
            if unit == None:
                raise ValueError("invalid unit type")
            return self.__unit2pnt(val, unit)

    # User functions
    def f(self, val, unit=None):
        """
        Convert string or value/unit pair to float
        """
        return self.__convert(val, unit)

    def i(self, val, unit=None):
        """
        Convert string or value/unit pair to integer
        """
        return int(round(self.__convert(val, unit)))

    def ppm(self, val):
        """
        Convert to ppm
        """
        return self.__pnt2unit(val, "PPM")

    def hz(self, val):
        """
        Convert to Hz
        """
        return self.__pnt2unit(val, "HZ")

    def percent(self, val):
        """
        Convert to percent
        """
        return self.__pnt2unit(val, "PERCENT")

    def seconds(self, val):
        """
        Convert to seconds
        """
        return self.__pnt2unit(val, "SEC")

    def sec(self, val):
        """
        Convert to seconds
        """
        return self.__pnt2unit(val, "SEC")

    def ms(self, val):
        """
        Convert to milliseconds (ms)
        """
        return self.__pnt2unit(val, "MS")

    def us(self, val):
        """
        Convert to microseconds (us)
        """
        return self.__pnt2unit(val, "US")

    def unit(self, val, unit):
        """
        Convert val points to unit
        """
        return self.__pnt2unit(val, unit)

    # limits and scales
    def percent_limits(self):
        """
        Return tuple of left and right edges in percent
        """
        return 0.0, 100.0

    def percent_scale(self):
        """
        Return array of percent values
        """
        return np.linspace(0.0, 100.0, self._size)

    def ppm_limits(self):
        """
        Return tuple of left and right edges in ppm
        """
        return self.ppm(0), self.ppm(self._size - 1)

    def ppm_scale(self):
        """
        Return array of ppm values
        """
        x0, x1 = self.ppm_limits()
        return np.linspace(x0, x1, self._size)

    def hz_limits(self):
        """
        Return tuple of left and right edges in Hz
        """
        return self.hz(0), self.hz(self._size - 1)

    def hz_scale(self):
        """
        Return array of Hz values
        """
        x0, x1 = self.hz_limits()
        return np.linspace(x0, x1, self._size)

    def sec_limits(self):
        """
        Return tuple of left and right edges in seconds
        """
        return self.sec(0), self.sec(self._size - 1)

    def sec_scale(self):
        """
        Return array of seconds values
        """
        x0, x1 = self.sec_limits()
        return np.linspace(x0, x1, self._size)

    def ms_limits(self):
        """
        Return tuple of left and right edges in milliseconds
        """
        return self.ms(0), self.ms(self._size - 1)

    def ms_scale(self):
        """
        Return array of seconds values
        """
        x0, x1 = self.ms_limits()
        return np.linspace(x0, x1, self._size)

    def us_limits(self):
        """
        Return tuple of left and right edges in milliseconds
        """
        return self.us(0), self.us(self._size - 1)

    def us_scale(self):
        """
        Return array of seconds values
        """
        x0, x1 = self.us_limits()
        return np.linspace(x0, x1, self._size)

    __call__ = i    # calling the object x is the same as x.i


def open_towrite(filename, overwrite=False):
    """
    Open filename for writing and return file object

    Function checks if file exists (and raises IOError if overwrite=False) and
    creates necessary directiories as needed.
    """
    # check if file exists and overwrite if False
    if os.path.exists(filename) and (overwrite == False):
        raise IOError("File exists, recall with overwrite=True")

    p, fn = os.path.split(filename)  # split into filename and path
    # create directories if needed
    if p != '' and os.path.exists(p) == False:
        os.makedirs(p)

    return open(filename, 'wb')

################################################
# numpy ndarray emulation and helper functions #
################################################

# iterators for ND array


def ndfrom_iter(shape, slices):
    ch = [range(lenx)[sX] for lenx, sX in zip(shape, slices)]
    return itertools.product(*ch)


def ndto_iter(shape, slices):
    ich = [range(len(range(lenx)[sX])) for lenx, sX in zip(shape, slices)]
    return itertools.product(*ich)


def ndtofrom_iter(shape, slices):
    ch = [range(lenx)[sX] for lenx, sX in zip(shape, slices)]
    ich = [range(len(i)) for i in ch]
    return zip(itertools.product(*ich), itertools.product(*ch))


def size_and_ndtofrom_iter(shape, slices):
    ch = [range(lenx)[sX] for lenx, sX in zip(shape, slices)]
    s = [len(i) for i in ch]
    ich = [range(i) for i in s]
    return s, zip(itertools.product(*ich), itertools.product(*ch))


# index2trace and trace2index functions


def index2trace_flat(shape, index):
    """
    Calculate trace number from shape and index of all indirect dimensions
    assuming a flat structure
    """
    # We need to perform:
    # index[0]*shape[1]*...shape[-1] + index[1]*shape[2]*...shape[-1] + ...
    # + index[-1]*shape[-1] + index[-1]
    # To do this we calculate the product of shape[X] elements and multiple
    # by the corresponding index element, index[-1] as added at the beginning
    a = index[-1]
    for i, v in enumerate(index[:-1]):
        mult = reduce(lambda x, y: x * y, shape[i + 1:])
        a = a + mult * v
    return a


def trace2index_flat(shape, ntrace):
    """
    Calculate the index of a trace assuming a flat structure
    """
    # algorithm is to take quotient/remainers of sizes in reverse
    q = ntrace  # seed quotient with remained
    index = []
    for s in shape[:0:-1]:  # loop from last size to 2nd size
        q, r = divmod(q, s)
        index.insert(0, r)
    index.insert(0, q)
    return tuple(index)


def index2trace_opp(shape, index):
    """
    Calculate trace number from shape and index of all indirect dimensions
    assuming a phase ordering opposite the time increments.
    """
    n = len(shape)
    # deal with the phase component
    phases = [v % 2 for v in index]
    nphase = index2trace_flat([2] * n, phases[::-1])
    # deal with the remainer
    pindex = [v // 2 for v in index]
    pshape = [i / 2 for i in shape]
    nbase = index2trace_flat(pshape, pindex)
    return nbase * 2 ** n + nphase


def trace2index_opp(shape, ntrace):
    """
    Calculate the index of a trace assuming opposite phase/time increment
    ordering
    """
    n = len(shape)
    q, r = divmod(ntrace, 2 ** n)
    to_add = list(trace2index_flat([2] * n, r))[::-1]
    pshape = [i / 2 for i in shape]
    base = list(trace2index_flat(pshape, q))
    total = [b * 2 + a for b, a in zip(base, to_add)]
    return tuple(total)


def index2trace_reg(shape, index):
    """
    Calculate trace number from shape and index of all indirect dimensions
    assuming the same  phase and time ordering.
    """
    n = len(shape)
    # deal with the phase component
    phases = [v % 2 for v in index]
    nphase = index2trace_flat([2] * n, phases)
    # deal with the remainer
    pindex = [v // 2 for v in index]
    pshape = [i / 2 for i in shape]
    nbase = index2trace_flat(pshape, pindex)
    return nbase * 2 ** n + nphase


def trace2index_reg(shape, ntrace):
    """
    Calculate the index of a trace assuming the same phase/time increment
    ordering
    """
    n = len(shape)
    q, r = divmod(ntrace, 2 ** n)
    to_add = list(trace2index_flat([2] * n, r))
    pshape = [i / 2 for i in shape]
    base = list(trace2index_flat(pshape, q))
    total = [b * 2 + a for b, a in zip(base, to_add)]
    return tuple(total)

#
# data_nd class
#
# inherited classes should define:
#
#    __init__ which sets up the object and defines at minimum
#
#       self.fshape shape of data on disk (shape when order = (0,1,2...)
#       self.order order of axes, default is (0,1,2,...)
#       self.dtype
#
#   self.__setdimandshape__ can be called to set self.dim and self.shape
#    if they are not set by __init__
#
#    __fgetitem__ which takes well formatted tuples of slices
#    and returns ndarray objects
#
#    __fcopy__ which creates a copy provided only self and order parameters
#


class data_nd(object):
    """
    Base class for building objects which emulate ndarray objects without
    loading data into memory.  These object have the following properties:

    * slicing operations return ndarray objects
    * can iterate over with expected results
    * transpose and swapaxes functions create a new data_nd object with the
      new axes ordering
    * has ndim, shape, and dtype attributes.

    Notes
    -----

    Classes which are use this class as a base should define the following
    methods:

    __init__ which must set up the object and defines at minimum:

        self.fshape : tuple
            Shape of the data on disk, the shape when order = (0, 1, 2, ..)
        self.order : tuple
            Ordering of the axes
        self.dtype : dtype
            Dtype of the emulated ndarray

        self.__setdimandshape__ should be called if self.dim and self.shape
            are not set up __init__.

    __fgetitem__ which takes a well formatted tuple of slices and returns a
        ndarray object with the selected data

    __fcopy__(self, order) which created a copy of the object with the new
        order

    """

    def __init__(self, order):
        pass

    def __setdimandshape__(self):
        """ Set the the ndim and shape attributes from fshape """
        # set ndim and shape
        self.ndim = len(self.fshape)
        self.shape = tuple([self.fshape[i] for i in self.order])

    def __copy__(self):
        """
        create a copy
        """
        return __fcopy(self, self.order)

    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        # convert the key into a list
        if type(key) != tuple:
            rlist = [key]
        else:
            rlist = list(key)

        # remove Ellipsis
        while Ellipsis in rlist:
            i = rlist.index(Ellipsis)
            rlist.pop(i)
            for j in range(self.ndim - len(rlist)):
                rlist.insert(i, slice(None))

        if len(rlist) > self.ndim:
            raise IndexError("invalid index")

        # replace none-slices with slice objects
        for i, v in enumerate(rlist):
            if type(v) != slice:

                # check for out of range indexes
                if v >= self.shape[i]:
                    raise IndexError("index(%s) out of range(0 <= index < %s) \
                    in dimension %s" % (v, self.shape[i] - 1, i))

                if v <= (-1 * self.shape[i] - 1):
                    raise IndexError("index(%s) out of range(0 <= index < %s) \
                    in dimension %s" % (v, self.shape[i] - 1, i))

                if v < 0:
                    w = self.shape[i] + v
                    rlist[i] = slice(w, w + 1, 1)
                else:
                    rlist[i] = slice(v, v + 1, 1)

        # pad the list with additional dimentions
        for i in range(len(rlist), self.ndim):
            rlist.append(slice(None))

        # reorder the slices into file order
        frlist = [rlist[self.order.index(i)] for i in range(self.ndim)]

        # get the data from disk
        data = self.__fgetitem__(tuple(frlist))

        # re-order the data
        if data.shape != (0,):
            return np.squeeze(data.transpose(self.order))
        else:
            return data

    def __len__(self):
        """
        x._len__ <==> len(x)
        """
        return self.shape[0]

    def __iter__(self):
        """ x.__iter__() <==> iter(x) """
        for index in xrange(0, self.shape[0]):
            yield self[index]

    def swapaxes(self, axis1, axis2):
        """
        Return object with `axis1` and `axis2` interchanged.
        """
        axis1, axis2 = int(axis1), int(axis2)
        if axis1 < 0:
            axis1 = self.ndim - axis1
        if axis2 < 0:
            axis2 = self.ndim - axis2
        if axis1 >= self.ndim:
            raise ValueError("bad axis1 argument to swapaxes")
        if axis2 >= self.ndim:
            raise ValueError("bad axis2 argument to swapaxes")

        order = list(self.order)
        order[axis1], order[axis2] = order[axis2], order[axis1]
        n = self.__fcopy__(order=order)
        return n

    def transpose(self, *axes):
        """
        Return object with axes transposed.

        Parameters
        ----------
        axes : None, tuple or ints, or `n` ints
            * None or no arguments: reverse order of the axes

            * tuple of ints: `i` in the `j`-th place in the tuples means the
            'i'-th axis becomes the new objects `j`-th axis.

            * `n` ints: same as an n-tuple.

        Returns
        -------
        out : data_nd object
            Object whose axes are permuted.

        """
        if axes == ():    # default is to switch order of axes
            axes = range(self.ndim)[::-1]

        if len(axes) == 1:    # if a single tuple is given unpack
            axes = axes[0]

        try:    # convert to integers
            axes = [int(i) for i in axes]
        except:
            raise TypeError("an integer is required")

        if len(axes) != self.ndim:   # check for to few/many axes
            raise ValueError("axes don't match array")

        # replace negatives axes values with positives
        for i, v in enumerate(axes):
            if v < 0:
                axes[i] = self.ndim + v

        # check for invalid axes
        for v in axes:
            if v >= self.ndim:
                raise ValueError("invalid axis for this array")

        # check for repeated axes
        if len(set(axes)) != self.ndim:
            raise ValueError("repeated axis in tranpose")

        # create a new data_nd object with transposed order
        return self.__fcopy__(order=tuple([self.order[i] for i in axes]))
