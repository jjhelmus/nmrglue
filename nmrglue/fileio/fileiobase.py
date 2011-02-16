"""
fileiobase provides general purpose NMR file IO functions and classes
used by multiple nmrglue.fileio modules

"""

import numpy as np
import os
import string



def create_blank_udic(ndim):
    """ 
    Create a blank universal dictionary for a ndim-D spectrum
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

        d["label"] = ["X","Y","Z","A"][i]   # name of dimension

        # encoding of dimension, ie states, tppi, etc.  The direct dimension
        # should be listed as direct.
        if i==ndim-1:
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
    object provides methods to convert between common nmr units
    """

    def __init__(self,size,cplx,sw,obs,car):
        """ 
        create and set up a unit_conversion object

        Parameters:
        * size    Data size in points (R|I)
        * cplx    Complex flag (True or False)
        * sw      Spectral Width in Hz
        * obs     Observation Frequency in MHz
        * car     Carrier frequency in Hz

        """

        # fundamental units
        self._size = size
        self._cplx = cplx
        self._sw   = sw
        self._obs  = obs
        self._car  = car

        # derived units (these are in ppm)
        self._delta = -self._sw/(self._size*self._obs)
        self._first = self._car/self._obs - self._delta*self._size/2.


    # individual unit conversion functions
    def __percent2pts(self,percent):
        return percent*(self._size-1)/100.0

    def __pts2percent(self,pts):
        return pts*100/(self._size-1.0)

    def __hz2pts(self,hz):
        return ((hz/self._obs)-self._first)/self._delta

    def __pts2hz(self,pts):
        return (pts*self._delta+self._first)*self._obs

    def __ppm2pts(self,ppm):
        return (ppm-self._first)/self._delta   
    
    def __pts2ppm(self,pts):
        return (pts*self._delta)+self._first

    # times based units: seconds, ms, and us
    def __sec2pts(self,sec):
        return sec*self._sw

    def __pts2sec(self,pts):
        return pts*1./self._sw

    def __ms2pts(self,ms):
        return ms*self._sw/1.e3

    def __pts2ms(self,pts):
        return pts*1.e3/self._sw

    def __us2pts(self,us):
        return us*self._sw/1.e6

    def __pts2us(self,pts):
        return pts*1.e6/self._sw


    # routers
    def __unit2pnt(self,val,units):
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

        if self._cplx:
            return pts+round(pts)
        else:
            return pts


    def __pnt2unit(self,val,units):
        """ 
        Convert points to units
        """
        units = units.upper()

        if self._cplx:
            val = val-round(val)

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


    def __str2pnt(self,s):
        """ 
        Convert string with units to points
        """

        units = s.strip(string.digits+string.whitespace+"."+"-").upper()
        val   = float(s.strip(string.ascii_letters+string.whitespace+"%"))

        return self.__unit2pnt(val,units)
    

    def __convert(self,val,unit=None):
        """
        Convert string or value/unit pair
        """

        if type(val) == str:
            return self.__str2pnt(val)
        else:
            if unit==None:
                raise ValueError("invalid unit type")
            return self.__unit2pnt(val,unit)

    # User functions
    def f(self,val,unit=None):
        """
        Convert string or value/unit pair to float
        """
        return self.__convert(val,unit)

    def i(self,val,unit=None):
        """
        convert string or value/unit pair to integer
        """
        return int(round(self.__convert(val,unit)))

    def ppm(self,val):
        """
        Convert to ppms
        """
        return self.__pnt2unit(val,"PPM")

    def hz(self,val):
        """
        Convert to Hz
        """
        return self.__pnt2unit(val,"HZ")

    def percent(self,val):
        """
        Convert to percent
        """
        return self.__pnt2unit(val,"PERCENT")

    def seconds(self,val):
        """ 
        Convert to seconds
        """
        return self.__pnt2unit(val,"SEC")

    def sec(self,val):
        """ 
        Convert to seconds
        """
        return self.__pnt2unit(val,"SEC")

    def ms(self,val):
        """ 
        Convert to milliseconds (ms)
        """
        return self.__pnt2unit(val,"MS")

    def us(self,val):
        """ 
        Convert to microseconds (us)
        """
        return self.__pnt2unit(val,"US")

    def unit(self,val,unit):
        """
        Convert val points to unit
        """
        return self.__pnt2unit(val,unit)


    # limits and scales
    def percent_limits(self):
        """ 
        Return tuple of left and right edges in percent
        """
        return 0.0,100.0
    
    def percent_scale(self):
        """
        Return array of percent values
        """
        return linspace(0.0,100.0,self._size)

    def ppm_limits(self):
        """
        Return tuple of left and right edges in ppm
        """
        return self.ppm(0),self.ppm(self._size-1)

    def ppm_scale(self):
        """
        Return array of ppm values
        """
        x0,x1 = self.ppm_limits()
        return np.linspace(x0,x1,self._size)

    def hz_limits(self):
        """ 
        Return tuple of left and right edges in Hz
        """
        return self.hz(0),self.hz(self._size-1)

    def hz_scale(self):
        """
        Return array of hz values
        """
        x0,x1 = self.hz_limits()
        return np.linspace(x0,x1,self._size)

    def sec_limits(self):
        """ 
        Return tuple of left and right edges in seconds
        """
        return self.sec(0),self.sec(self._size-1)

    def sec_scale(self):
        """
        Return array of seconds values
        """
        x0,x1 = self.sec_limits()
        return np.linspace(x0,x1,self._size)

    def ms_limits(self):
        """ 
        Return tuple of left and right edges in milliseconds
        """
        return self.ms(0),self.ms(self._size-1)

    def ms_scale(self):
        """
        Return array of seconds values
        """
        x0,x1 = self.ms_limits()
        return np.linspace(x0,x1,self._size)

    def us_limits(self):
        """ 
        Return tuple of left and right edges in milliseconds
        """
        return self.us(0),self.us(self._size-1)

    def us_scale(self):
        """
        Return array of seconds values
        """
        x0,x1 = self.us_limits()
        return np.linspace(x0,x1,self._size)

    __call__ = i    # calling the object x is the same as x.i


def open_towrite(filename,overwrite=False):
    """ 
    Open filename for writing and return file object

    Function checks if file exists (and raises IOError if overwrite=False) and
    creates necessary directiories as needed.

    """

    # check if file exists and overwrite if False
    if os.path.exists(filename) and (overwrite==False):
        raise IOError,"File exists, recall with overwrite=True"
        return

    p,fn = os.path.split(filename)  # split into filename and path
    # create directories if needed
    if p != '' and os.path.exists(p) == False:
        os.makedirs(p)

    return open(filename,'w')


#
#data_* class primatives
#
#inherited classes should define:
#
#    __init__ which sets up the object and defines:
#        
#        self.lenX
#        self.lenY
#        self.lenZ if data_3d or data_4d
#        self.lenA if data_4d
#        self.order list of "x","y",{"z","a"}
#        self.shape
#        self.dtype
#        self.ndim
#
#    __fgetitem__ which takes well formatted tuples of slices
#    and returns ndarray objects
#
#    __fcopy__ which creates a copy provided only self and order parameters
#


class data_2d(object):
    """
    data_2d emulates numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects
    * can iterate over with expected results
    * transpose and swapaxes functions create a new data_2d object with the
      new axes ordering
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,order):
        """
        Create and set up a data_2d object
        """
        pass


    def __copy__(self):
        """
        Create a copy
        """

        return __fcopy__(self,self.order)


    def __getitem__(self,key):
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
            for j in range(2-len(rlist)):
                rlist.insert(i,slice(None))

        if len(rlist) > 2:
            raise IndexError,"invalid index"

        # replace integers with slices
        for i,v in enumerate(rlist):
            if type(v) == int:
                
                # check for out of range indexes
                if v >= self.shape[i]:
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v <= (-1*self.shape[i]-1):
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v < 0:
                    w  = self.shape[i]+v
                    rlist[i] = slice(w,w+1,1)
                else:
                    rlist[i] = slice(v,v+1,1)
            
        # pad the list with additional dimentions
        for i in range(len(rlist),2):
            rlist.append(slice(None))

        # reorder the slices into z,y,x
        sy = rlist[self.order.index("y")]
        sx = rlist[self.order.index("x")]

        # get the data
        data = self.__fgetitem__( (sy,sx) )

        # reorder the data
        if data.shape != (0,):
            a = [ ["y","x"].index(n) for n in self.order ]
            return np.squeeze(data.transpose(a))
        else:
            data

    def __len__(self):
        """
        x._len__ <==> len(x)
        """
        return self.shape[0]

    def __iter__(self):
        for index in xrange(0,self.shape[0]):
            yield self[index]

    def swapaxes(self,axis1,axis2):
        """ 
        Return a data_2d object with axis1 and axis2 interchanged

        Parameters:
        * axis1 First axis
        * axis2 Second axis

        """

        axis1,axis2 = int(axis1),int(axis2)

        if axis1 < 0:
            axis1 = 2-axis1
        if axis2 < 0:
            axis2 = 2-axis2
        if axis1 >= 2:
            raise ValueError,"bad axis1 argument to swapaxes"
        if axis2 >= 2:
            raise ValueError,"bad axis2 argument to swapaxes"

        order = list(self.order)
        order[axis1],order[axis2] = order[axis2],order[axis1]
        n = self.__fcopy__(order=order)

        return n

    def transpose(self,(axis1,axis2)=(1,0)):
        """
        Transpose data
        """
        ax1,ax2 = int(axis1),int(axis2)

        if ax1 < 0:
            ax1 = 2-ax1
        if ax2 < 0:
            ax2 = 2-ax2

        if ax1 == ax2:
            raise ValueError, "repeated axis in transpose"

        if ax1>=2 or ax2>=2:
            raise ValueError, "invalid axis for this array"

        order = list(self.order)
        new_order = [ order[ax1],order[ax2] ]
        n = self.__fcopy__(order=new_order)
        
        return n



class data_3d(object):
    """
    data_3d emulates numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects
    * can iterate over with expected results
    * transpose and swapaxes functions create a new fid_3d object with the
      new axes ordering
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,order):

        pass

    def __copy__(self):
        """ 
        create a copy
        """
        return __fcopy(self,self.order)

    def __getitem__(self,key):
        """
        x.__getitem__(y) <==> x[y]
        """

        # formats the input into a formated 

        # convert the key into a list
        if type(key) != tuple:
            rlist = [key]
        else:
            rlist = list(key)

        # remove Ellipsis
        while Ellipsis in rlist:
            i = rlist.index(Ellipsis)
            rlist.pop(i)
            for j in range(3-len(rlist)):
                rlist.insert(i,slice(None))

        if len(rlist) > 3:
            raise IndexError,"invalid index"

        # replace integers with slices
        for i,v in enumerate(rlist):
            if type(v) == int:
                
                # check for out of range indexes
                if v >= self.shape[i]:
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v <= (-1*self.shape[i]-1):
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v < 0:
                    w  = self.shape[i]+v
                    rlist[i] = slice(w,w+1,1)
                else:
                    rlist[i] = slice(v,v+1,1)
            
        # pad the list with additional dimentions
        for i in range(len(rlist),3):
            rlist.append(slice(None))

        # reorder the slices into z,y,x
        sz = rlist[self.order.index("z")]
        sy = rlist[self.order.index("y")]
        sx = rlist[self.order.index("x")]

        # get the data
        data = self.__fgetitem__( (sz,sy,sx) )

        # reorder the data
        if data.shape != (0,):
            a = [ ["z","y","x"].index(n) for n in self.order ]
            return np.squeeze(data.transpose(a))
        else:
            data

    def __len__(self):
        """
        x._len__ <==> len(x)
        """
        return self.shape[0]

    def __iter__(self):
        for index in xrange(0,self.shape[0]):
            yield self[index]

    def swapaxes(self,axis1,axis2):
        """
        Return fid_3d object with axis1 and axis2 interchanged

        Parameters:
        * axis1 First axis
        * axis2 Second axis

        """

        axis1,axis2 = int(axis1),int(axis2)

        if axis1 < 0:
            axis1 = 3-axis1
        if axis2 < 0:
            axis2 = 3-axis2
        if axis1 >= 3:
            raise ValueError,"bad axis1 argument to swapaxes"
        if axis2 >= 3:
            raise ValueError,"bad axis2 argument to swapaxes"

        order = list(self.order)
        order[axis1],order[axis2] = order[axis2],order[axis1]
        n = self.__fcopy__(order=order)

        return n

    def transpose(self,(axis1,axis2,axis3)=(2,1,0)):

        ax1,ax2,ax3 = int(axis1),int(axis2),int(axis3)

        if ax1 < 0:
            ax1 = 3-ax1
        if ax2 < 0:
            ax2 = 3-ax2
        if ax3 < 0:
            ax3 = 3-ax2

        if ax1 == ax2 or ax1 == ax3 or ax2 == ax3:
            raise ValueError, "repeated axis in transpose"

        if ax1>=3 or ax2>=3 or ax3>=3:
            raise ValueError, "invalid axis for this array"

        order = list(self.order)
        new_order = [ order[ax1],order[ax2],order[ax3] ]
        n = self.__fcopy__(order=new_order)
        
        return n


class data_4d(object):
    """
    data_4d emulates numpy.ndarray object without loading data into memory

    * slicing operations return ndarray objects
    * can iterate over with expected results
    * transpose and swapaxes functions create a new data_4d object with the
      new axes ordering
    * has ndim, shape, and dtype attributes.

    """

    def __init__(self,order):

        pass

    def __copy__(self):
        """ 
        create a copy
        """
        return __fcopy(self,self.order)

    def __getitem__(self,key):
        """
        x.__getitem__(y) <==> x[y]
        """

        # formats the input into a formated 

        # convert the key into a list
        if type(key) != tuple:
            rlist = [key]
        else:
            rlist = list(key)

        # remove Ellipsis
        while Ellipsis in rlist:
            i = rlist.index(Ellipsis)
            rlist.pop(i)
            for j in range(4-len(rlist)):
                rlist.insert(i,slice(None))

        if len(rlist) > 4:
            raise IndexError,"invalid index"

        # replace integers with slices
        for i,v in enumerate(rlist):
            if type(v) == int:
                
                # check for out of range indexes
                if v >= self.shape[i]:
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v <= (-1*self.shape[i]-1):
                    raise IndexError,"index(%s) out of range(0<=index<%s) \
                    in dimension %s" % (v,self.shape[i]-1,i)

                if v < 0:
                    w  = self.shape[i]+v
                    rlist[i] = slice(w,w+1,1)
                else:
                    rlist[i] = slice(v,v+1,1)
            
        # pad the list with additional dimentions
        for i in range(len(rlist),4):
            rlist.append(slice(None))

        # reorder the slices into z,y,x
        sa = rlist[self.order.index("a")]
        sz = rlist[self.order.index("z")]
        sy = rlist[self.order.index("y")]
        sx = rlist[self.order.index("x")]

        # get the data
        data = self.__fgetitem__( (sa,sz,sy,sx) )

        # reorder the data
        if data.shape != (0,):
            a = [ ["a","z","y","x"].index(n) for n in self.order ]
            return np.squeeze(data.transpose(a))
        else:
            data

    def __len__(self):
        """
        x._len__ <==> len(x)
        """
        return self.shape[0]

    def __iter__(self):
        for index in xrange(0,self.shape[0]):
            yield self[index]

    def swapaxes(self,axis1,axis2):
        """
        Return fid_3d object with axis1 and axis2 interchanged

        Parameters:
        * axis1 First axis
        * axis2 Second axis

        """

        axis1,axis2 = int(axis1),int(axis2)

        if axis1 < 0:
            axis1 = 4-axis1
        if axis2 < 0:
            axis2 = 4-axis2
        if axis1 >= 4:
            raise ValueError,"bad axis1 argument to swapaxes"
        if axis2 >= 4:
            raise ValueError,"bad axis2 argument to swapaxes"

        order = list(self.order)
        order[axis1],order[axis2] = order[axis2],order[axis1]
        n = self.__fcopy__(order=order)

        return n

    def transpose(self,(axis1,axis2,axis3,axis4)=(3,2,1,0)):

        ax1,ax2,ax3,ax4 = int(axis1),int(axis2),int(axis3),int(axis4)

        if ax1 < 0:
            ax1 = 4-ax1
        if ax2 < 0:
            ax2 = 4-ax2
        if ax3 < 0:
            ax3 = 4-ax2
        if ax4 < 0:
            ax4 = 4-ax4



        if (ax1==ax2 or ax1==ax3 or ax1==ax4 or 
            ax2==ax3 or ax2==ax4 or ax3==ax4) :
            raise ValueError, "repeated axis in transpose"

        if ax1>=4 or ax2>=4 or ax3>=4 or ax4>=4:
            raise ValueError, "invalid axis for this array"

        order = list(self.order)
        new_order = [ order[ax1],order[ax2],order[ax3],order[ax4] ]
        n = self.__fcopy__(order=new_order)
        
        return n
