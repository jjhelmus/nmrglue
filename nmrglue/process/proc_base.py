"""
Base NMR spectral processing functions

proc_base
=========

Provides basic NMR spectral processing functions.  These functions act on the -1 (last)
axis.

Documentation is available in the docstrings and at http://XXX

To do:  refactor using code from pipe_proc
        documentation and cleanup
        replace numpy.fft with scipy.fftpack requires rev 5355
        rewrite ha function using FWHT

History
(jjh) 2009.10.02 refactoring
(jjh) 2009.05.05 forked from process.py

"""

# standard library modules

# external modules
import numpy as np

# nmrglue modules


pi = np.pi


# Hadamard Transform functions

def int2bin(n,digits=8):
    """ integer to binary string"""
    return "".join([str((n >> y) & 1) for y in range(digits-1, -1, -1)])

def bin2int(s):
    """ binary string to integer"""
    o = 0
    k = len(s)-1
    for i,v in enumerate(s):
        o = o + int(v)*2**(k-i)
    return o

def gray(n):
    """
    Calculate n-bit gray code
    """

    g = [0,1]
    
    for i in range(1,int(n)):
        mg = g+g[::-1]   # mirror the current code
        first = [0]*2**(i) + [2**(i)]*2**(i)    # first bit 0/2**i for mirror
        g = [mg[j] + first[j] for j in range(2**(i+1))]

    return g

def hadamard(n, dtype=int):
    """hadamard(n) returns Hadamard matrix of order n, where n must be power of 2

    Source: Scipy Dev Wiki by: Ivo
    http://projects.scipy.org/scipy/attachment/ticket/675/hadamard.py
    """
   
    ord     = np.log2(n)
    int_ord = int(ord)
  
    if ord != int_ord:
        raise ValueError("n must be an integer, and n must be power of 2")
  
    H = np.array([1],dtype=dtype)
    if n == 1:
        return H
  
    # Sylvester's construction
    for i in range(0,int_ord): 
        H = np.c_[np.r_[H, H], np.r_[H, -H]]
  
    return H

def ha(data):
    """ 
    Hadamard transform
    """
    # What we really want to do is Fast Walsh-Hadamard Transform (FWHT)
    # with sequency/Walsh ordering ie FWHT_w
    # See http://en.wikipedia.org/wiki/Walsh_matrix
    # http://en.wikipedia.org/wiki/Fast_Hadamard_transform


    # The current implementation is a proof of concept and EXTEMEMLY SLOW



    # determind the order and final size of input vectors
    ord = int(np.ceil(np.log2(data.shape[-1])))  # Walsh/Hadamard order 
    max = 2**ord

    # zero fill to power of 2
    pad = max - data.shape[-1]
    zdata = zf(data,pad)

    # Multiple each vector by the hadamard matrix
    nat = np.zeros(zdata.shape,dtype=zdata.dtype)
    H = hadamard(max)
    # XXX lots of optimization can be done here XXX
    nat = np.dot(zdata,H)
    nat = np.array(nat,dtype=data.dtype)

    # Bit-Reveral  Permutation 
    s = [int2bin(x,digits=ord)[::-1] for x in range(max)]
    brp = [bin2int(x) for x in s]
    brp_data = np.take(nat,brp,axis=-1)

    # Gray code permutation (bit-inverse)
    gp = gray(ord)
    gp_data = np.take(brp_data,gp,axis=-1)

    return gp_data

# Packing/Unpacking Functions

def unpack_complex(data):
    """
    unpacks complex array into real array

    Note:   interleaves values ie. (r,i,r,i....)
    """
    
    s = list(data.shape)
    s[-1] = s[-1]*2
    udata = np.zeros(s,dtype=data.real.dtype)
    udata[...,::2] = data.real
    udata[...,1::2] = data.imag

    return udata

def pack_complex(data):
    """
    packs real array into complex array 

    Note: Complex values must be interleaved
          dtype should be set outside this function.
    """

    return data[...,::2]+data[...,1::2]*1.j

def decode_States(data):
    """
    decode data collected using States

    Note: dtype should be set outside this function
    """
    # this is really just a hook to call pack_complex
    return pack_complex(data)

def decode_TPPI(data):
    """
    decode data collected using TPPI
    """

    raise NotImplemented

# Apodization Functions

def em(data,lb=0.0,inv=False):
    """
    exponential apodization on data

    em(x) = exp(-pi*i*lb)

    Parameters:
    data    ndarray containing data to process
    lb      Exponential line broadening
    inv     Invert apodization

    Returns ndarray object
    """

    apod = np.exp(-pi*np.arange(data.shape[-1])*lb,sig=data.dtype)
    if inv:
        apod = 1/apod   # invert the window

    return apod*data

def gm(data,g1=0.0,g2=0.0,g3=0.0,inv=False):
    """ 
    lorentz-to-gauss apodization
    
    gm(x_i) = exp(e - g*g)
    
    where:  e = pi*i*g1
            g = 0.6*pi*g2*(g3*(size-1)-i)

    Parameters:
    data    ndarray containing data to process
    g1      Inverse Exponential Width
    g2      Gaussian Broaden Width
    g3      Location of Gauss Maximum

    Returns ndarray object
    """
    size = data.shape[-1]
    e = pi*np.arange(size)*g1
    g = 0.6*pi*g2*(g3*(size-1) - np.arange(size))
    apod = np.exp(e-g*g, sig = data.dtype)
    if inv:
        apod = 1/apod

    return apod*data

def gmb(data,a=0.0,b=0.0,inv=False):
    """ 
    Gaussian apodization 

    gmb(x_i) = exp(-a*i - b*i*i)

    Parameters
    a   exponential term
    b   gaussian term
    inv invert apodization

    returns ndarray object
    """
        
    size = data.shape[-1]
    apod = np.exp(-a*np.arange(size)-b*np.arange(size)**2,sig=data.dtype)
    if inv:
        apod = 1/apod

    return apod*data


def jmod(data,e=0.0,off=0.0,end=0.0,inv=False):
    """ 
    Exponentially Damped J-Modulation

    jmod(x_i) = exp (-e)*sin( pi*off + pi*(end-off)*i/(size-1))

    Parameters
    e   exponential term
    off modulation start
    end modulation end
    inv invert apodization

    returns ndarray object
    """

    size = data.shape[-1]
    apod = np.exp(-e*np.arange(size),sig=data.dtype) * np.sin(pi*off+pi*(end-off)*np.arange(size)/(size-1),sig=data.dtype)
    if inv:
        apod = 1/apod

    return apod*data


def sp(data,off=0,end=1.0,pow=1.0,inv=False):
        """ 
        shifted sine bell apodization on data

        Parameters:
        data    ndarray containing data to process
        start   Start value of Sine apodization (0.0 -> 1.0)
        end     End value of Sine apodizations (0.0 -> 1.0 )
        pow     Power of Sine bell (1.0)

        Returns ndarray object
        """

        size = data.shape[-1]

        # SP apodization given by: 
        # sin( (pi*off + pi*(end-off)*i/(size-1) )^pow

        apod = np.power(np.sin(pi*off+pi*(end-off)*np.arange(size)/(size-1),sig=data.dtype) ,pow,sig=data.dtype)
        if inv:
            apod = 1/apod
        return apod*data

def tm(data,t1=0.0,t2=0.0,inv=False):
    """
    trapezoid apodization

    t1 and t2 define number of points to build 0.0 to 1.0 linear ramp

    Parameters:
    data    ndarray contraining data to process
    t1      left ramp length
    t2      right ramp length
    inv     invert apodization

    returns ndarray object
    """

    size = data.shape[-1] 
    apod = np.concatenate( (np.linspace(0,1,t1),np.ones(size-t1-t2),np.linspace(1,0,t2)) )
    apod = np.array(apod,dtype=data.dtype)  # type cast array
    if inv:
        apod = 1/apod

    return apod*data

def tri(data,loc="auto",lHi=0.0,rHi=0.0,inv=False):
    """
    triangle window

    triangle window from lHi to 1.0  to rHi centered at loc

    Parameters:
    data    ndarray containing data to process
    loc     location of apex
    lHi     Left hand size starting height
    rHi     Right hand size starting height

    returns ndarray object
    """

    size = data.shape[-1]
    if loc=="auto":
        loc = size/2
    apod = np.concatenate( (np.linspace(lHi,1.0,loc),np.linspace(1.0,rHi,size-loc+1)[1:] ) )
    apod = np.array(apod,dtype=data.dtype)
    if inv:
        apod = 1/apod

    return apod*data


# Shift/Roll Functions

def rs(data,pts=0.0):
    """
    right shift and fill with zero
    
    Parameters:
    data    ndarray object
    pts     number of points to right shift

    returns ndarray object

    Note: fills shifted points with zeros
    """

    data = np.roll(data,int(pts),axis=-1)
    data[...,:int(pts)] = 0
    return data

def ls(data,pts=0.0):
    """
    right shift and fill with zero
    
    Parameters:
    data    ndarray object
    pts     number of points to right shift

    returns ndarray object

    Note: fills shifted points with zeros
    """

    data = np.roll(data,-int(pts),axis=-1)
    data[...,-int(pts):] = 0
    return data

def roll(data,pts=0.0,neg=False):
    """

    Roll axis optionally negating shifted points
    """
    
    data = np.roll(data,int(pts),axis=-1)

    if neg:
        if pts > 0:
            data[...,:pts] = -data[...,:pts]
        else:
            data[...,pts:] = -data[...,pts:]

    return data



def zf_inter(data,zf=1):
    """
    zero fill zf zeros between points

    Parameters
    data    ndarray object
    pad     number zeros between points

    returns ndarray objects
    """

    size = list(data.shape)
    size[-1] = (zf+1)*size[-1]
    z = np.zeros(size,dtype=data.dtype)
    z[...,::zf+1] = data[...,:] 

    return z

def zf(data,pad=0,mid=False):
    """
    zero fill data with pad zeros

    Parameters
    data  ndarray object
    pad   number of zeros to add
    mid   set to True to zero fill in middle of data

    returns ndarray object
    """
    
    size = list(data.shape)
    size[-1] = pad
    z = np.zeros(size,dtype=data.dtype)

    if mid:
        h = int(data.shape[-1]/2.0)
        return np.concatenate( (data[...,:h],z,data[...,h:]),axis=-1)

    else:
        return np.concatenate( (data,z),axis=-1)


def nmr_reorder(v):
    """ reorder spectrum after transform  

    swap halves and reverses
    """

    s = v.shape[-1]
    return np.append(v[...,s/2::-1],v[...,s:s/2:-1],axis=-1)
    

def rft(x):
    """ real Fourier transform """

    s = x.shape[-1]
    xp = np.zeros(x.shape,dtype="complex64")
    xp[...,1:s/2] = x[...,1:-1:2]+x[...,2::2]*1.j
    xp[...,0] = x[...,0]/2.
    xp[...,s/2] = x[...,-1]/2.

    return np.array(nmr_reorder(np.fft.fft(2*xp,axis=-1).real),dtype="float32")


def irft(xp):
    """inverse real Fourier transform """
    s = xp.shape[-1]
    xp = np.fft.ifft(nmr_reorder(xp))   # re-order, inverse FT
    
    # output results
    x = np.zeros(xp.shape,dtype="float32")

    # unpack ifft data
    x[...,1:-1:2] = xp[...,1:s/2].real
    x[...,2::2]   = xp[...,1:s/2].imag
    x[...,0]      = xp[...,0].real
    x[...,-1]     = xp[...,s/2].real

    return x


def icomplexft(data):
    """ inverse nmr fft on complex data """
    # fft and recast to correct dtype
    size = data.shape[-1]
    data = np.array(np.fft.ifft(data),dtype=data.dtype) 
    # keep zero-freq term first, reverse rest of spectrum and sign alt.
    data = np.roll(data,-1,axis=-1)[...,::-1]
    data[...,1::2] = data[...,1::2]*-1.
       
    return data


def complexft(data):
    """ nmr fft on complex data """
    # fft and recast to correct dtype
    data = np.array(np.fft.fft(data),dtype=data.dtype)
    # swap halves and reverse
    # XXX might be able to speed this up with fancy indexing
    # check out ../numpy/fft/helper.py for examples
    size = data.shape[-1]
    data = np.append(data[...,size/2::-1],data[...,size:size/2:-1],axis=-1) 
       
    return data

def fsh(data,pts):
    """ Frequency Shift data 
    
    pts is number of points to right shift data, use negative value to
    left shift data
    """

    s= float(data.shape[-1])

    #inverse fft -> first order phase correction -> fft
    #idata = icomplexft(data)
    #pdata = np.exp(-2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*icomplexft(data)
    #data = complexft(pdata)

    # fast version 
    data = complexft(np.exp(-2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*
            icomplexft(data))


    return data
    


#def nmr_dct(data):
#    """ NMR ordered Discrete Cosine ('real ft') Transform"""
#
#    # XXX rewrite this when DCT wrapper to fftpack included in scipy 0.8
#    data = np.array(transforms.dct(data),dtype=data.dtype)
#    # swap halves and reverse
#    size = data.shape[-1]
#    data = np.append(data[...,size/2::-1],data[...,size:size/2:-1],axis=-1)
#
#    return data
#
#def nmr_idct(data):
#    """ NMR ordered Discrete Cosine ('real ft') Transform"""
#
#    # XXX rewrite this when DCT wrapper to fftpack included in scipy 0.8
#    size = data.shape[-1]
#    data = np.array(transforms.idct(data),dtype=data.dtype)
#    # keep zero-freq term first, reverse rest of spectrum and sign alt.
#    data = np.roll(data,-1,axis=-1)[...,::-1]
#    data[...,1::2] = data[...,1::2]*-1.
#
#    return data

def ps(data,p0=0.0,p1=0.0,inv=False):
    """ 
    linear phase correction
    
    Note: p0, p1 are in degrees
    """

    p0 = p0*pi/180. # convert to radians
    p1 = p1*pi/180. 

    size = data.shape[-1]
    apod = np.exp(1.0j*(p0+(p1*np.arange(size)/size) ),sig=data.dtype)

    if inv:
        apod = 1/apod
    return apod*data

def ps_exp(data,p0=0.0,tc=0.0,inv=False):
    """
    exponential phase correction

    Note: inv does not work
    """

    p0 = p0*pi/180. # convert to radians

    size = data.shape[-1]
    apod = np.exp(1.0j*(p0*np.exp(-tc*np.arange(size)/size)),sig=data.dtype)

    if inv:
        apod = 1/apod

    return apod*data


def smo(data,n):
    """
    Smooth data - takes average of +/-n points
    """

    n  = int(n)
    s = data.shape[-1]

    # a is the accumulator
    a = np.copy(data)
    
    for i in range(1,n+1):
        a = a+rs(data,i)+ls(data,i)
   
    # divide the interior by 2*n+1 to get mean
    a[...,n:-n] = a[...,n:-n]/(2*n+1)

    # divide the left edges by 2n+1-i where i is the distance from the interior
    for i in range(1,n+1):
        a[...,n-i] = a[...,n-i] / (2.*n+1-i)

    # divide the right edge similarly
    for i in range(-n,0):
        a[...,i] = a[...,i] / (n-i)

    return a 

# --- below here old code






# FFT related functions

def decode_States_old(data,axis=-1):
        """ create complex array from States data 
        
        Removes imaginary data from data before creating States
        complex array

        """
        data = real(data)       # discard imaginary portion of the data array

        # swap axes if needed
        if axis != 0:
                data = data.swapaxes(axis,0)

        # create the complex matrix 
        complex= data[::2] + data[1::2]*1.j
        
        # return the array with the axes if needed
        if axis !=0:
                return complex.swapaxes(axis,0)
        else:
                return complex

def decode_TPPI_old(data,axis=-1):
    # XXX write
    pass

    return

def decode_Hypercomplex_old(data,axis=-1):
    # XXX write
    pass


def complexft_old(data,axis=-1):
    """ nmr fft on complex data """

    # XXX rewrite this to work on the -1 axis
    # this works best with the axis to fourier transform in the 0th axis
    if axis != 0:
        data = data.swapaxes(axis,0)

    rawft = np.fft.fft(data,axis=0)                # calculate the raw fft
    rawft = np.array(rawft,dtype=data.dtype)
    shift = np.fft.fftshift(rawft,axes=[0])        # shift to frequencies
        
    # the Nyquist component has been stored as y[0] and we need is to 
    # be y[-1] so we correct this
    # XXX array.roll when avialable will make this nicer (numpy 1.2!)

    shiftprime = np.append(shift[1:],shift[np.newaxis,0],axis=0)

    rev = shiftprime[::-1]      # reverse the spectrum

    # swap the axes back if needed
    if axis!=0:
        return rev.swapaxes(axis,0) 
    else:
        return rev              
