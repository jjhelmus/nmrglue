"""
A large collection of NMR spectral processing functions which operate on the 
last (-1) 2D arrays.  These functions are wrapped by other processing modules.
All parameter are assumed to be in units of points unless otherwise noted.

"""

import numpy as np
import scipy
import scipy.signal

pi = np.pi


#########################
# Apodization functions #
#########################

def em(data,lb=0.0,inv=False):
    """
    Exponential Apodization

    Functional form:
        em(x_i) = exp(-pi*i*lb)

    Parameters:

    * data  Array of spectral data.
    * lp    Exponential line broadening.
    * inv   Set True for inverse apodization.

    """
    apod = np.exp(-pi*np.arange(data.shape[-1])*lb,sig=data.dtype)
    if inv:
        apod = 1/apod   # invert apodization

    return apod*data

def gm(data,g1=0.0,g2=0.0,g3=0.0,inv=False):
    """ 
    Lorentz-to-Gauss Apodization
    
    Functional form:
        gm(x_i) = exp(e - g*g)
    
    Where:  e = pi*i*g1
            g = 0.6*pi*g2*(g3*(size-1)-i)

    Parameters:

    * data  Array of spectral data.
    * g1    Inverse exponential width.
    * g2    Gaussian broaden width.
    * g3    Location of gauss maximum.
    * inv   Set True for inverse apodization.

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
    Modified Gaussian Apodization 

    Functional form:
        gmb(x_i) = exp(-a*i - b*i*i)

    Parameters:

    * data  Array of spectral data.
    * a     Exponential term.
    * b     Gaussian term.
    * inv   Set True for inverse apodization.

    returns ndarray object
    """ 
    size = data.shape[-1]
    apod = np.exp(-a*np.arange(size)-b*np.arange(size)**2,sig=data.dtype)
    if inv:
        apod = 1/apod

    return apod*data

def jmod(data,e=0.0,off=0.0,end=0.0,inv=False):
    """ 
    Exponentially Damped J-Modulation Apodization

    Functional form:
        jmod(x_i) = exp(-e)*sin( pi*off + pi*(end-off)*i/(size-1))

    Parameters:

    * data  Array of spectral data.
    * e     Exponential term.
    * off   Start of  modulation in fraction of pi radians.
    * end   End of modulation in fraction of pi radians.
    * inv   Set True for inverse apodization

    """
    size = data.shape[-1]
    apod = ( np.exp(-e*np.arange(size),sig=data.dtype) * np.sin(pi*off+pi*
           (end-off)*np.arange(size)/(size-1),sig=data.dtype) )
    if inv:
        apod = 1/apod

    return apod*data

def sp(data,off=0,end=1.0,pow=1.0,inv=False):
    """ 
    Shifted Sine-Bell Apodization

    Functional form:
            
        sp(x_i) = sin( (pi*off + pi*(end-off)*i/(size-1) )**pow

    Parameters:

    * data  Array of spectral data.
    * start Start of Sine-Bell as percent of vector (0.0 -> 1.0)
    * end   End of Sine-Bell as percent of vector (0.0 -> 1.0 )
    * pow   Power of Sine-Bell
    * inv   Set True for inverse apodization.

    """
    size = data.shape[-1]

    apod = np.power(np.sin(pi*off+pi*(end-off)*np.arange(size)/(size-1),
           sig=data.dtype),pow,sig=data.dtype)
    if inv:
        apod = 1/apod
    return apod*data

sine = sp

def tm(data,t1=0.0,t2=0.0,inv=False):
    """
    Trapezoid Apodization

    Functional form:
        0:t1        linear increase from 0.0 to 1.0
        t1:size-t2  1.0
        -t2:        linear decrease from 1.0 to 0.0
    
    Parameters:

    * data  Array of spectral data.
    * t1    Left ramp length.
    * t2    Right ramp length.
    * inv   Set True for inverse apodization.
    
    """
    size = data.shape[-1] 
    apod = np.concatenate( (np.linspace(0,1,t1),np.ones(size-t1-t2),
                            np.linspace(1,0,t2)) )
    apod = np.array(apod,dtype=data.dtype)  # type cast array
    if inv:
        apod = 1/apod

    return apod*data

def tri(data,loc="auto",lHi=0.0,rHi=0.0,inv=False):
    """
    Triangle Apodization

    Functional form:
        0:loc   linear increase/decrease from lHi to 1.0   
        loc:    linear increase/decrease from 1.0 to rHi

    Parameters:

    * data  Array of spectral data.
    * loc   Location of apex, "auto" sets to middle.
    * lHi   Left side starting height.
    * rHi   Right side starting height.
    * inv   Set True for inverse apodization.

    """
    size = data.shape[-1]
    if loc=="auto":
        loc = int(size/2.)
    apod = np.concatenate( (np.linspace(lHi,1.0,loc),
                            np.linspace(1.0,rHi,size-loc+1)[1:] ) )
    apod = np.array(apod,dtype=data.dtype)  # type cast
    if inv:
        apod = 1/apod

    return apod*data

###################
# Shift functions #
###################

def rs(data,pts=0.0):
    """
    Right Shift and Zero Fill
    
    Use roll or cs when zero fill not desired.
    
    Parameters:
    
    * data  Array of spectral data.
    * pts   Number of points to right shift.

    """
    data = np.roll(data,int(pts),axis=-1)
    data[...,:int(pts)] = 0
    return data

def ls(data,pts=0.0):
    """
    Left shift and fill with zero
    
    Use roll or cs when zero fill not desired.

    Parameters:

    * data  Array of spectral data.
    * pts   Number of points to right shift.

    """
    data = np.roll(data,-int(pts),axis=-1)
    data[...,-int(pts):] = 0
    return data


def cs(data,pts=0.0,neg=False):
    """
    Circular Shift 

    Parameters:

    * data  Array of spectral data.
    * pts   Number of points to right (+) or left (-) shift.
    * neg   Set True to Negate shifted points

    """

    return roll(data,pts,neg)

def roll(data,pts=0.0,neg=False):
    """
    Roll Axis 

    Parameters:

    * data  Array of spectral data.
    * pts   Number of points to right (+) or left (-) shift.
    * neg   Set True to Negate shifted points
    
    """
    data = np.roll(data,int(pts),axis=-1)
    if neg:
        if pts > 0:
            data[...,:pts] = -data[...,:pts]
        else:
            data[...,pts:] = -data[...,pts:]

    return data

def fsh(data,pts):
    """ 
    Frequency Shift by Fourier Transform (ifft->phase->fft)

    Positive pts shift spectrum to the right, negative to the left.

    Parameters:

    * data  Array of spectral data.
    * pts   Number of points to frequency shift

    """
    s= float(data.shape[-1])

    # inverse fft -> first order phase correction -> fft
    # idata = icomplexft(data)
    # pdata = np.exp(-2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*icomplexft(data)
    # data = complexft(pdata)

    # inplace version 
    return fft_positive(np.exp(-2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*
           ifft_negative(data))

def fsh2(data,pts):
    """
    Frequency Shift by Fourier Transform (fft->phase->ifft)

    Positive pts shift spectrum to the right, negative to the left.
    Odd values of pts change sign of spectrum, even values keep it the same.

    Parameters:

    * data  Array of spectral data.
    * pts   Number of points to frequency shift

    """
    s= float(data.shape[-1])
    return fft_positive(np.exp(2.j*pi*pts*np.arange(s)/s,sig=data.dtype)*
           ifft_negative(data)) 


##############
# Transforms #
##############

def nmr_reorder(data):
    """ 
    Reorder spectrum after FT transform to NMR order (swap halves and reverse)
    """
    s = data.shape[-1]
    return np.append(data[...,int(s/2)::-1],data[...,s:int(s/2):-1],axis=-1)

def swap_halves(data):
    """ 
    Swap the halves of a spectrum
    """
    s = data.shape[-1]
    return np.append(data[...,int(s/2):],data[...,:int(s/2)],axis=-1)

# Fourier based Transforms
def rft(x):
    """ 
    Real Fourier Transform 
    """
    # XXX figure out what exactly this is doing...
    s = x.shape[-1]
    xp = np.zeros(x.shape,dtype="complex64")
    xp[...,1:s/2] = x[...,1:-1:2]+x[...,2::2]*1.j
    xp[...,0] = x[...,0]/2.
    xp[...,s/2] = x[...,-1]/2.
    return np.array(nmr_reorder(np.fft.fft(2*xp,axis=-1).real),dtype="float32")

def irft(xp):
    """
    Inverse Real Fourier transform 
    """
    # XXX figure out what exactly this is doing
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


# Fourier transforms (keep these ones)
def fft(data):
    """ 
    Fourier transform, NMR ordering of results.

    There are a number of definitions of the discrete Fourier transform
    the version used in this function is:

    A_k =  \sum_{m=0}^{n-1} a_m \exp\left\{-2\pi i{mk \over n}\right\}
       \qquad k = 0,\ldots,n-1.

    With the inverse DFT in the ifft function defined as:

    .. math::
    a_m = \frac{1}{n}\sum_{k=0}^{n-1}A_k\exp\left\{2\pi i{mk\over n}\right\}
       \qquad n = 0,\ldots,n-1.

    Two alternative definitions are also supported by nmrglue. First one in
    which both the sum in the fft and ifft are multiplied by
    :math: \frac{1}{\sqrt{n}} which results in a pair of transforms in which 
    the total power contained in the the signals perform and after the 
    transforms are equal.  This is the type transforms used in the 
    Rowland NMR Toolkit. This type of transform is performed by the `fft_norm` 
    and `ifft_norm` functions. 

    The second definition changes the sign of the exponent to be positive while
    keeping the normalization factors the same.  This type of transform is
    performed by the NMRPipe processing package and the functions 
    `fft_positive` and `ifft_positive`.

    """
    return np.fft.fftshift(np.fft.fft(data, axis = -1).astype(data.dtype), -1)

def fft_norm(data):
    """
    Fourier transform, total power preserved, NMR ordering of results

    This is similar to the transform performed by RNMRTK's FFT function
    """
    return fft(data) / np.sqrt(float(data.shape[-1]))

def fft_positive(data):
    """ 
    Fourier transform with positive exponential, NMR ordering of results

    This is similar to the transform performed by NMRPipe's FFT function
    """
    # a positive exponential is the same as a IFFT, but we need to undo
    # the 1/N scaling
    s = float(data.shape[-1])
    return np.fft.fftshift(np.fft.ifft(data, axis=-1).astype(data.dtype),-1)*s

def ifft(data):
    """ 
    Inverse fourier transform, NMR ordering of results.
    """
    return np.fft.ifft(np.fft.ifftshift(data, -1), axis=-1).astype(data.dtype)

def ifft_norm(data):
    """ 
    Inverse fourier transform, total power preserved, NMR ordering of results

    This is similar to the transform performed by RNMRTK's IFFT function.
    """
    return ifft(data) * np.sqrt(float(data.shape[-1]))

def ifft_positive(data):
    """
    Inverse fourier transform with positive exponential, NMR ordered results.

    This is similar to the transform performed by NMRPipe's FFT function with
    the -inv flag
    """
    # a inverse fft with positive exponential in the FFT definition is the
    # same as a FFT with negative exponentials, but with a 1/N scaling factor
    s = 1.0 / float(data.shape[-1])
    return np.fft.fft(np.fft.ifftshift(data, -1), axis=-1).astype(data.dtype)*s 

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
    """
    Returns Hadamard matrix of order n (a power of 2)

    Source: Scipy Dev Wiki
    Author: Ivo
    http://projects.scipy.org/scipy/attachment/ticket/675/hadamard.py
    """
   
    ord = np.log2(n)
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
    Hadamard Transform

    This function is very slow.  Implement a Fast Walsh-Hadamard Transform
    with sequency/Walsh ordering (FWHT_w) for faster tranforms.

    See:
    
    http://en.wikipedia.org/wiki/Walsh_matrix
    http://en.wikipedia.org/wiki/Fast_Hadamard_transform

    """
    # implementation is a proof of concept and EXTEMEMLY SLOW

    # determind the order and final size of input vectors
    ord = int(np.ceil(np.log2(data.shape[-1])))  # Walsh/Hadamard order 
    max = 2**ord

    # zero fill to power of 2
    pad = max - data.shape[-1]
    zdata = zf(data,pad)

    # Multiple each vector by the hadamard matrix
    nat = np.zeros(zdata.shape,dtype=zdata.dtype)
    H = hadamard(max)
    nat = np.dot(zdata,H)
    nat = np.array(nat,dtype=data.dtype)

    # Bit-Reversal Permutation 
    s = [int2bin(x,digits=ord)[::-1] for x in range(max)]
    brp = [bin2int(x) for x in s]
    brp_data = np.take(nat,brp,axis=-1)

    # Gray code permutation (bit-inverse)
    gp = gray(ord)
    gp_data = np.take(brp_data,gp,axis=-1)

    return gp_data

def ht(data,N=None):
    """
    Hilbert Transform

    Reconstruct Imaginary Data Via Hilbert Transform

    Parameter:
    
    * data  Array of spectral data.
    * N     Number of Fourier components.

    """
  
    # XXX come back and fix this when a sane version of scipy.signal.hilbert
    # is included with scipy 0.8
    
    # create an empty output array
    fac = N/data.shape[-1]
    z = np.empty(data.shape,dtype=(data.flat[0]+data.flat[1]*1.j).dtype)
    if data.ndim == 1:
        z[:] = scipy.signal.hilbert(data.real,N)[:data.shape[-1]]*fac
    else:
        for i,vec in enumerate(data):
            z[i] = scipy.signal.hilbert(vec.real,N)[:data.shape[-1]]*fac

    # correct the real data as sometimes it changes
    z.real = data.real

    return z

##########################
# Standard NMR Functions #
##########################

def di(data):
    """
    Delete Imaginary Data
    """
    return data.real

def ps(data,p0=0.0,p1=0.0,inv=False):
    """ 
    Linear Phase Correction

    Parameters:

    * data  Array of spectral data.
    * p0    Zero order phase in degrees.
    * p1    First order phase in degrees.
    * inv   Set True for inverse phase correction

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
    Exponential Phase Correction

    Parameters:

    * data  Array of spectral data.
    * p0    Zero order phase in degrees.
    * tc    Exponential decay constant.

    """
    p0 = p0*pi/180. # convert to radians
    size = data.shape[-1]
    apod = np.exp(1.0j*(p0*np.exp(-tc*np.arange(size)/size)),sig=data.dtype)

    if inv:
        apod = 1/apod
    return apod*data

def tp(data):
    """ 
    Transpose Data

    * data  Array of spectral data.
    * hyper Set True if hypercomplex data.

    """
    return data.transpose()

ytp = tp
xy2yx = tp


def tp_hyper(data):
    """ 
    Hypercomplex Tranpose 

    Use when both dimension are complex
    """
    return c2ri(ri2c(data).transpose())

def zf_inter(data,pts=1):
    """
    Zero Fill between points

    Parameters:
    * data    Array of spectral data.
    * pts     number zeros between points

    """

    size = list(data.shape)
    size[-1] = (pts+1)*size[-1]
    z = np.zeros(size,dtype=data.dtype)
    z[...,::pts+1] = data[...,:] 
    return z

def zf_pad(data,pad=0,mid=False):
    """
    Zero Fill with pad zeros

    Parameters:
    * data  Array of spectral data.
    * pad   Number of zeros to add.
    * mid   Set to True to zero fill in middle of data.

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

zf = zf_pad

def zf_double(data,n,mid=False):
    """
    Zero fill to double size n times

    Parameters:

    * data  Array of spectral data.
    * n     Number of times to double data
    * mid   Set to True to zero fill in middle of data.
    
    """
    return zf_pad(data,int((data.shape[-1]*2**n)-data.shape[-1]),mid)
    
def zf_size(data,size,mid=False):
    """
    Zero fill to given size

    Parameters:
    * data  Array of spectral data.
    * size  Final size of data.
    * mid   Set to True to zero fill in middle of data.

    """
    return zf_pad(data,int(size-data.shape[-1]),mid)

def zf_auto(data,mid=False):
    """ 
    Zero fill to next largest power of two

    Parameters:

    * data  Array of spectral data.
    * mid   Set to True to zero fill in middle of data.

    """
    zf_size(data,np.ceil(np.log(data.shape[-1])/np.log(2)),mid)

####################
# Basic Untilities #
####################

# Add Constant

def add(data,r=0.0,i=0.0,c=0.0):
    """
    Add Constant

    Parameters:
    
    * data  Array of spectral data.
    * r     Constant to add to read channel.
    * i     Constant to add to imaginary channel.
    * c     Constant to add to both channels.
    
    """
    data.real = data.real + r + c
    data.imag = data.imag + i + c
    return data

def add_ri(data):
    """ 
    Add real and imaginary data to real channel
    """
    return data.real+data.imag

# Derivative

def dx(data):
    """
    Derivative by central difference

    Edges are takes as difference between nearest points

    """
    z = np.empty_like(data)
    z[...,0] = data[...,1] - data[...,0]    # first point
    z[...,-1]= data[...,-1] - data[...,-2]  # last point
    z[...,1:-1]  = data[...,2:] - data[...,:-2] # interior 
    return z

# Extract Region

def ext(data,x0=None,xn=None,y0=None,yn=None):
    """
    Extract Region

    Parameters:
    
    * data  Array of spectral data.
    * x0    X-axis extract region start.
    * xn    X-axis extract region end.
    * y0    Y-axis extract region start.
    * yn    Y-axis extract region end.
    
    """
    return data[y1:yn,x1:xn]

def ext_left(data):
    """
    Extract left half of spectrum
    """
    return data[...,0:int(data.shape[-1]/2.)]

def ext_right(data):
    """
    Extract right half of spectrum
    """
    return data[...,int(data.shape[-1]/2.):]

def ext_mid(data):
    """
    Extract middle of spectrum
    """
    return data[...,int(data.shape[-1]*1./4.):int(shape[-1]*3./4.)]

# Integrate

def integ(data):
    """ 
    Integrate by cumulative sum
    """
    return np.cumsum(data,axis=-1)

# Modulus Calculation

def mc(data):
    """ 
    Modulus calculation

    Calculated sqrt(real^2+imag^2)
    """
    return np.sqrt(data.real**2+data.imag**2)

def mc_pow(data):
    """
    Modulus calculation (square of)

    Calculated as: real^2+imag^2
    """
    return data.real**2+data.imag**2

# Mirror

def mir_left(data):
    """
    Append Mirror to left 
    """
    return np.append(data,data[...,::-1],axis=-1)

def mir_right(data):
    """
    Append Mirror to right
    """
    return np.append(data[...,::-1],data,axis=-1)

def mir_center(data):
    """
    Append Mirror to center
    """
    s = data.shape[-1] 
    return np.concatenate((data[...,s/2:],data,data[...,:s/2]),axis=-1)

def mir_center_onepoint(data):
    """
    Append mirror to center with one point shift (negate appended imag data)
    """
    s = int(data.shape[-1])
    data =  np.concatenate( (data[...,s-1:0:-1],data),axis=-1)
    if np.iscomplexobj(data):
        data.imag[...,:s-1] = -data.imag[...,:s-1]
    return data

# Multiply by a constant

def mult(data,r=1.0,i=1.0,c=1.0):
    """
    Multiply by a Constant

    Parameters:
    
    * data  Array of spectral data.
    * r     Constant to multiply real channel by.
    * i     Constant to multiply imaginary channel by.
    * c     Constant to multiply both channels by.
    
    """
    data.real = data.real*r*c
    if np.iscomplexobj(data):
        data.imag = data.imag*i*c
    return data

# Reverse

def rev(data):
    """
    Reverse Data
    """
    return data[...,::-1]

# Set to a Constant

def set(data,c):
    """
    Set Data to a Constant

    Parameters:
        
    * data  Array of spectral data.
    * c     Constant to set data to (may be complex)

    """
    data[...,:]=c
    return data

def set_complex(data,v):
    """
    Set real and imaginary portions of data to value

    Parameters:
    
    * data  Array of spectral data.
    * v     Constant to set data to (real number)
    
    """
    data.real = v
    if np.iscomplexobj(data):
        data.imag=v
    return data

def set_real(data,v):
    """
    Set real portion of data to a constant.

    Parameters:

    * data  Array of spectral data.
    * v     Constant to set data to (real number).

    """
    data.real = v
    return data

def set_imag(data,v):
    """
    Set imaginary portion of data to a constant.

    Parameters:

    * data  Array of spectral data.
    * v     Constant to set data to (real number).

    """
    if np.iscomplexobj(data):
        data.imag=v
    return data

# Shuffle Utilities

def ri2c(data):
    """
    Interleave real and imaginary data into a real array
    """
    s = list(data.shape)
    s[-1] = s[-1]*2
    n = np.empty(s,data.real.dtype)
    n[...,::2] = data.real
    n[...,1::2]= data.imag
    return n

def interleave_complex(data):
    """ 
    Unpack complex data into interleaved real,imaginary array
    """
    return ri2c(data)

def unpack_complex(data):
    """
    Unpacks complex array into real array (interleaves values)
    """
    return ri2c(data)    

def c2ri(data):
    """
    Seperate interleaved real,imaginary data into complex array

    Assumes data is real only, ignores imaginary portion of data

    """
    # make a 1,1 array to determind dtype
    temp = np.array(data.flat[0]+data.flat[1]*1j)
    s = list(data.shape)
    s[-1] = int(s[-1]/2)
    n = np.empty(s,temp.dtype)
    del(temp)
    n.real = data.real[...,::2]
    n.imag = data.real[...,1::2]
    return n

def seperate_interleaved(data):
    """
    Seperate interleaved real,imaginary data into complex array
    """
    return c2ri(data)

def pack_complex(data):
    """
    Packs interleaved real array into complex array
    """
    return c2ri(data)

def decode_States(data):
    """
    Decode data collected using States (seperate interleaved data)
    """
    return c2ri(data)


def ri2rr(data):
    """
    Append imaginary data to end of real data, returning a  real array
    """
    s = list(data.shape)
    half = int(s[-1])
    s[-1] = half*2
    n = np.empty(s,data.real.dtype)
    n[...,:half] = data.real
    n[...,half:] = data.imag
    return n

append_imag = ri2rr

def rr2ri(data):
    """
    Unappend real and imaginary data returing a complex array
    """
    # make a 1,1 array to determind dtype
    temp = np.array(data.flat[0]+data.flat[1]*1.j)
    s = list(data.shape)
    half = int(s[-1] / 2.0)
    s[-1] = half
    n = np.empty(s,temp.dtype)
    del(temp)
    n.real = data[...,:half]
    n.imag = data[...,half:]
    return n

unappend_imag = rr2ri

def exlr(data):
    """
    Exchange left and right halves of array
    """
    half = int(data.shape[-1]/2)
    n = np.empty_like(data)
    n[...,:half]=data[...,half:]
    n[...,half:]=data[...,:half]
    return n

exchange_lr = exlr

def rolr(data):
    """
    Rotate left and right halves of array
    """
    half = int(data.shape[-1]/2)
    n = np.empty_like(data)
    n[...,:half] = data[...,(half-1)::-1]
    n[...,half:] = data[...,:(half-1):-1]
    return n

rotate_lr = rolr

def swap(data):
    """ 
    Swap real and imaginary data
    """
    n = np.empty_like(data)
    n.real = data.imag
    n.imag = data.real
    return n

swap_ri = swap

def bswap(data):
    """
    Byteswap data
    """
    return data.byteswap()

byte_swap = bswap

# Sign Manipulation Utilities

def neg_left(data):
    """
    Negate left half
    """
    data[...,0:int(data.shape[-1]/2.)] = -data[...,0:int(data.shape[-1]/2.)]
    return data

def neg_right(data):
    """ 
    Negate right half
    """
    data[...,int(data.shape[-1]/2.):] = -data[...,int(data.shape[-1]/2.):]
    return data

def neg_middle(data):
    """
    Negate middle half
    """
    data[...,int(data.shape[-1]*1./4.):int(shape[-1]*3./4.)] = ( 
    -data[...,int(data.shape[-1]*1./4.):int(shape[-1]*3./4.)] )
    return data

def neg_edges(data):
    """
    Negate edge half of spectra
    """
    data[...,:int(data.shape[-1]*1./4)]= -data[...,:int(data.shape[-1]*1./4)]
    data[...,int(data.shape[-1]*3./4):]= -data[...,int(data.shape[-1]*3./4):]
    return data

def neg_all(data):
    """
    Negate data
    """
    return -data

def neg_real(data):
    """
    Negate real data
    """
    data.real = -data.real
    return data

def neg_imag(data):
    """
    Negate imaginary data
    """
    data.imag = -data.imag
    return data

def neg_even(data):
    """
    Negate even points
    """
    data[...,::2] = -data[...,::2]
    return data

def neg_odd(data):
    """
    Negate odd points
    """
    data[...,1::2] = -data[...,1::2]
    return data

def neg_alt(data):
    """
    Negate alternate (odd) points
    """
    return neg_odd(data)

def abs(data):
    """
    Replace data with absolute value of data (abs of real,imag seperately)
    """
    data.real = np.abs(data.real)
    data.imag = np.abs(data.imag)
    return data

def sign(data):
    """
    Replace data with sign (-1 or 1) of data (seperately on each channel)
    """
    data.real = np.sign(data.real)
    data.imag = np.sign(data.imag)
    return data

##################
# Misc Functions #
##################

# Coadd data

def coadd(data,clist,axis=-1):
    """
    Coadd data

    Reduce data along axis by applying blocks after multiplying by 
    coefficients in clist.  Incomplete blocks are discarded.   

    Parameters:
    
    * data  Array of spectral data.
    * clist List of Coefficients
    * axis  Axis to reduce (0=y,1=-1=x)
    
    """

    # there is probably a efficient way to do this with tile and inner
    # or scipy.ndimage.generic_filter

    # algorith creates a empty array, then fills it element wise
    # with each factor from clist times the blocks selected

    s = list(data.shape)    # data shape
    k = len(clist)          # length of coefficient list

    if axis == 1 or axis == -1:   # 'x' axis
        s[-1] = np.floor(float(s[-1])/k)
        n = np.zeros(s,dtype=data.dtype)
        m = s[-1] * k   # last element read
        for i in range(k):
            n = n + clist[i]*data[...,i:m:k]
    else:   # 'y' axis
        s[0] = np.floor(float(s[0])/k)
        n = np.zeros(s,dtype=data.dtype)
        m = s[0] * k
        for i in range(k):
            n = n + clist[i]*data[i:m:k]

    return n
        
coad = coadd

# Image Processing

def thres(data,thres=0.0):
    """
    Mark values less than thres as invalid (for use with filters)

    Parameter:

    * data  Array of spectral data.
    * thres Threshold value.

    """
    return np.ma.masked_less(data,thres)

def conv(data,kern=[1.],m="wrap",c=0.0):
    """
    Convolute data with kernel

    Parameters:

    * data  Array of spectral data.
    * kern  list or array describing convolution kernel
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode.

    """
    kern = np.array(kern)
    data.real = scipy.ndimage.convolve(data.real,weights=kern,mode=w,cval=c)
    data.imag = scipy.ndimage.convolve(data.imag,weights=kern,mode=w,cval=c)
    return data

convolute = conv

def corr(data,kern=[1.],m="wrap",c=0.0):
    """
    Correlate data with kernel (weights)

    Parameters:

    * data  Array of spectral data.
    * kern  list or array describing correlation kernel (weights)
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode.

    """
    kern = np.array(kern)
    data.real = scipy.ndimage.correlate(data.real,weights=kern,mode=w,cval=c)
    data.imag = scipy.ndimage.correlate(data.imag,weights=kern,mode=w,cval=c)
    return data

correlate = corr

def filter_median(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply median filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    data.real=scipy.ndimage.median_filter(data.real,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.median_filter(data.imag,size=s,mode=m,cval=c)
    return data

def filter_min(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply minimum filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    data.real=scipy.ndimage.minimum_filter(data.real,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.mimimum_filter(data.imag,size=s,mode=m,cval=c)
    return data

def filter_max(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply maximum filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    data.real=scipy.ndimage.maximum_filter(data.real,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.maximum_filter(data.imag,size=s,mode=m,cval=c)
    return data

def filter_percentile(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply a percentile filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    data.real=scipy.ndimage.percentile_filter(data.real,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.percentile_filter(data.imag,size=s,mode=m,cval=c)
    return data

def filter_rank(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply a rank filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    data.real=scipy.ndimage.rank_filter(data.real,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.rank_filter(data.imag,size=s,mode=m,cval=c)
    return data

# These filter are much slower as they use the generic filter functions...

def filter_amin(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply absolute minimum filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = amin_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_amax(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply absolute maximum filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = amax_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_range(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply range filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = range_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_avg(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply average filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = avg_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_dev(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply Standard Deviation filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = std_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_sum(data,s=(1,1),m="wrap",c=0.0):
    """
    Apply summation filter to data (real and imaginary seperately)

    Parameters:

    * data  Array of spectral data.
    * s     tuple defining shape or size taken for each step of the filter.
    * m     Defines how edges are determinded ('reflect','constant','nearest',
            'mirror','wrap').  Filter mode parameter.
    * c     Constant Value for use in 'constant' mode

    """
    flt = sum_flt
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

def filter_generic(data,filter,s=(1,1),m="wrap",c=0.0):
    """
    Apply generic filter to data (real and imaginary seperately)

    Parameters:

    * data   Array of spectral data.
    * filter python function which takes an array and returns a single value.
    * s      tuple defining shape or size taken for each step of the filter.
    * m      Defines how edges are determinded ('reflect','constant','nearest',
             'mirror','wrap').  Filter mode parameter.
    * c      Constant Value for use in 'constant' mode

    """
    flt = filter
    data.real=scipy.ndimage.generic_filter(data.real,flt,size=s,mode=m,cval=c) 
    data.imag=scipy.ndimage.generic_filter(data.imag,flt,size=s,mode=m,cval=c)
    return data

# filter functions

def amin_flt(arr): return  arr[np.abs(arr).argmin()]

def amax_flt(arr): return arr[np.abs(arr).argmax()]

def range_flt(arr): return arr.max() - arr.min()

def avg_flt(arr): return arr.avg()

def std_flt(arr): return arr.std()

def sum_flt(arr): return arr.sum()

# Scale Quad Artifacts

def qart(data,a=0.0,f=0.0):
    """
    Scale Quad Artifacts

    Functional Form:
        R' = R
        I' = (1+a)*I + f*R

    Parameters:

    * data  Array of spectral data.
    * a     Amplitude adjustment.
    * f     Phase adjustment.

    """
    data.imag = (1+a)*data.imag + f*data.real
    return data

def qart_auto(data):
    """
    Scale Quad artifacts by values from Gram-Schmidt orthogonalization
    """
    a,f = gram_schmidt(data)
    return qart(data,a,f)

def gram_schmidt(data):
    """
    calculate Gram-Schmidt orthogonalization parameters.
    """
    # method similar to Hock and Stern, "NMR Data Processing" p.61
    
    # sum of correlation between data.real and data.imag
    C = (data.real*data.imag).sum()

    # total power in real channel
    R = (data.real*data.real).sum()

    # remove correlation from imag channel
    idata = data.imag-(C/R)*data.real

    # total power in uncorrelated imag channel
    S = (idata*idata).sum()

    # imag(data'') = R/S*imag(data')
    # imag(data')  = imag(data)-C/R * real(data)
    # therefore:
    # imag(data'') = R/S*imag(data) - R*C/(S*R) * real(data) 
    # so A = R/S, B=-C/(S)
    return(R/S,-C/S)

# Complex Mixing

def qmix(data,carr):
    """
    Mix input and output channels provided coefficient array

    Parameters:
        
    * data  Array of spectral data.
    * carr  Array of coefficients for mixing
    
    """
    carr = np.array(carr).transpose()
    ic = carr.shape[1]  # input channels
    oc = carr.shape[0]  # output channels

    if data.shape[0]%ic!=0 or data.shape[0]%oc!=0:
        raise ValueError("Coefficient array does not evenly divide data")

    # create an empty blank output array
    s = list(data.shape)
    s[0] = s[0]*float(oc)/float(ic)
    n = np.empty(s,data.dtype)

    # remix each block
    for i in xrange(int(data.shape[0]/float(ic))):
        block = data[i*ic:(i+1)*ic]
        n[i*oc:(i+1)*oc]=np.dot(carr,block)

    return n

# Smooth and Center

def smo(data,n):
    """
    Smooth data

    Parameters:

    * data  Array of spectral data.
    * n     Smoothing window (+/- n points)
    """

    # XXX this can probably be accomplished by a median_filter
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

def center(data,n):
    """ 
    Center data

    Parameters:

    * data  Array of spectral data.
    * n     Centering window (+/- n points)

    """
    return data-smo(data,n)

def zd(data,window,x0=0.0,slope=1.0):
    """
    Zero Diagonal band with generic window functions

    Parameters:

    * data   Array of spectral data.
    * window Window function to apply to diagonal band
    * wide   Diagonal band half width in points.
    * x0     Diagonal starting location in points.
    * slope  Diagonal slope.
    
    """
    width = len(window)     # full width
    wide = (width-1.)/2     # half width        
    rows = data.shape[0]    # rows in data
    cols = data.shape[-1]   # columns in data
    c_start = x0+slope      # start of center diagonal band

    # last row to apply window to is last row or where we run off the grid
    max_r=int( min(rows, np.floor( (cols-c_start+wide)/slope)+1) )

    # apply window to band row by row
    for r in xrange(max_r): # r from 0 to max_r-1
        w_min = 0           # window min
        w_max = width       # window max

        c_mid = int(r*slope+(c_start))        # middle of diagonal band
        c_min = c_mid-wide
        c_max = c_mid+wide+1

        if c_min < 0:
            w_min = int(-c_min)
            c_min = 0
        if c_max > cols:
            w_max = int(w_max-(c_max-cols))
            c_max = cols

        w_size = c_max - c_min
        data[r,c_min:c_max] = data[r,c_min:c_max] * window[w_min:w_max]

    return data

def zd_boxcar(data,wide=1.0,x0=0.0,slope=1.0):
    """
    Zero Diagonal band with boxcar function

    Parameters:

    * data  Array of spectral data.
    * wide  Diagonal band half width in points.
    * x0    Diagonal starting location in points.
    * slope Diagonal slope.
    
    """
    window = np.zeros(2*wide+1)
    return zd(data,window,x0=x0,slope=slope)

def zd_triangle(data,wide=1.0,x0=0.0,slope=1.0):
    """
    Zero Diagonal band with triangle function

    Parameters:

    * data  Array of spectral data.
    * wide  Diagonal band half width in points.
    * x0    Diagonal starting location in points.
    * slope Diagonal slope.
    
    """
    window = np.append(np.linspace(1,0,wide+1),np.linspace(0,1,wide+1)[1:])
    return zd(data,window,x0=x0,slope=slope)

def zd_sinebell(data,wide=1.0,x0=0.0,slope=1.0):
    """
    Zero Diagonal band with sinebell function

    Parameters:

    * data  Array of spectral data.
    * wide  Diagonal band half width in points.
    * x0    Diagonal starting location in points.
    * slope Diagonal slope.
    
    """
    window = 1-np.sin(np.linspace(0,pi,2*wide+1))
    return zd(data,window,x0=x0,slope=slope)

def zd_gaussian(data,wide=1.0,x0=0.0,slope=1.0,g=1):
    """
    Zero Diagonal band with gaussian function

    Parameters:

    * data  Array of spectral data.
    * wide  Diagonal band half width in points.
    * x0    Diagonal starting location in points.
    * slope Diagonal slope.
    * g     Gauss width.
    
    """
    tln2 = np.sqrt(2*np.log(2))
    window = 1-scipy.signal.gaussian(2*wide+1,g/tln2)
    return zd(data,window,x0=x0,slope=slope)
