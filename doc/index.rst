Overview
========

:Release: |version|
:Date: |today|



nmrglue is a module for working with NMR data in `Python <http://python.org>`_.
When used with the `Numpy <http://numpy.scipy.org/>`_, 
`Scipy <http://scipy.org/>`_, and 
`matplotlib <http://matplotlib.sourceforge.net/>`_ packages nmrglue provides a 
robust interpreted environment for processing, analysing, and inspecting NMR
data.  

What can nmrglue do?
--------------------

nmrglue has the ability to read, write and convert between a number of common 
NMR file formats including Varian, Bruker, NMRPipe, and Sparky files.  The 
files, which are represented in python as dictionaries of spectral parameters 
and `Numpy <http://numpy.scipy.org/>`_ array objects, can be easily examined, 
modified amd processed as desired.

nmrglue provides a number of common functions for processing NMR data such as
apodization, spectral shifting, Fourier and other transformations, 
baseline smoothing and flattening, and linear prediction.  In addition new 
processing schemes can be implemented easily using the nmrglue provided 
functions and the multitude of numerical routines provided by the 
`Numpy <http://numpy.scipy.org/>`_ and `Scipy <http://scipy.org/>`_ packages.

When used in conjunction with the 
`matplotlib <http://matplotlib.sourceforge.net/>`_ (or other) python plotting 
library nmrglue can be used to create publication quality figures of NMR 
spectrum or examine data interactively.  

.. image:: /examples/el/plotting/1d_spectrum/spectrum.png
    :scale: 50 %

.. image:: /examples/el/plotting/2d_spectrum/spectrum.png
    :scale: 50 %


nmrglue can be used to analysis NMR data, with routines to perform peak 
picking, multidimensional lineshape fitting (peak fitting), and 
peak intergration.  New analysis methods can be rapidly developed and tested
in Python or by integrating Fortran and C/C++ code.


documentation
-------------

.. toctree::
    :maxdepth: 2

    tutorial
    reference/index
    examples/index

