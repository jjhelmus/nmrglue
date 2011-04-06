nmrglue 0.2 
===========

What is nmrglue?
----------------

nmrglue is a module for working with NMR data in Python. When used with the
numpy, scipy, and matplotlib packages nmrglue provides a robust interpreted
environment for processing, analysing, and inspecting NMR data.

What can nmrglue do?
--------------------

nmrglue has the ability to read, write and convert between a number of common
NMR file formats including Varian, Bruker, NMRPipe, and Sparky files. The
files, which are represented in python as dictionaries of spectral parameters
and Numpy array objects, can be easily examined, modified and processed as
desired. 

nmrglue provides a number of common functions for processing NMR data such as
apodization, spectral shifting, Fourier and other transformations, baseline
smoothing and flattening, and linear prediction. In addition new processing
schemes can be implemented easily using the nmrglue provided functions and the
multitude of numerical routines provided by the Numpy and Scipy packages. 

When used in conjunction with the matplotlib (or other) python plotting
library nmrglue can be used to create publication quality figures of NMR
spectrum or examine data interactively.

nmrglue can be used to analysis NMR data, with routines to perform peak
picking, multidimensional lineshape fitting (peak fitting), and peak
integration. New analysis methods can be rapidly developed and tested in
Python or by integrating Fortran and C/C++ code.

Citing nmrglue
--------------

The article describing nmrglue is still in preparation. For the time being
please cite nmrglue as: 

J. J. Helmus and C.P. Jaroniec, nmrglue, http://code.google.com/p/nmrglue, The
Ohio State University.
