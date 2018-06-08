TEST
=======
nmrglue 
=======

What is nmrglue?
----------------

nmrglue is a module for working with NMR data in Python. When used with the
numpy, scipy, and matplotlib packages nmrglue provides a robust interpreted
environment for processing, analyzing, and inspecting NMR data.

Important Links
---------------

* Landing page: http://www.nmrglue.com
* Documentation: http://nmrglue.readthedocs.org/en/latest/index.html
* Examples: http://nmrglue.readthedocs.org/en/latest/examples/index.html
* Mailing List: https://groups.google.com/forum/#!forum/nmrglue-discuss
* Source code: https://github.com/jjhelmus/nmrglue

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

To get started, see our online documentation:
http://nmrglue.readthedocs.org/en/latest/index.html

Citing nmrglue
--------------

If you find nmrglue useful in your research please cite the package as:

J.J. Helmus, C.P. Jaroniec, Nmrglue: An open source Python package for the
analysis of multidimensional NMR data, J. Biomol. NMR 2013, 55, 355-367,
http://dx.doi.org/10.1007/s10858-013-9718-x.
