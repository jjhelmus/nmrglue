nmrglue 0.1 
===========

This is the first public release of nmrglue.  Please feel free to test it and
provide feedback.  For install instructions please see INSTALL.txt

What is nmrglue?
----------------

nmrglue is a Python module accessing a number of common NMR file formats in
python as numpy array objects and a python dictionary of spectral parameters.
nmrglue provides functions to read, writing, process and convert Bruker, 
NMRPipe, Sparky, Varian, and its own internal file format.

Where to get it
---------------

Information on nmrglue can be found at code.google.com/p/nmrglue.  The website
contains links to detailed documentation, code examples, and links to sample
data.  

Requirements
------------

nmrglue requires python 2.6, numpy 1.3.0 or later, scipy 0.7 or later and 
h5py 1.2 or later if reading and writing to the glue format is desired. 
