nmrglue examples
================

This directory contains a number of examples which show how nmrglue can be used
to modify, convert, process, analyze, and visualize NMR data.  See the online
documentation for a description of the various examples.  

Getting the data required for the examples
------------------------------------------

The NMR data required for a number of these examples is large and not included 
in the git repository.  All of the required data can be download XXX

The data for individual tests can be 


Generating the example zip files
--------------------------------

To make the zip files of the various example run the ``make_
This creates zip files in the ``zip_files`` directory.  These can be uploaded 
to the nmrglue google code website by hand or using the googlecode_upload.py
command line utility after the previous file is deleted or depreciated using a 
command such as:

googlecode_upload.py -s "Convert Agilent to NMRPipe 1D example data" \
-p nmrglue -l Example-Data zip_files/example_agilent2pipe_1d.zip

