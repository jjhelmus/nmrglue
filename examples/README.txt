nmrglue examples
================

This directory contains a number of examples which show how nmrglue can be used
to modify, convert, process, analyze, and visualize NMR data.  See the online
documentation for a description of the various examples.  


Getting the data required for the examples
------------------------------------------

The NMR data required for a number of these examples is large and not included 
in the git repository.  

To obtain all of the required data for the examples:

* Download and unpack the test data archive (test_data_VERSION.zip) into
  a data directory, run the included conversion scripts (see README)

* Edit the ``make_links.sh`` script to point to this directory and run to make
  symbolic links to the test data.

* Download the ``all_none_test_example_data.zip`` file and unzip the archive in
  the example directory.

The data for a single test can be also be downloaded as `individual zip files
<http://code.google.com/p/nmrglue/downloads/list?q=label:Example-Data>`_.


Generating the example zip files
--------------------------------

To make the zip files of the various examples run the
``make_example_zip_files.com`` script.  This creates zip files in the 
``zip_files`` directory.  These can be uploaded to the nmrglue google code 
website by hand or using the googlecode_upload.py command line utility after 
the previous file is deleted or depreciated using a command such as:

googlecode_upload.py -s "Convert Agilent to NMRPipe 1D example data" \
-p nmrglue -l Example-Data zip_files/example_agilent2pipe_1d.zip

