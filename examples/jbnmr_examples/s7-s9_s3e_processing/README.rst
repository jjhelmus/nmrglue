Processing S3E filtered data example
====================================

Introduction
------------

This example is taken from Listing S7, S8 and S9 in the 2013 JBNMR nmrglue
paper.  In this example a 2D Agilent/Varian data set collect using a S3E filter
is separated (`seperate_s3e.py`), converted to NMRPipe format (`Sparky file
(`data.ucsf`) is converted to a NMRPipe file ('convert.py') and finally
processed (`xy_s3e.py`).


Instructions
------------

Execute `python seperate_s3e.py` to separate the S3E sum and difference
spectra from data set in the Agilent/Varian `fid` file.  This creates the files
fid_dif and fid_sum.

Execute `python convert.py` to convert these two files to NMRPipe format.  This
step creates the files `test_sum.fid` and `test_dif.fid`.

Execute `python xy_s3e.py` to process and combine the sum and different
spectra. This step creates the `test.ft2` file.
