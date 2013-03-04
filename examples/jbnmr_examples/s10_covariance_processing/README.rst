Covariance Processing example
=============================

Introduction
------------

This example is taken from Listing S10 in the 2013 JBNMR nmrglue paper.  In
this example covaraince processing is performed on a 2D NMRPipe file.


Instructions
------------

The `test.ft` file is provided in the archive.  To create this from time
domain data use the NMRPipe script `x.com`

Execute `python cov_process.py` to perform the covariance processing on the
`test.ft` file.  The file `test.ft2` is created.

Execute `python cov_plot.py` to perform the covariance processing and plot the
results.  The output, `covariance_figure.png` is presented as Figure 5 in the
article.
