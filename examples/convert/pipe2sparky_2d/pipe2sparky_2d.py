#! /usr/bin/env python

import nmrglue as ng

# read in the NMRPipe data
dic,data = ng.pipe.read("nmrpipe_2d/test.ft2")

# Set the spectral parameters
udic = ng.pipe.guess_udic(dic, data)

# create the converter object and initialize with NMRPipe data
C = ng.convert.converter()
C.from_pipe(dic, data, udic)

# create Sparky data and then write it out
ng.sparky.write("sparky_2d.ucsf", *C.to_sparky(), overwrite=True)
