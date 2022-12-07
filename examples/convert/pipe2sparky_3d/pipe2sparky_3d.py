#! /usr/bin/env python

import nmrglue as ng

# read in the NMRPipe data (either line will work)
dic, data = ng.pipe.read("nmrpipe_3d/ft/test%03d.ft3")
#dic, data = ng.pipe.read_lowmem("nmrpipe_3d/ft/test%03d.ft3")

# Set the spectral parameters
udic = ng.pipe.guess_udic(dic, data)

# create the converter object and initialize with NMRPipe data
C = ng.convert.converter()
C.from_pipe(dic, data, udic)

# create Sparky data and then write it out
ng.sparky.write("sparky_3d.ucsf", *C.to_sparky(), overwrite=True)
