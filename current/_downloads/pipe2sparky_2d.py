#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.pipe.read("../common_data/2d_pipe/test.ft2")

# Set the parameters 
u = ng.pipe.guess_udic(dic,data)

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_pipe(dic,data,u)

# create pipe data and then write it out
ng.sparky.write("2d_sparky.ucsf",*C.to_sparky(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
sdic,sdata = ng.sparky.read("2d_sparky.ucsf")
sdic2,sdata2 = ng.sparky.read("../common_data/2d_sparky/data.ucsf")
print ng.misc.pair_similar(sdic,sdata,sdic2,sdata2,verb=True)
