#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.varian.read_fid("../common_data/1d_varian/fid")

# Set the parameters (many of these can be found in the procpar file
u = ng.varian.guess_udic(dic,data)
u[0]['size']     = 1500
u[0]['complex']  = True
u[0]['encoding'] = 'direct'
u[0]['sw']       = 50000.0
u[0]['obs']      = 125.681
u[0]['car']      = 99.0*125.681
u[0]['label']    = 'C13'

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_varian(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("1d_pipe.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read("1d_pipe.fid")
pdic2,pdata2 = ng.pipe.read("../common_data/1d_varian/test.fid")
print ng.misc.pair_similar(pdic,pdata,pdic2,pdata2,verb=True)
