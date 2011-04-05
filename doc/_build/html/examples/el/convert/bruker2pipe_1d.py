#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.bruker.read("../common_data/1d_bruker/")

# Set the parameters (many of these can be found in the procpar file
u = ng.bruker.guess_udic(dic,data)
u[0]['size']     = 2048
u[0]['complex']  = True
u[0]['encoding'] = 'direct'
u[0]['sw']       = 10000.000
u[0]['obs']      = 600.133
u[0]['car']      = 4.773*600.133
u[0]['label']    = '1H'

# create the converter object and initilize with bruker data
C = ng.convert.converter()
C.from_bruker(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("1d_pipe.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read("1d_pipe.fid")
pdic2,pdata2 = ng.pipe.read("../common_data/1d_bruker/test.fid")
print ng.misc.pair_similar(pdic,pdata,pdic2,pdata2,verb=True)
