#! /usr/bin/env python

import nmrglue as ng

# read in the varian data (any of the follow lines will work)
dic,data=ng.bruker.read_lowmem("../common_data/3d_bruker")

# Set the parameters (many of these can be found in the procpar file
u = ng.bruker.guess_udic(dic,data)
# Direct Dimension
u[2]['size']     = 768
u[2]['complex']  = True
u[2]['encoding'] = 'direct'
u[2]['sw']       = 11061.947
u[2]['obs']      = 800.134
u[2]['car']      = 4.784*800.134
u[2]['label']    = '1H'

# First indirect dimension
u[1]['size']     = 128
u[1]['complex']  = True
u[1]['encoding'] = 'states'
u[1]['sw']       = 2500.000
u[1]['obs']      = 81.086
u[1]['car']      = 119.787*81.086
u[1]['label']    = '15N'

# Second indirect dimension
u[0]['size']     = 116
u[0]['complex']  = True
u[0]['encoding'] = 'states'
u[0]['sw']       = 5555.556
u[0]['obs']      = 201.204
u[0]['car']      = 55.743*201.204
u[0]['label']    = '13C'


# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_bruker(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("./data/3d_pipe%03d.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe (only check the 5th slice of 3D)
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read_lowmem("./data/3d_pipe%03d.fid")
pdic2,pdata2 = ng.pipe.read_lowmem("../common_data/3d_bruker/fid/test%03d.fid")
print ng.misc.pair_similar(pdic,pdata[5],pdic2,pdata2[5],verb=True)
