#! /usr/bin/env python

import nmrglue as ng

# read in the varian data (any of the follow lines will work)
#dic,data=ng.varian.read_fid_lowmem_3D("../common_data/3d_varian/fid",(128,88))
#dic,data=ng.varian.read_fid("../common_data/3d_varian/fid",False,(128,88))
dic,data=ng.varian.read_fid_lowmem("../common_data/3d_varian/fid",(128,88))

# Set the parameters (many of these can be found in the procpar file
u = ng.varian.guess_udic(dic,data)
# Direct Dimension
u[2]['size']     = 1250
u[2]['complex']  = True
u[2]['encoding'] = 'direct'
u[2]['sw']       = 50000.0
u[2]['obs']      = 125.676
u[2]['car']      = 56.0*125.676
u[2]['label']    = 'CX'

# First indirect dimension
u[1]['size']     = 88
u[1]['complex']  = True
u[1]['encoding'] = 'states'
u[1]['sw']       = 2777.778
u[1]['obs']      = 50.648
u[1]['car']      = 120.0*50.648
u[1]['label']    = 'N'

# Second indirect dimension
u[0]['size']     = 128
u[0]['complex']  = True
u[0]['encoding'] = 'states'
u[0]['sw']       = 5555.556
u[0]['obs']      = 125.676
u[0]['car']      = 56.0*125.676
u[0]['label']    = 'CA'


# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_varian(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("./data/3d_pipe%03d.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe (only check the 5th slice of 3D)
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read_lowmem("./data/3d_pipe%03d.fid")
pdic2,pdata2 = ng.pipe.read_lowmem("../common_data/3d_varian/data/test%03d.fid")
print ng.misc.pair_similar(pdic,pdata[5],pdic2,pdata2[5],verb=True)
