#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.varian.read_fid("../common_data/2d_varian_tppi/fid")

# Set the parameters (many of these can be found in the procpar file
u = ng.varian.guess_udic(dic,data)
# Direct Dimsion          
u[1]['size']    = 1400      # this is the R|I size
u[1]['complex'] = True     
u[1]['encoding']= 'direct' 
u[1]['sw']      = 50000.0  
u[1]['obs']     = 125.681  
u[1]['car']     = 101.274*125.681  
u[1]['label']   = 'C13x'
#Indirect Dimension
u[0]['size']    = 600       # this should be the R+I size
u[0]['complex'] = False 
u[0]['encoding']= 'tppi'
u[0]['sw']      = 33333.333
u[0]['obs']     = 125.681
u[0]['car']     = 101.274*125.681
u[0]['label']   = 'C13y'           

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_varian(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("2d_pipe_tppi.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read("2d_pipe_tppi.fid")
pdic2,pdata2 = ng.pipe.read("../common_data/2d_varian_tppi/test.fid")
print ng.misc.pair_similar(pdic,pdata,pdic2,pdata2,verb=True)
