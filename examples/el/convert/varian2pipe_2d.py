#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.varian.read_fid("../common_data/2d_varian/fid")

# Set the parameters (many of these can be found in the procpar file
u = ng.varian.guess_udic(dic,data)
# Direct Dimsion          
u[1]['size']    = 1500      # this is the R|I size
u[1]['complex'] = True     
u[1]['encoding']= 'direct' 
u[1]['sw']      = 50000.0  
u[1]['obs']     = 125.691  
u[1]['car']     = 55.0*125.691  
u[1]['label']   = '13C'    
#Indirect Dimension
u[0]['size']    = 332       # this should be the R+I size
u[0]['complex'] = True 
u[0]['encoding']= 'states'
u[0]['sw']      = 5555.556
u[0]['obs']     = 50.648
u[0]['car']     = 120.0*50.648
u[0]['label']   = '15N'            

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_varian(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("2d_pipe.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read("2d_pipe.fid")
pdic2,pdata2 = ng.pipe.read("../common_data/2d_varian/test.fid")
print ng.misc.pair_similar(pdic,pdata,pdic2,pdata2,verb=True)
