#! /usr/bin/env python

import nmrglue as ng

# read in the varian data
dic,data = ng.bruker.read("../common_data/2d_bruker/")

# Set the parameters (many of these can be found in the procpar file
u = ng.bruker.guess_udic(dic,data)
# Direct Dimsion          
u[1]['size']    = 768      # this is the R|I size
u[1]['complex'] = True     
u[1]['encoding']= 'direct' 
u[1]['sw']      = 11061.947  
u[1]['obs']     = 800.134
u[1]['car']     = 4.773*800.134  
u[1]['label']   = '1H'    
#Indirect Dimension
u[0]['size']    = 600       # this should be the R+I size
u[0]['complex'] = True 
u[0]['encoding']= 'states'
u[0]['sw']      = 4000.000
u[0]['obs']     = 201.204
u[0]['car']     = 58.742*201.204
u[0]['label']   = '13C'            

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_bruker(dic,data,u)

# create pipe data and then write it out
ng.pipe.write("2d_pipe.fid",*C.to_pipe(),overwrite=True)

# check the conversion against NMRPipe
print "Conversion complete, listing differences between files:"
pdic,pdata = ng.pipe.read("2d_pipe.fid")
pdic2,pdata2 = ng.pipe.read("../common_data/2d_bruker/test.fid")
print ng.misc.pair_similar(pdic,pdata,pdic2,pdata2,verb=True)
