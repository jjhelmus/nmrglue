#! /usr/bin/env python

import nmrglue as ng

# read in the Bruker data
dic,data = ng.bruker.read("bruker_2d")

# Set the spectral parameters
u = ng.bruker.guess_udic(dic, data)

# Direct Dimsion                       #Indirect Dimension
u[1]['size']     = 768              ;  u[0]['size']     = 600
u[1]['complex']  = True             ;  u[0]['complex']  = True 
u[1]['encoding'] = 'direct'         ;  u[0]['encoding'] = 'states'
u[1]['sw']       = 11061.947        ;  u[0]['sw']       = 4000.000
u[1]['obs']      = 800.134          ;  u[0]['obs']      = 201.204
u[1]['car']      = 4.773 * 800.134  ;  u[0]['car']      = 58.742 * 201.204
u[1]['label']    = '1H'             ;  u[0]['label']    = '13C'            

# create the converter object and initilize with Bruker data
C = ng.convert.converter()
C.from_bruker(dic, data, u)

# create NMRPipe data and then write it out
ng.pipe.write("2d_pipe.fid", *C.to_pipe(), overwrite=True)
