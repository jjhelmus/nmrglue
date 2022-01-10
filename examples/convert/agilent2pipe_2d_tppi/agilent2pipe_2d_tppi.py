#! /usr/bin/env python

import nmrglue as ng

# read in the Agilent data
dic,data = ng.varian.read("agilent_2d_tppi")

# Set the spectral parameters.
u = ng.varian.guess_udic(dic,data)

# Direct Dimsion                        #Indirect Dimension
u[1]['size']     = 1400                 ; u[0]['size']     = 600
u[1]['complex']  = True                 ; u[0]['complex']  = False 
u[1]['encoding'] = 'direct'             ; u[0]['encoding'] = 'tppi'
u[1]['sw']       = 50000.0              ; u[0]['sw']       = 33333.333
u[1]['obs']      = 125.681              ; u[0]['obs']      = 125.681
u[1]['car']      = 101.274 * 125.681    ; u[0]['car']      = 101.274 * 125.681
u[1]['label']    = 'C13x'               ; u[0]['label']    = 'C13y'           

# create the converter object and initialize with Agilent data
C = ng.convert.converter()
C.from_varian(dic, data, u)

# create NMRPipe data and then write it out.
ng.pipe.write("2d_pipe_tppi.fid", *C.to_pipe(), overwrite=True)
