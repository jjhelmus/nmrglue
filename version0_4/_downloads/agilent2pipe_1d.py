#! /usr/bin/env python

import nmrglue as ng

# read in the Agilent data
dic, data = ng.varian.read('agilent_1d')

# Set the spectral parameters.
udic = ng.varian.guess_udic(dic, data)
udic[0]['size']     = 1500
udic[0]['complex']  = True
udic[0]['encoding'] = 'direct'
udic[0]['sw']       = 50000.0
udic[0]['obs']      = 125.681
udic[0]['car']      = 99.0*125.681
udic[0]['label']    = 'C13'

# create the converter object and load in Agilent data
C = ng.convert.converter()
C.from_varian(dic, data, udic)

# create NMRPipe data and then write it out
ng.pipe.write("1d_pipe.fid", *C.to_pipe(), overwrite=True)

