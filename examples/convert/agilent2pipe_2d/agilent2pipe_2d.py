#! /usr/bin/env python

import nmrglue as ng

# read in the Agilent data
dic, data = ng.varian.read("agilent_2d")

# Set the spectral parameters
udic = ng.varian.guess_udic(dic, data)

# Direct dimension                  # Indirect dimension
udic[1]['size']     = 1500          ; udic[0]['size']     = 332 
udic[1]['complex']  = True          ; udic[0]['complex']  = True
udic[1]['encoding'] = 'direct'      ; udic[0]['encoding'] = 'states'
udic[1]['sw']       = 50000.0       ; udic[0]['sw']       = 5555.556
udic[1]['obs']      = 125.691       ; udic[0]['obs']      = 50.648
udic[1]['car']      = 55.0 * 125.691; udic[0]['car']      = 120.0 * 50.648
udic[1]['label']    = '13C'         ; udic[0]['label']    = '15N'    

# create the converter object and initilize with Agilent data
C = ng.convert.converter()
C.from_varian(dic, data, udic)

# create NMRPipe data and then write it out
ng.pipe.write("2d_pipe.fid", *C.to_pipe(), overwrite=True)
