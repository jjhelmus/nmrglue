#! /usr/bin/env python

import nmrglue as ng

# read in the Agilent data (any of the follow lines will work)
#dic, data=ng.varian.read("agilent_3d")
dic, data=ng.varian.read_lowmem("agilent_3d")

# Set the spectral parameters
udic = ng.varian.guess_udic(dic, data)

# Direct Dimension
udic[2]['size']     = 1250
udic[2]['complex']  = True
udic[2]['encoding'] = 'direct'
udic[2]['sw']       = 50000.0
udic[2]['obs']      = 125.676
udic[2]['car']      = 56.0 * 125.676
udic[2]['label']    = 'CX'

# First indirect dimension
udic[1]['size']     = 88
udic[1]['complex']  = True
udic[1]['encoding'] = 'states'
udic[1]['sw']       = 2777.778
udic[1]['obs']      = 50.648
udic[1]['car']      = 120.0 * 50.648
udic[1]['label']    = 'N'

# Second indirect dimension
udic[0]['size']     = 128
udic[0]['complex']  = True
udic[0]['encoding'] = 'states'
udic[0]['sw']       = 5555.556
udic[0]['obs']      = 125.676
udic[0]['car']      = 56.0 * 125.676
udic[0]['label']    = 'CA'


# create the converter object and initilize with Agilent data
C = ng.convert.converter()
C.from_varian(dic, data, udic)

# create NMRPipe data and then write it out
ng.pipe.write("./data/3d_pipe%03d.fid", *C.to_pipe(), overwrite=True)
