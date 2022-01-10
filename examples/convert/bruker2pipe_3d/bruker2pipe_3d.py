#! /usr/bin/env python

import nmrglue as ng

# read in the Bruker data
dic, data = ng.bruker.read("bruker_3d")
#dic, data = ng.bruker.read_lowmem("bruker_3d")

# Set the spectral parameters
udic = ng.bruker.guess_udic(dic, data)

# Direct Dimension
udic[2]['size']     = 768
udic[2]['complex']  = True
udic[2]['encoding'] = 'direct'
udic[2]['sw']       = 11061.947
udic[2]['obs']      = 800.134
udic[2]['car']      = 4.784 * 800.134
udic[2]['label']    = '1H'

# First indirect dimension
udic[1]['size']     = 128
udic[1]['complex']  = True
udic[1]['encoding'] = 'states'
udic[1]['sw']       = 2500.000
udic[1]['obs']      = 81.086
udic[1]['car']      = 119.787 * 81.086
udic[1]['label']    = '15N'

# Second indirect dimension
udic[0]['size']     = 116
udic[0]['complex']  = True
udic[0]['encoding'] = 'states'
udic[0]['sw']       = 5555.556
udic[0]['obs']      = 201.204
udic[0]['car']      = 55.743 * 201.204
udic[0]['label']    = '13C'


# create the converter object and initialize with Bruker data
C = ng.convert.converter()
C.from_bruker(dic, data, udic)

# create NMRPipe data and then write it out
ng.pipe.write("./data/3d_pipe%03d.fid", *C.to_pipe(), overwrite=True)
