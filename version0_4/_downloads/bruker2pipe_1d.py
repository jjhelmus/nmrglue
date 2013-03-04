#! /usr/bin/env python

import nmrglue as ng

# read in the Bruker data
dic, data = ng.bruker.read("bruker_1d")

# Set the spectral parameters.
udic = ng.bruker.guess_udic(dic, data)
udic[0]['size']     = 2048
udic[0]['complex']  = True
udic[0]['encoding'] = 'direct'
udic[0]['sw']       = 10000.000
udic[0]['obs']      = 600.133
udic[0]['car']      = 4.773 * 600.133
udic[0]['label']    = '1H'

# create the converter object and initilize with Bruker data
C = ng.convert.converter()
C.from_bruker(dic, data, udic)

# create NMRPipe data and then write it out
ng.pipe.write("1d_pipe.fid", *C.to_pipe(), overwrite=True)
