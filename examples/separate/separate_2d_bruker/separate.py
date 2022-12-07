#! /usr/bin/env python
# Separate 2D data sets from an arrayed data set, directories will be names
# 1, 2, 3, ... 23

import nmrglue as ng

# read in the NMR data
dic, data = ng.bruker.read('arrayed_data.dir', shape=(7360, 640), cplex=True)

array_size = 23

# loop over the arrayed data, separating and saving each 2D
for i in xrange(array_size):
    dir_name = str(i+1)
    print "Creating directory:", dir_name
    ng.bruker.write(dir_name, dic, data[i::array_size], overwrite=True)
