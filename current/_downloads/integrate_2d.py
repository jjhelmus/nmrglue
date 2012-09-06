#! /usr/bin/env python
# Example script to show integration of a 2D spectrum

import nmrglue as ng
import numpy as np

# read in the data from a NMRPipe file
dic, data = ng.pipe.read("nmrpipe_2d/test.ft2")

# read in the integration limits
peak_list = np.recfromtxt("limits.in")

# prepare the output file
f = open("volumes.out",'w')
f.write("# Name\tVolume\n")

# loop over the integration limits
for name, x0, y0, x1, y1 in peak_list:

    if x0 > x1:
        x0, x1 = x1, x0
    if y0 > y1:
        y0, y1 = y1, y0

    vol = data[y0:y1 + 1, x0:x1 + 1].sum()
    f.write("%s\t%.3f\n"%(name, vol))

# close the output file
f.close()
