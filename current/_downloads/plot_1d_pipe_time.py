#! /usr/bin/env python
# Create a 1D plot of NMRPipe data

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

# read in the data from a NMRPipe file
dic,data = ng.pipe.read("nmrpipe_1d/test.fid")

# make a unit conversion object
uc = ng.pipe.make_uc(dic, data)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(uc.ms_scale(), data, 'k-')

# decorate axes
ax.set_yticklabels([])
ax.set_title("Protein 1D FID")
ax.set_xlabel("Time (ms)")
ax.set_ylim(-100000, 100000)

# save the figure
fig.savefig("fid.png") # this can be to .pdf, .ps, etc
