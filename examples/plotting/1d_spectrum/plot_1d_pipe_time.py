#! /usr/bin/env python
# Create a 1D plot of NMRPipe data

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

# read in the data from a NMRPipe file
dic,data = ng.pipe.read("../../common_data/1d_pipe/test.fid")

# determind the time scale
sw = dic["FDF2SW"]
times = np.array([x*1.e3/sw for x in range(data.size)])

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(times,data,'k-')

# decorate axes
ax.set_yticklabels([])
ax.set_title("Protein 1D FID")
ax.set_xlabel("Time (ms)")
ax.set_ylim(-100000,100000)

# save the figure
fig.savefig("fid.png") # change this to .pdf, .ps, etc
