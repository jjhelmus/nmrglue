#! /usr/bin/env python
# Create a 1D interactive plot from NMRPipe data

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

# read in the data from a NMRPipe file
dic,data = ng.pipe.read("../../common_data/1d_pipe/test.ft")

# determind the ppm scale
uc = ng.pipe.make_uc(dic,data)
ppm_scale = np.linspace(uc.ppm(0),uc.ppm(data.size-1),data.size)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ppm_scale,data,'k-')

# decorate axes
ax.set_yticklabels([])
ax.set_title("Protein 1D Spectrum")
ax.set_xlabel("13C ppm")
ax.set_xlim(200,0)
ax.set_ylim(-80000,2500000)

# start interactive session, script ends when windows closes
plt.show()
