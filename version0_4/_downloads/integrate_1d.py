#! /usr/bin/env python
# Example scipt to show integration of a 1D spectrum

import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt

# read in the data from a NMRPipe file
dic, data = ng.pipe.read("1d_data.ft")
length = data.shape[0]

# read in the integration limits
peak_list = np.recfromtxt("limits.in")

# determind the ppm scale
uc = ng.pipe.make_uc(dic, data)
ppm_scale = uc.ppm_scale()

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ppm_scale, data, 'k-')

# prepare the output file
f = open("area.out", 'w')
f.write("#Name\tStart\tStop\tArea\n")

# loop over the integration limits
for name, start, end in peak_list:
    min = uc(start, "ppm")
    max = uc(end, "ppm")
    if min > max:
        min, max = max, min

    # extract the peak
    peak = data[min:max + 1]
    peak_scale = ppm_scale[min:max + 1]

    # plot the integration lines, limits and name of peaks
    ax.plot(peak_scale, peak.cumsum() / 100. + peak.max(), 'g-')
    ax.plot(peak_scale, [0] * len(peak_scale), 'r-')
    ax.text(peak_scale[0], 0.5 * peak.sum() / 100. + peak.max(), name,
                fontsize=8)

    # write out the integration info
    tup = (name, peak_scale[0], peak_scale[-1], peak.sum())
    f.write("%s\t%.3f\t%.3f\t%E\n" % tup)

# close the output file and save the plot
f.close()
ax.set_xlim(ppm_scale[0], ppm_scale[-1])
fig.savefig("plot.png")
