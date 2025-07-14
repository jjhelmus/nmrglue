#! /usr/bin/env python
# Create a interactive contour plot of a 2D NMRPipe spectrum

import nmrglue as ng
import matplotlib.pyplot as plt
import matplotlib.cm

# plot parameters
cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)
contour_start = 30000           # contour level start value
contour_num = 20                # number of contour levels
contour_factor = 1.20           # scaling factor between contour levels

# calculate contour levels
cl = [contour_start * contour_factor ** x for x in range(contour_num)]

# read in the data from a NMRPipe file
dic, data = ng.pipe.read("nmrpipe_2d/test.ft2")

# make ppm scales
uc_13c = ng.pipe.make_uc(dic, data, dim=1)
ppm_13c_0, ppm_13c_1 = uc_13c.ppm_limits()
uc_15n = ng.pipe.make_uc(dic, data, dim=0)
ppm_15n_0, ppm_15n_1 = uc_15n.ppm_limits()

# create the figure
fig = plt.figure()
ax = fig.add_subplot(111)

# plot the contours
ax.contour(data, cl, cmap=cmap,
        extent=(ppm_13c_0, ppm_13c_1, ppm_15n_0, ppm_15n_1))

# decorate the axes
ax.set_ylabel("15N (ppm)")
ax.set_xlabel("13C (ppm)")
ax.set_title("Protein 2D NCa Spectrum")
ax.set_xlim(70, 40)
ax.set_ylim(135, 100)

# start interactive session, script ends when window is closed
plt.show()
