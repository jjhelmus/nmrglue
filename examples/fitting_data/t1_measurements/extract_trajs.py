#! /usr/bin/env python
# Script to extract trajectories from a series a 2D spectrum.

import nmrglue as ng
import numpy as np

# read in the integration limits and list of spectra
peak_list = np.genfromtxt("boxes.in", dtype=None, encoding="UTF-8")
spectra_list = np.genfromtxt("spectra.in", dtype=None, encoding="UTF-8")

# prepare the trajs records array
num_spec = spectra_list.size
num_peaks = peak_list.size
elist = [np.empty(num_spec, dtype="float") for i in range(num_peaks)]
trajs = np.rec.array(elist, names=list(peak_list.f0))

# loop over the spectra
for sn, spectra in enumerate(spectra_list):

    # read in the data from a NMRPipe file
    print("Extracting from:", spectra)
    dic, data = ng.pipe.read(spectra)

    # loop over the integration limits
    for name, x0, y0, x1, y1 in peak_list:

        # integrate the region and save in trajs record array
        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0
        trajs[name][sn] = data[y0:y1 + 1, x0:x1 + 1].sum()

# normalize each trajectory
for peak in trajs.dtype.names:
    trajs[peak] = trajs[peak] / trajs[peak].max()

# save the trajectories records array to disk
np.save("traj.npy", trajs)
