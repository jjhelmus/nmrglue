import nmrglue as ng
import numpy as np

# read the integration limits and list of spectra
peak_list = np.genfromtxt("boxes.in", names=True)
spectra_list = np.genfromtxt("spectra.in")

# create an array to hold the trajectories
trajectories = np.empty((peak_list.size, spectra_list.size), dtype='float')

# loop over the spectra
for sn, spectra in enumerate(spectra_list):

    # read in the spectra data
    print "Extracting peak intensities from:", spectra
    dic, data = ng.pipe.read(spectra)

    # loop over the integration limits
    for i, (name, x0, y0, x1, y1) in enumerate(peak_list):

        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0

        # integrate the region and save in trajectories array
        trajectories[i][sn] = data[y0:y1 + 1, x0:x1 + 1].sum()

# write out the trajectories for each peak
for itraj, peak_traj in enumerate(trajectories):
    peak_traj /= peak_traj.max()    # normalize the peak's trajectory
    fname = peak_list.peak_label[itraj] + '.dat'
    f = open(fname, 'w')
    for v in peak_traj:
        f.write(str(v) + '\n')
    f.close()
