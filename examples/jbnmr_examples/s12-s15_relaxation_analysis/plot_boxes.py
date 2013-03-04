import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
import matplotlib.cm

# plot parameters
xpad = 5                        # padding around peak box on x-axis
ypad = 5                        # padding around peak box on y-axis
cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)

# contour levels
cl = 30000 * 1.20 ** np.arange(20)

# read in the box limits and list of spectra
peak_list = np.recfromtxt("boxes.in", names=True)
spectra_list = np.recfromtxt("spectra.in")

# loop over the spectra 
for spec_number, spectra in enumerate(spectra_list):

    # read in the spectral data
    dic, data = ng.pipe.read(spectra)

    # loop over the peaks
    for peak, x0, y0, x1, y1 in peak_list:

        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0

        # slice the data around the peak
        slice = data[y0 - ypad:y1 + 1 + ypad, x0 - xpad:x1 + 1 + xpad]

        # create the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # plot the contours
        print "Plotting:", peak, spec_number
        extent = (x0 - xpad + 1, x1 + xpad - 1, y0 - ypad + 1, y1 + ypad - 1)
        ax.contour(slice, cl, cmap=cmap, extent=extent)

        # draw a box around the peak
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k--')

        # draw lighter boxes at +/- 1 point
        ax.plot([x0 - 1, x1 + 1, x1 + 1, x0 - 1, x0 - 1],
                [y0 - 1, y0 - 1, y1 + 1, y1 + 1, y0 - 1], 'k--', alpha=0.35)
        ax.plot([x0 + 1, x1 - 1, x1 - 1, x0 + 1, x0 + 1],
                [y0 + 1, y0 + 1, y1 - 1, y1 - 1, y0 + 1], 'k--', alpha=0.35)

        # set the title, save the figure
        ax.set_title('Peak: %s Spectrum: %i'%(peak, spec_number))
        fig.savefig('peak_%s_spectrum_%i'%(peak, spec_number))
        del(fig)
