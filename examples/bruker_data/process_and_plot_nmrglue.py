#! /usr/bin/env python

import nmrglue as ng
import matplotlib.pyplot as plt

# read in the bruker formatted data
dic, data = ng.bruker.read('expnmr_00001_1')

# remove the digital filter
data = ng.bruker.remove_digital_filter(dic, data)

# process the spectrum
data = ng.proc_base.zf_size(data, 32768)    # zero fill to 32768 points
data = ng.proc_base.fft(data)               # Fourier transform
data = ng.proc_base.ps(data, p0=-88.0)      # phase correction
data = ng.proc_base.di(data)                # discard the imaginaries
data = ng.proc_base.rev(data)               # reverse the data

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data[20000:25000])
fig.savefig('figure_nmrglue.png')
