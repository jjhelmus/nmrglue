#! /usr/bin/env python

import nmrglue as ng
import matplotlib.pyplot as plt

# read in the data
dic, data = ng.pipe.read('test.ft2')

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data[20000:25000])
fig.savefig('figure_nmrpipe.png')
