#! /usr/bin/env python

import nmrglue as ng
import matplotlib.pyplot as plt

# read in the data
data_dir = "/Users/kjf/git/nmrglue/data/bruker_processed_1d/1/pdata/1"
dic, data = ng.bruker.read_pdata(data_dir)

udic = ng.bruker.guess_udic(dic, data)
uc = ng.fileiobase.uc_from_udic(udic)
ppm_scale = uc.ppm_scale()

# plot the spectrum
fig = plt.figure()
plt.plot(ppm_scale, data)
plt.show()
