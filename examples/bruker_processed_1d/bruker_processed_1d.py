#! /usr/bin/env python

import nmrglue as ng
import matplotlib.pyplot as plt

# read in the data
data_dir = "data/bruker_processed_1d/1/pdata/1"

# From pre-procced data.
dic, data = ng.bruker.read_pdata(data_dir, scale_data=True)

udic = ng.bruker.guess_udic(dic, data)
uc = ng.fileiobase.uc_from_udic(udic)
ppm_scale = uc.ppm_scale()

# From FID
dic1, data1 = ng.bruker.read(data_dir)

# remove the digital filter, this data if from an analog spectrometer.
# data = ng.bruker.remove_digital_filter(dic, data)

# process the spectrum
data1 = ng.proc_base.ls(data1, 1)             # left shift
data1 = ng.proc_base.gm(data1, g2=1/2.8e3)    # To match proc data
data1 = ng.proc_base.zf_size(data1, 1024*32)  # zero fill
data1 = ng.proc_base.fft(data1)               # FT
data1 = ng.proc_base.ps(data1, p0=93)         # phase is 180 off bruker
data1 = ng.proc_base.di(data1)                # discard
data1 = ng.proc_base.rev(data1)               # reverse the data

udic1 = ng.bruker.guess_udic(dic1, data1)
uc1 = ng.fileiobase.uc_from_udic(udic1)
ppm_scale1 = uc1.ppm_scale()

# plot the spectrum
fig = plt.figure()
plt.hold(True)
plt.plot(ppm_scale, data)
plt.plot(ppm_scale1, data1)
plt.hold(False)
plt.xlim([50, -50])
plt.show()
