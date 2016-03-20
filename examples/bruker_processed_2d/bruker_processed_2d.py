#! /usr/bin/env python
"""
Read and Plot pre-processed data.
"""
import os
import nmrglue as ng
import matplotlib.pyplot as plt

# Find the data
file_path = os.path.dirname(__file__)
data_file = os.path.join(file_path, 'data', 'bruker_exp', '1', 'pdata', '1')

# Read the data
dic, data = ng.bruker.read_pdata(data_file)

udic = ng.bruker.guess_udic(dic, data, strip_fake=True)

# make ppm scales
uc_13c = ng.fileiobase.uc_from_udic(udic, dim=1)
ppm_13c = uc_13c.ppm_scale()
ppm_13c_0, ppm_13c_1 = uc_13c.ppm_limits()

uc_15n = ng.fileiobase.uc_from_udic(udic, dim=0)
ppm_15n = uc_15n.ppm_scale()
ppm_15n_0, ppm_15n_1 = uc_15n.ppm_limits()

plt.figure()
plt.contour(data, extent=(ppm_13c_0, ppm_13c_1, ppm_15n_0, ppm_15n_1))
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.show()
