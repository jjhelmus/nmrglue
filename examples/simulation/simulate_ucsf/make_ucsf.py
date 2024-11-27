#! /usr/bin/env python

import nmrglue as ng
import numpy as np

# create a sparky dictionary
# A dictionary from a existing Sparky ucsf file can be found using:
# ng.sparky.guess_udic(*ng.sparky.read('filename.ucsf'))
udic = {
    'ndim': 2,
    0: {'car': 7329.0,
        'complex': False,
        'encoding': 'states',
        'freq': True,
        'label': '15N',
        'obs': 60.8,
        'size': 512,
        'sw': 1523.43,
        'time': False},
    1: {'car': 5403.570418865944,
        'complex': False,
        'encoding': 'direct',
        'freq': True,
        'label': '1H',
        'obs': 600.0,
        'size': 1024,
        'sw': 3606.5,
        'time': False}
}

dic = ng.sparky.create_dic(udic)
data = np.empty((512, 1024), dtype='float32')

# read in the peak list
peak_list = np.genfromtxt('peaks.txt', names=True)
npeaks = len(peak_list)

# convert the peak list from PPM to points
uc_15N = ng.sparky.make_uc(dic, None, 0)
uc_1H = ng.sparky.make_uc(dic, None, 1)

lw_15N = 5.0    # 15N dimension linewidth in points
lw_1H = 5.0     # 1H dimension linewidth in points

params = []
for ppm_15N, ppm_1H in peak_list:
    pts_15N = uc_15N.f(ppm_15N, 'ppm')
    pts_1H = uc_1H.f(ppm_1H, 'ppm')
    params.append([(pts_15N, lw_15N), (pts_1H, lw_1H)])

# simulate the spectrum
shape = (512, 1024)      # size should match the dictionary size
lineshapes = ('g', 'g')  # gaussian in both dimensions
amps = [100.0] * npeaks
data = ng.linesh.sim_NDregion(shape, lineshapes, params, amps)

# save the spectrum
ng.sparky.write("test.ucsf", dic, data.astype('float32'), overwrite=True)
