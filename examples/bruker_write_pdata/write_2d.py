"""
Writing 2D data processed with Nmrpipe back into TopSpin format

"""
import os
from os.path import dirname as up
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

# set data paths
DATA_DIR = os.path.join(up(up(up(os.path.abspath(__file__)))),
                        'data', 'bruker_2d')

# read in data processed using TOPSPIN
# this is required only to get access to the dictionary
# as a reference starting point to write the data back out
# in a way TopSpin can read it
tdic, _ = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '1'))

# read in data processed using nmrpipe
# alternately, this data can be processed with nmrglue itself
ndic, ndata = ng.pipe.read(os.path.join(DATA_DIR, 'test.ft2'))

# Update dictionary parameters to match the data
# Note that none of the other parameters will macth what is
# written in the procs files and the procs file will not correspond to
# the actual data
tdic['procs']['SI'] = 2048
tdic['proc2s']['SI'] = 1024

# this writes only a '2rr' file which is read back correctly in topspin
ng.bruker.write_pdata(DATA_DIR, tdic, ndata.real, pdata_folder=10,
                      write_procs=True, overwrite=True)

# Read the processed data back in for plotting
data = {}
_, data['Topspin'] = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '1'))
_, data['NMRPipe'] = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '10'))

# PLOT
fig, ax = plt.subplots(figsize=(10, 5), ncols=2)

for axis, (k, v) in zip(ax, data.items()):
    noise = np.std(v.real[-50:, -50:])
    base = np.average(v.real[-50:, -50:])
    clev = 50*noise * 2 ** np.arange(4)
    axis.contour(v.real-base, levels=clev, linewidths=0.5, colors='k')
    axis.set_title(k)
    axis.set_xlim(200, 800)
    axis.set_ylim(200, 800)
plt.savefig('Comparison-2D.png')
