import os
from os.path import dirname as up
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt

# PATH TO DATA
# DATA_DIR = os.path.join(up(up(up(os.path.abspath(__file__)))), 'data', 'bruker_1d')
DATA_DIR = '/home/kaustubh/Desktop/bruker-testing/yellow/1'

# TOPSPIN
## read in data processed using tospin
dic, data = ng.bruker.read(DATA_DIR) 
tdic, (treal, timag) = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '1'),
                                            bin_files=['1r', '1i'])

# NMRPIPE
# read in data processed using nmrpipe
pdic, pdata = ng.pipe.read(os.path.join(DATA_DIR, 'test.ft1')) 

# write out the real and imag data in topspin format in pdata folder 10
for dataset, filename in zip([pdata.real, pdata.imag], ['1r', '1i']):
    ng.bruker.write_pdata(DATA_DIR, dic, dataset, bin_file=filename, 
                          pdata_folder=10, write_procs=True, overwrite=True)

# read back in the nmrpipe processed files for comparison
pdic, (preal, pimag) = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '10'),
                                            bin_files=['1r', '1i'])

# NMRGLUE 
## read in raw FID and process using nmrglue
ndic, ndata = ng.bruker.read(DATA_DIR)
ndata = ng.proc_base.zf_size(ndata, size=8192)
ndata = ng.proc_base.fft(ndata)
ndata = ng.bruker.remove_digital_filter(ndic, ndata, post_proc=True)
ndata = np.roll(ndata, 1)
ndata = ng.proc_base.rev(ndata)
ndata = ng.proc_base.ps(ndata, p0=101)

# write out the real and imag data in topspin format in pdata folder 11
for dataset, filename in zip([ndata.real, ndata.imag], ['1r', '1i']):
    ng.bruker.write_pdata(DATA_DIR, dic, dataset, bin_file=filename, 
                          pdata_folder=11, write_procs=True, overwrite=True)

## read back in the nmrglue processed files for comparison
ndic, (nreal, nimag) = ng.bruker.read_pdata(os.path.join(DATA_DIR, 'pdata', '11'),
                                            bin_files=['1r', '1i'])

# dataset and title lists
dataset = [treal, timag, preal, pimag, nreal, nimag]
titles = ['Topspin-REAL', 'Topspin-IMAG', 'Nmrpipe-REAL', 
          'Nmrpipe-IMAG', 'Nmrglue-REAL', 'Nmrglue-IMAG']

# PLOT
fig, ax = plt.subplots(figsize=(8, 5), nrows=3, ncols=2, sharex='col', sharey=False)
for axis, dataset, title in zip(ax.flat, dataset, titles): 
    axis.plot(dataset, linewidth=0.7)
    axis.set_title(title)
    axis.set_yticklabels([])
    # axis.set_xlim(1000, 5200)
plt.tight_layout()
plt.savefig('Comparison-1D.png')
