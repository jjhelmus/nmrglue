#! /usr/bin/env python
# Unarray 2D Agilent/Varian data creating directories based on array parameter
# found in procpar file.

import nmrglue as ng

# read in the NMR data
dic, data = ng.varian.read('arrayed_data.dir')

# set the new size of the separated data
dic['nblocks'] = data.shape[0]

arrayed_param = dic['procpar']['array']['values'][0]

# loop over the echo times, separating and saving each 2D
for i, array_val in enumerate(dic['procpar'][arrayed_param]['values']):
    dir_name = arrayed_param + '_' + array_val + '.fid'
    print "Creating directory:", dir_name
    ng.varian.write(dir_name, dic, data[:, i, :], overwrite=True)
