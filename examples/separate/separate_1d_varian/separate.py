#! /usr/bin/env python
# Separate 1D Agilent/Varian data creating directories based on the
# array parameter found in the procpar file.

import nmrglue as ng

dic, data = ng.varian.read('arrayed_data.dir')

dic['nblocks'] = 1

arrayed_param = dic['procpar']['array']['values'][0]

# loop over the echo times, separating and saving each 1D
for i, array_val in enumerate(dic['procpar'][arrayed_param]['values']):
    dir_name = arrayed_param + '_' + array_val + '.fid'
    print "Creating directory:", dir_name
    ng.varian.write(dir_name, dic, data[i], overwrite=True)
