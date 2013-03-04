import nmrglue as ng
import numpy as np

# read the NMR data, forcing the data to be two dimensional
dic, data = ng.varian.read('arrayed_data.fid', as_2d=True)

# set the new size of the separated data
array_size = len(dic['procpar']['nredor']['values'])
out_shape = int(data.shape[0] / array_size), data.shape[1]
dic['nblocks'] = out_shape[0]

# loop over the redor multiples, separating and saving each 2D
for i, nredor in enumerate(dic['procpar']['nredor']['values']):
    dir_name = 'nredor_' + nredor + '.fid'
    print "Creating directory:", dir_name
    sdata = np.empty(out_shape, dtype=data.dtype)
    sdata[::2] = data[2 * i::2 * array_size]
    sdata[1::2] = data[2 * i + 1::2 * array_size]
    ng.varian.write(dir_name, dic, sdata, overwrite=True)
