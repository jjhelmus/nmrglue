import nmrglue as ng

# read in the NMR data
dic, data = ng.varian.read('arrayed_data.fid')

# set the new size of the separated data
dic['nblocks'] = data.shape[0]

# loop over the echo times, separating and saving each 2D
for i, techo in enumerate(dic['procpar']['techo']['values']):
    dir_name = 'techo_' + techo + '.fid'
    print "Creating directory:", dir_name
    ng.varian.write(dir_name, dic, data[:, i, :], overwrite=True)
