import nmrglue as ng

# read in the sum data set
dic, data = ng.varian.read('.', fid_file='fid_sum', as_2d=True)

# set the spectral parameters
udic = ng.varian.guess_udic(dic, data)
udic[1]['size']     = 1500             ; udic[0]['size']     = 256
udic[1]['complex']  = True             ; udic[0]['complex']  = True
udic[1]['encoding'] = 'direct'         ; udic[0]['encoding'] = 'states'
udic[1]['sw']       = 50000.000        ; udic[0]['sw']       = 5000.0
udic[1]['obs']      = 125.690          ; udic[0]['obs']      = 50.648
udic[1]['car']      = 174.538 * 125.690; udic[0]['car']      = 119.727 * 50.648
udic[1]['label']    = 'C13'            ; udic[0]['label']    = 'N15'

# convert to NMRPipe format
C = ng.convert.converter()
C.from_varian(dic, data, udic)
pdic, pdata = C.to_pipe()

# write out the NMRPipe file
ng.pipe.write("test_sum.fid", pdic, pdata, overwrite=True)

# repeat for the difference data set
dic, data = ng.varian.read('.', fid_file='fid_dif', as_2d=True)
C = ng.convert.converter()
C.from_varian(dic, data, udic)
pdic, pdata = C.to_pipe()
ng.pipe.write("test_dif.fid", pdic, pdata, overwrite=True)
