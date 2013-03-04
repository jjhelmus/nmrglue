import nmrglue as ng
dic, data = ng.varian.read('.', as_2d=True)
dic['nblocks'] /= 2
A = data[::2]
B = data[1::2]
ng.varian.write_fid('fid_sum', dic, A + B, overwrite=True)
ng.varian.write_fid('fid_dif', dic, A - B, overwrite=True)
