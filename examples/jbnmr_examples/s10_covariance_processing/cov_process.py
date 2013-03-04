import nmrglue as ng
import numpy as np

# open the data
dic, data = ng.pipe.read("test.ft")

# compute the covariance
C = np.cov(data.T).astype('float32')

# update the spectral parameter of the indirect dimension
dic['FDF1FTFLAG'] = dic['FDF2FTFLAG']
dic['FDF1ORIG'] = dic['FDF2ORIG']
dic['FDF1SW'] = dic['FDF2SW']
dic["FDSPECNUM"] = C.shape[1]

# write out the covariance spectrum
ng.pipe.write("test.ft2", dic, C, overwrite=True)
