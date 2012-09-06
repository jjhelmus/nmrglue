#! /usr/bin/env python

import nmrglue as ng
import numpy as np
import glob

# create a list of files to coadd
flist = glob.glob("test*.fid")
flist.sort()

# initilize the new data
dic, data = ng.pipe.read(flist[0])
coadd_data = np.zeros_like(data)
coadd_dic = dict(dic)

# loop over the files, adding them to the coadded array
for f in flist:
    print "Coadding file:", f
    dic, data = ng.pipe.read(f)
    coadd_data += data / len(flist)

# write out the file
print "Writing out coadded data"
ng.pipe.write("coadded.fid", coadd_dic, coadd_data, overwrite=True)
