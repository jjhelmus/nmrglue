#! /usr/bin/env python

import nmrglue as ng
import numpy as np
import glob
import os.path

# create a list of directories to coadd
dlist = glob.glob("run*.dir")
dlist.sort()

# create a list of 2D files in the first directory
flist = glob.glob(os.path.join(dlist[0], "*.fid", "test.fid"))
flist.sort()

# loop over the files
for base_fname in flist:

    # initialize the new data
    dic, data = ng.pipe.read(base_fname)
    coadd_data = np.zeros_like(data)
    coadd_dic = dict(dic)

    # loop over the pseudo-3D directories
    for d in dlist:

        # the file names is found by replace the directory name
        f = base_fname.replace(dlist[0], d, 1)
        print "Coadding file:", f
        dic, data = ng.pipe.read(f)
        coadd_data += data / len(dlist)

    # write out the file
    of = base_fname.replace(dlist[0], "coadded_data.dir", 1)
    print "Writing out:", of
    ng.pipe.write(of, coadd_dic, coadd_data, overwrite=True)
