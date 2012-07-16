#! /usr/bin/env python

# import the necessary functions
from nmrglue import *
import numpy as np
import glob

# create a list of files to coadd
flist = glob.glob("test*.fid")
flist.sort()

# initilize the new data
dic,data = pipe.read(flist[0])
coadd_data = np.zeros_like(data)
coadd_dic  = dict(dic)

# loop over files and add them coadded array
for f in flist:
    print "Reading file:",f
    dic,data = pipe.read(f)
    coadd_data = coadd_data + data/len(flist)

# write out the file
print "Writing out file"
pipe.write("coadded.fid",coadd_dic,coadd_data,True)
