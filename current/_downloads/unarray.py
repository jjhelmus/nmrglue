#! /usr/bin/env python
# Unarray 1D varian data creating directories based on array parameter found in
# procpar file.  Directories are named param_value.fid

import nmrglue as ng
import os
import shutil

# the directory name of the arrayed data
dirname = "arrayed_data.dir"

# read the fid file and procpar file
fdic,data = ng.varian.read_fid_lowmem(os.path.join(dirname,"fid"))
dic = ng.varian.read_procpar(os.path.join(dirname,"procpar"))

# create directory names based of procpar array values
dir_pre = dic["array"]["values"][-1]
dir_mid = "_"
dir_values = dic[dir_pre]["values"]
dir_post = ".fid"
dir_names = ["".join([dir_pre,dir_mid,v,dir_post]) for v in dir_values]

# update the file header (assumes equal number of blocks in each slice)
fdic["nblocks"] = int(round(fdic["nblocks"]/len(dir_names)))

# unarray the data
for i,dir in enumerate(dir_names):
    
    # create the directory
    print "Creating:",dir
    if os.path.exists(dir)==False:
        os.makedirs(dir)
    
    # copy procpar file to the new directory
    src = os.path.join(dirname,"procpar")
    dst = os.path.join(dir,"procpar")
    shutil.copyfile(src,dst)

    # slice the data 
    slice = data[i::len(dir_names)] 
    
    # filename
    fname = os.path.join(dir,"fid")

    # write out the file
    ng.varian.write_fid(fname,fdic,slice,overwrite=True)
