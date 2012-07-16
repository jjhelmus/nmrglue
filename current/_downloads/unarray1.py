#! /usr/bin/env python
# Unarray 2D bruker data into directories names 1,2...

import nmrglue as ng
import os
import shutil

# these parameters must determined by the user
shape = (7360,640) # shape of data
array_size = 23     # number of interleaved data sets

# the directory name of the arrayed data
dirname = "arrayed_data.dir"

import numpy as np
from nmrglue import *

# read the binary data and reshape
dic,data = ng.bruker.read_binary(os.path.join(dirname,"ser"),shape=shape)

for i in xrange(array_size):
    
    # create the directory
    dir = str(i+1)
    print "Creating:",dir
    if os.path.exists(dir)==False:
        os.makedirs(dir)

    # write the ser file
    fname = os.path.join(dir,"ser")
    ng.bruker.write_binary(fname,{},data[i::array_size],overwrite=True)

    # copy over additional files (adjust as needed)
    shutil.copy(os.path.join(dirname,"acqu"),os.path.join(dir,"acqu"))
    shutil.copy(os.path.join(dirname,"acqu2"),os.path.join(dir,"acqu2"))
    shutil.copy(os.path.join(dirname,"acqu2s"),os.path.join(dir,"acqu2s"))
    shutil.copy(os.path.join(dirname,"acqus"),os.path.join(dir,"acqus"))
    shutil.copy(os.path.join(dirname,"pulseprogram"),
                os.path.join(dir,"pulseprogram"))

