#! /usr/bin/env python

import nmrglue as ng

# read in the file
dic, data = ng.pipe.read("nmrpipe_1d/test.fid")

# process the direct dimension
dic, data = ng.pipe_proc.sp(dic, data, off=0.35, end=0.98, pow=2, c=1.0)
dic, data = ng.pipe_proc.zf(dic, data, auto=True)
dic, data = ng.pipe_proc.ft(dic, data, auto=True)
dic, data = ng.pipe_proc.ps(dic, data, p0=-17.7, p1=-36.0)
dic, data = ng.pipe_proc.di(dic, data)

# write out processed data
ng.pipe.write("1d_pipe.ft", dic, data, overwrite=True)
