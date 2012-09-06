#! /usr/bin/env python

import nmrglue as ng

# read in the file
dic, data = ng.pipe.read("nmrpipe_2d_tppi/test.fid")

# process the direct dimension
dic, data = ng.pipe_proc.sp(dic, data, off=0.35, end=0.98, pow=2, c=1.0)
dic, data = ng.pipe_proc.zf(dic, data, auto=True)
dic, data = ng.pipe_proc.ft(dic, data, auto=True)
dic, data = ng.pipe_proc.ps(dic, data, p0=151.0, p1=0.0)
dic, data = ng.pipe_proc.di(dic, data)

# process the indirect dimension
dic, data = ng.pipe_proc.tp(dic, data)
dic, data = ng.pipe_proc.sp(dic, data, off=0.35, end=0.98, pow=2, c=0.5)
dic, data = ng.pipe_proc.zf(dic, data, auto=True)
dic, data = ng.pipe_proc.ft(dic, data, auto=True)
dic, data = ng.pipe_proc.ps(dic, data, p0=0.0, p1=0.0)
dic, data = ng.pipe_proc.di(dic, data)
dic, data = ng.pipe_proc.rev(dic, data, sw=True)
dic, data = ng.pipe_proc.tp(dic, data)

# write out processed data
ng.pipe.write("2d_pipe_tppi.ft2", dic, data, overwrite=True)
