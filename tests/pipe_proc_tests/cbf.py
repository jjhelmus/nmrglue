#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.cbf(d,a)
pipe.write("cbf1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.cbf(d,a,last=30)
pipe.write("cbf2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.cbf(d,a,reg=slice(1,100))
pipe.write("cbf3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.cbf(d,a,slice=slice(9,20))
pipe.write("cbf4.glue",d,a,overwrite=True)
