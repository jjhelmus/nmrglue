#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("1D_time.fid")
d,a = p.ha(d,a)
pipe.write("ha1.glue",d,a,overwrite=True)

d,a = pipe.read("1D_time.fid")
d,a = p.ha(d,a,inv=True)
pipe.write("ha2.glue",d,a,overwrite=True)
