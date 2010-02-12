#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.tri(d,a,loc=500,lHi=0.5,rHi=0.8,inv=True)
pipe.write("tri.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.tri(d,a,loc=750,lHi=0.1,rHi=0.5)
pipe.write("tri2.glue",d,a,overwrite=True)



