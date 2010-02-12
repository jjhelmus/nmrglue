#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.qart(d,a,a=1.0,f=0.5)
pipe.write("qart.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.qart(d,a,a=0.8,f=1.2)
pipe.write("qart2.glue",d,a,overwrite=True)
