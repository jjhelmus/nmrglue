#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# EM Testing
d,a = pipe.read("time_complex.fid")
d,a = p.em(d,a,lb=40,c=0.5,start=100,size=750,one=True)
pipe.write("em.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.em(d,a,lb=20,c=1.5,size=900)
pipe.write("em2.glue",d,a,overwrite=True)
