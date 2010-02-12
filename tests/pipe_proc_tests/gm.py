#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.gm(d,a,g1=5.0,g2=2.0,g3=0.0,c=1.0)
pipe.write("gm.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.gm(d,a,g1=2.0,g2=5.0,g3=0.5,c=1.5,start=100,size=654,one=True)
pipe.write("gm2.glue",d,a,overwrite=True)


d,a = pipe.read("time_complex.fid")
d,a = p.gm(d,a,g1=2.0,g2=5.0,g3=0.5,c=1.5,start=100)
pipe.write("gm3.glue",d,a,overwrite=True)
