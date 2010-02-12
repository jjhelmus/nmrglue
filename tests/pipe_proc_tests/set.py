#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.set(d,a,r=1.0,i=2.0)
pipe.write("set1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.set(d,a,c=3.0)
pipe.write("set2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.set(d,a,r=8.0,i=-2.0,c=10.0)
pipe.write("set3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.set(d,a,r=1.0,i=2.0,x1=50,xn=700)
pipe.write("set4.glue",d,a,overwrite=True)
