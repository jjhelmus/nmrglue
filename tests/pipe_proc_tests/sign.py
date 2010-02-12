#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,ri=True)
pipe.write("s1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,r=True)
pipe.write("s2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,i=True)
pipe.write("s3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,left=True)
pipe.write("s4.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,right=True)
pipe.write("s5.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,alt=True)
pipe.write("s6.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,abs=True)
pipe.write("s7.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.sign(d,a,sign=True)
pipe.write("s8.glue",d,a,overwrite=True)

