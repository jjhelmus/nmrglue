#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="ri2c")
pipe.write("s1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="c2ri")
pipe.write("s2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="ri2rr")
pipe.write("s3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="exlr")
pipe.write("s4.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="rolr")
pipe.write("s5.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="swap")
pipe.write("s6.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="inv")
pipe.write("s7.glue",d,a,overwrite=True)
