#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.zf(d,a,zf=2)
pipe.write("zf.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.zf(d,a,pad=200)
pipe.write("zf2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.zf(d,a,size=8000)
pipe.write("zf3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.zf(d,a,size=4096,mid=True)
pipe.write("zf4.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.zf(d,a,size=4096,inter=True)
pipe.write("zf5.glue",d,a,overwrite=True)
