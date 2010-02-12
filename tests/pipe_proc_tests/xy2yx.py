#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.xy2yx(d,a,auto=True)
pipe.write("xy.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.xy2yx(d,a,hyper=True)
pipe.write("xy2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.xy2yx(d,a,nohyper=True)
pipe.write("xy3.glue",d,a,overwrite=True)
