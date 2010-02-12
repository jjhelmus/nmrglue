#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p


d,a = pipe.read("time_complex.fid")
d,a = p.mc(d,a,mode="pow")
pipe.write("mc.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mc(d,a,mode="mod")
pipe.write("mc2.glue",d,a,overwrite=True)
