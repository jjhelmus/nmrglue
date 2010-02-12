#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.coad(d,a,cList=[1,1],axis="x",time=True)
pipe.write("co1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.coad(d,a,cList=[1,0,-5,8],axis="y",time=True)
pipe.write("co2.glue",d,a,overwrite=True)
