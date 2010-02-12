#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.jmod(d,a,cos=True,j=10.0,lb=5.0,c=1.0)
pipe.write("jmod.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.jmod(d,a,sin=True,j=18.0,lb=1.4,c=0.5,start=100,size=800,one=True)
pipe.write("jmod2.glue",d,a,overwrite=True)
