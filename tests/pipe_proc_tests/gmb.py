#! /usr/bin/env python
""" Create files for gmb unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.gmb(d, a, lb=2.0, gb=0.5, c=1.0, inv=True)
pipe.write("gmb1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.gmb(d, a, lb=10, gb=0.2, c=0.5, start=20, size=800)
pipe.write("gmb2.glue", d, a, overwrite=True)
