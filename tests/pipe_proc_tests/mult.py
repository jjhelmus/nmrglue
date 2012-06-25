#! /usr/bin/env python
""" tests for MULT function """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.mult(d, a, c=1.5, inv=True)
pipe.write("mult1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.mult(d, a, r=2.0, i=1.2)
pipe.write("mult2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.mult(d, a, r=1.5, i=0.9, x1=90, xn=600)
pipe.write("mult3.glue", d, a, overwrite=True)

