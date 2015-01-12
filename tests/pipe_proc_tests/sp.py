#! /usr/bin/env python
""" Create files for sp unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.sp(d, a, off=0.35, end=0.98, pow=2.0, c=1.0)
pipe.write("sp1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sp(d, a, off=0.10, end=0.75, pow=1.0, c=0.5, size=200, one=True)
pipe.write("sp2.glue", d, a, overwrite=True)
