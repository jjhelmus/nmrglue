#! /usr/bin/env python
""" Create files for dx unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.dx(d, a)
pipe.write("dx1.glue", d, a, overwrite=True)
