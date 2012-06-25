#! /usr/bin/env python
""" Create files for integ unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.integ(d, a)
pipe.write("integ1.glue", d, a, overwrite=True)
