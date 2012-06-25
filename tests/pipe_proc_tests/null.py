#! /usr/bin/env python
""" Create files for null unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.null(d, a)
pipe.write("null1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.null(d, a)
d, a = p.di(d, a)
pipe.write("null2.glue", d, a, overwrite=True)
