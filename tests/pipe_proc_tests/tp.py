#! /usr/bin/env python
""" Create files for tp unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, auto=True)
pipe.write("tp1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, hyper=True)
pipe.write("tp2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, nohyper=True)
pipe.write("tp3.glue", d, a, overwrite=True)
