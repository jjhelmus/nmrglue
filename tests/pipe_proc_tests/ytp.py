#! /usr/bin/env python
""" Create files for ytp unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.ytp(d, a, auto=True)
pipe.write("ytp1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ytp(d, a, hyper=True)
pipe.write("ytp2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ytp(d, a, nohyper=True)
pipe.write("ytp3.glue", d, a, overwrite=True)
