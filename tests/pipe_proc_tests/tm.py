#! /usr/bin/env python
""" Create files for tm unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.tm(d, a, t1=100, t2=200, c=1.0, inv=True)
pipe.write("tm1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.tm(d, a, t1=10, t2=1400, c=1.5, start=5, size=1490)
pipe.write("tm2.glue", d, a, overwrite=True)
