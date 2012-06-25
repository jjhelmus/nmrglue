#! /usr/bin/env python
""" Create files for smo unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.smo(d, a)
pipe.write("smo1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.smo(d, a, n=5)
pipe.write("smo2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.smo(d, a, center=True)
pipe.write("smo3.glue", d, a, overwrite=True)
