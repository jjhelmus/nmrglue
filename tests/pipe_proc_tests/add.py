#! /usr/bin/env python
""" Create files for add unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, r=2.0, i=-1.0)
pipe.write("add1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, c=10.0)
pipe.write("add2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, c=10.0, x1=10, xn=400)
pipe.write("add3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, ri=True, x1=50, xn=300)
pipe.write("add4.glue", d, a, overwrite=True)
