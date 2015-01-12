#! /usr/bin/env python
""" Create files for sign unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, ri=True)
pipe.write("sign1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, r=True)
pipe.write("sign2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, i=True)
pipe.write("sign3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, left=True)
pipe.write("sign4.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, right=True)
pipe.write("sign5.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, alt=True)
pipe.write("sign6.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, abs=True)
pipe.write("sign7.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.sign(d, a, sign=True)
pipe.write("sign8.glue", d, a, overwrite=True)
