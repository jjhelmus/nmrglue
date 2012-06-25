#! /usr/bin/env python
""" Create files for ft unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a)
pipe.write("ft1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, real=True)
pipe.write("ft2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, inv=True)
pipe.write("ft3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, alt=True)
pipe.write("ft4.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, neg=True)
pipe.write("ft5.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, null=True)
pipe.write("ft6.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, bruk=True)
pipe.write("ft7.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ft(d, a, auto=True)
pipe.write("ft8.glue", d, a, overwrite=True)
