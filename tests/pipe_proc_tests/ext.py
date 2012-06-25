#! /usr/bin/env python
""" Create files for ext unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, left=True, sw=True)
pipe.write("ext1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, right=True)
pipe.write("ext2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, mid=True)
pipe.write("ext3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, x1=1, xn=100)
pipe.write("ext4.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, xn=200, y1=50, yn=75)
pipe.write("ext5.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, x1=5, xn=200, pow2=True)
pipe.write("ext6.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ext(d, a, x1=5, xn=200, round=10)
pipe.write("ext7.glue", d, a, overwrite=True)
