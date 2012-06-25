#! /usr/bin/env python
""" Create files for zd unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.zd(d, a, wide=1.0, x0=20, slope=2, func=0)
pipe.write("zd1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.zd(d, a, wide=5.0, x0=10, slope=5, func=1)
pipe.write("zd2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.zd(d, a, wide=3.0, x0=100, slope=1, func=2)
pipe.write("zd3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.zd(d, a, wide=8.0, x0=15, slope=3, func=3, g=20)
pipe.write("zd4.glue", d, a, overwrite=True)
