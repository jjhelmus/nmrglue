#! /usr/bin/env python
""" Create files for unit_conversion unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# process 2D complex data with Gaussian and cos^2 windows
d, a = pipe.read("time_complex.fid")
d, a = p.gmb(d, a, gb=0.1, lb=-8, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, auto=True)
d, a = p.ps(d, a, p0=-36, p1=0)
d, a = p.di(d, a)
d, a = p.tp(d, a)
d, a = p.sp(d, a, off=0.5, pow=2, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, auto=True)
d, a = p.ps(d, a, p0=-7, p1=0)
d, a = p.di(d, a)
pipe.write("2d_complex_processing1.glue", d, a, overwrite=True)
