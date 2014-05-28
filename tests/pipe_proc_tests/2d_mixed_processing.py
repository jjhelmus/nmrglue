#! /usr/bin/env python
""" Create files for unit_conversion unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# process 2D mixed mode data
d, a = pipe.read("time_real.fid")
d, a = p.gmb(d, a, gb=0.1, lb=-8, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, alt=True) # BUG glue seems to double the data...?
#d, a = p.ps(d, a, p0=0, p1=0)
#d, a = p.tp(d, a, hyper=True)
#d, a = p.sp(d, a, off=0.5, pow=2, c=0.5)
#d, a = p.zf(d, a, auto=True)
#d, a = p.ft(d, a, auto=True)
#d, a = p.ps(d, a, p0=0, p1=0)
#d, a = p.di(d, a)
pipe.write("2d_mixed_processing1.glue", d, a, overwrite=True)

# process 2D mixed mode data
d, a = pipe.read("time_real.fid")
d, a = p.em(d, a, lb=8)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, auto=True)
d, a = p.di(d, a)
d, a = p.tp(d, a)
d, a = p.sp(d, a, off=0.5, pow=1, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, auto=True)
d, a = p.mc(d, a)
pipe.write("2d_mixed_processing2.glue", d, a, overwrite=True)
