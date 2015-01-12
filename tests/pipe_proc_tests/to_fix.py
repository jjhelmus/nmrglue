#! /usr/bin/env python
""" Tests which show a differences between NMRPipe's and nmrglue's processing
functions and a fix is desired.  """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# FSH
# the first three tests will fail because MIN/MAX values are off a little
d, a = pipe.read("1D_freq_real.dat")
d, a = p.fsh(d, a, dir="ls", pts=1)
pipe.write("fsh1.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.fsh(d, a, dir="rs", pts=8)
pipe.write("fsh2.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.fsh(d, a, dir="rs", pts=6.7)
pipe.write("fsh3.glue", d, a, overwrite=True)

# this fails because NMRPipe performs a Hilbert transform? If this is true it
# should be moved to the known_fail.py script.
d, a = pipe.read("1D_freq_complex.dat")
d, a = p.fsh(d, a, dir="ls", pts=9.5)
pipe.write("fsh4.glue", d, a, overwrite=True)

# RFT
# these are off by small amounts, mostly min/max values
d, a = pipe.read("1D_freq_complex.dat")
d, a = p.rft(d, a)
pipe.write("rft9.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.rft(d, a)
pipe.write("rft10.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.rft(d, a)
pipe.write("rft11.glue", d, a, overwrite=True)

# HT
# ps90-180 mode doesn't match
d, a = pipe.read("1D_freq_real.dat")
d, a = p.ht(d, a, mode="ps90-180")
pipe.write("ht4.glue", d, a, overwrite=True)

# Integration tests
# process 2D mixed mode data
d, a = pipe.read("time_real.fid")
d, a = p.gmb(d, a, gb=0.1, lb=-8, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, alt=True)  # BUG glue seems to double the data...?
d, a = p.ps(d, a, p0=0, p1=0)
d, a = p.tp(d, a, hyper=True)
d, a = p.sp(d, a, off=0.5, pow=2, c=0.5)
d, a = p.zf(d, a, auto=True)
d, a = p.ft(d, a, auto=True)
d, a = p.ps(d, a, p0=0, p1=0)
d, a = p.di(d, a)
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
