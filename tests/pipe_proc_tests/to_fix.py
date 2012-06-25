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
