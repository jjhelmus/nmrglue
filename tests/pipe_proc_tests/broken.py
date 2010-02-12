#! /usr/bin/env python

# A Collection of pipe_proc tests which do not match NMRPipe

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# FSH

# the first three tests will fail because MIN/MAX values are off a little
d,a = pipe.read("1D_freq_real.dat")
d,a = p.fsh(d,a,dir="ls",pts=1)
pipe.write("fsh1.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.fsh(d,a,dir="rs",pts=8)
pipe.write("fsh2.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.fsh(d,a,dir="rs",pts=6.7)
pipe.write("fsh3.glue",d,a,overwrite=True)

# this fails because NMRPipe performs a Hilbert transform?
d,a = pipe.read("1D_freq_complex.dat")
d,a = p.fsh(d,a,dir="ls",pts=9.5)
pipe.write("fsh4.glue",d,a,overwrite=True)

# RFT 
# these are off by small amounts, mostly min/max values
d,a = pipe.read("1D_freq_complex.dat")
d,a = p.rft(d,a)
pipe.write("rft9.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.rft(d,a)
pipe.write("rft10.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.rft(d,a)
pipe.write("rft11.glue",d,a,overwrite=True)

# HT
# ps90-180 mode doesn't match
d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a,mode="ps90-180")
pipe.write("ht4.glue",d,a,overwrite=True)

# ADD
# the following show the difference way NMRPipe and pipe_proc handle
# when c and r or i (or both) are defined.
d,a = pipe.read("time_complex.fid")
d,a = p.add(d,a,c=10.0,r=5.0,i=3.0)
pipe.write("add5.glue",d,a,overwrite=True)

# the following show the difference way NMRPipe and pipe_proc handle
# when c and r or i (or both) are defined.
d,a = pipe.read("time_complex.fid")
d,a = p.add(d,a,c=1.8,r=2.0,i=1.2)
pipe.write("mult4.glue",d,a,overwrite=True)

# SHUF
# NMRPipe does not ignore imaginary data when passed -rr2ri which can
# results in a double sized file, pipe_proc refuses to do this instead
# choosing to ignore imaginary data.
d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="rr2ri")
pipe.write("s8.glue",d,a,overwrite=True)

# this test will match pipe in the data portion
# but NMRPipe mis-sets the min, max parameters
d,a = pipe.read("time_real.fid")
d,a = p.shuf(d,a,mode="rr2ri")
pipe.write("s9.glue",d,a,overwrite=True)

# This function can result in nan in the array (as does pipe) and 
# therefore will not return true in np.allclose()
# Also NMRPipe does not set the scaleflag parameter
d,a = pipe.read("time_complex.fid")
d,a = p.shuf(d,a,mode="bswap")
pipe.write("s10.glue",d,a,overwrite=True)

# r2i and i2r and not implemented in prop_par 
# as pipe does not implement integer format
#d,a = pipe.read("time_complex.fid")
#d,a = p.shuf(d,a,mode="r2i")
#pipe.write("s11.glue",d,a,overwrite=True)

#d,a = pipe.read("time_complex.fid")
#d,a = p.shuf(d,a,mode="i2r")
#pipe.write("s12.glue",d,a,overwrite=True)

# DEV
# NMRPipe DEV functions goes into a infinite loop, oviously we don't
# want to reproduce this.
#d,a = pipe.read("time_complex.fid")
#d,a = p.dev(d,a)
#pipe.write("dev.glue",d,a,overwrite=True)
