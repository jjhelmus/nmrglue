#! /usr/bin/env python
""" A collection of test which are known to fail.  These tests fail because of
design differences between NMRPipe and nmrglue, they should not be fixed. """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# ADD
# the following show the difference way NMRPipe and pipe_proc handle
# when c and r or i (or both) are defined.
d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, c=10.0, r=5.0, i=3.0)
pipe.write("add5.glue", d, a, overwrite=True)

# MULT
# the following show the difference way NMRPipe and pipe_proc handle
# when c and r or i (or both) are defined.
d, a = pipe.read("time_complex.fid")
d, a = p.add(d, a, c=1.8, r=2.0, i=1.2)
pipe.write("mult4.glue", d, a, overwrite=True)

# SHUF
# NMRPipe does not ignore imaginary data when passed -rr2ri which can
# results in a double sized file, pipe_proc refuses to do this instead
# choosing to ignore the imaginary data.
d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="rr2ri")
pipe.write("shuf8.glue", d, a, overwrite=True)

# this test will match pipe in the data portion, but NMRPipe mis-sets the
# FDMIN and FDDISPMIN parameters
d, a = pipe.read("time_real.fid")
d, a = p.shuf(d, a, mode="rr2ri")
pipe.write("shuf9.glue", d, a, overwrite=True)

# This function can result in nan in the array (as does pipe) and
# therefore will not return true in np.allclose()
# Also NMRPipe does not set the scaleflag parameter
d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="bswap")
pipe.write("shuf10.glue", d, a, overwrite=True)

# r2i and i2r and not implemented in prop_par
# as pipe does not implement integer format
# d, a = pipe.read("time_complex.fid")
# d, a = p.shuf(d, a, mode="r2i")
# pipe.write("shuf11.glue", d, a, overwrite=True)

# d, a = pipe.read("time_complex.fid")
# d, a = p.shuf(d, a, mode="i2r")
# pipe.write("shuf12.glue", d, a, overwrite=True)

# DEV
# NMRPipe DEV functions goes into a infinite loop, obviously we don't
# want to reproduce this.
# d, a = pipe.read("time_complex.fid")
# d, a = p.dev(d, a)
# pipe.write("dev1.glue", d, a, overwrite=True)
