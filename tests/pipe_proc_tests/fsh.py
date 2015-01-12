#! /usr/bin/env python
""" Create files for fsh unit test """

# import nmrglue.fileio.pipe as pipe
# import nmrglue.process.pipe_proc as p

# the first three tests will fail because MIN/MAX values are off a little
# d, a = pipe.read("1D_freq_real.dat")
# d, a = p.fsh(d, a, dir="ls", pts=1)
# pipe.write("fsh1.glue", d, a, overwrite=True)

# d, a = pipe.read("1D_freq_real.dat")
# d, a = p.fsh(d, a, dir="rs", pts=8)
# pipe.write("fsh2.glue", d, a, overwrite=True)

# d, a = pipe.read("1D_freq_real.dat")
# d, a = p.fsh(d, a, dir="rs", pts=6.7)
# pipe.write("fsh3.glue", d, a, overwrite=True)

# this fails because NMRPipe performs a Hilbert transform?
# d, a = pipe.read("1D_freq_complex.dat")
# d, a = p.fsh(d, a, dir="ls", pts=9.5)
# pipe.write("fsh4.glue", d, a, overwrite=True)
