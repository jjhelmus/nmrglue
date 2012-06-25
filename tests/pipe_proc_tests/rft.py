#! /usr/bin/env python
""" Creat files for rft unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("1D_time_real.fid")
d, a = p.rft(d, a)
pipe.write("rft1.glue", d, a, overwrite=True)

d, a = pipe.read("1D_time.fid")
d, a = p.rft(d, a)
pipe.write("rft2.glue", d, a, overwrite=True)

d, a = pipe.read("time_real.fid")
d, a = p.rft(d, a)
pipe.write("rft3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.rft(d, a)
pipe.write("rft4.glue", d, a, overwrite=True)

d, a = pipe.read("1D_time_real.fid")
d, a = p.rft(d, a, inv=True)
pipe.write("rft5.glue", d, a, overwrite=True)

d, a = pipe.read("1D_time.fid")
d, a = p.rft(d, a, inv=True)
pipe.write("rft6.glue", d, a, overwrite=True)

d, a = pipe.read("time_real.fid")
d, a = p.rft(d, a, inv=True)
pipe.write("rft7.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.rft(d, a, inv=True)
pipe.write("rft8.glue", d, a, overwrite=True)

# Frequency domain
d, a = pipe.read("1D_freq_complex.dat")
d, a = p.rft(d, a, inv=True)
pipe.write("rft12.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.rft(d, a, inv=True)
pipe.write("rft13.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.rft(d, a, inv=True)
pipe.write("rft14.glue", d, a, overwrite=True)
