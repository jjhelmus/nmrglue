#! /usr/bin/env python
""" Create files for ps unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.ps(d, a, p0=11.1, p1=16.0, noup=True)
pipe.write("ps1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ps(d, a, p0=11.1, p1=16.0, noup=True)
pipe.write("ps2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ps(d, a, p0=0.1, p1=0.5, tc=1.1, exp=True, noup=True)
pipe.write("ps3.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_complex.dat")
d, a = p.ps(d, a, p0=-27.0, p1=0.0)
pipe.write("ps4.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.ps(d, a, p0=-27.0, p1=0.0, ht=True)
pipe.write("ps5.glue", d, a, overwrite=True)

d, a = pipe.read("1D_freq_real.dat")
d, a = p.ps(d, a, p0=-27.0, p1=0.0, ht=True, zf=True)
pipe.write("ps6.glue", d, a, overwrite=True)
