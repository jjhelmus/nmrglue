#! /usr/bin/env python
""" Create files for cbf unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.cs(d, a, dir="rs", pts=5.0, neg=True, sw=True)
pipe.write("cs1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.cs(d, a, dir="rs", pts=-2.0, sw=True)
pipe.write("cs2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.cs(d, a, dir="ls", pts=3.0, neg=True, sw=True)
pipe.write("cs3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.cs(d, a, dir="ls", pts=-8.0, neg=True, sw=True)
pipe.write("cs4.glue", d, a, overwrite=True)

# freq domain
d, a = pipe.read("freq_real.ft2")
d, a = p.cs(d, a, dir="ls", pts=7.0, neg=True, sw=True)
pipe.write("cs5.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.cs(d, a, dir="ls", pts=-3.0, sw=True)
pipe.write("cs6.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.cs(d, a, dir="rs", pts=9.0, neg=True, sw=True)
pipe.write("cs7.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.cs(d, a, dir="rs", pts=3.0, sw=True)
pipe.write("cs8.glue", d, a, overwrite=True)
