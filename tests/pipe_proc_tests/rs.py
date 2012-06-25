#! /usr/bin/env python
""" Create files for rs unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.rs(d, a, rs=2.0, sw=True)
pipe.write("rs1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.rs(d, a, rs=-3.0, sw=True)
pipe.write("rs2.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.rs(d, a, rs=2.0, sw=True)
pipe.write("rs3.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.rs(d, a, rs=17.0, sw=True)
pipe.write("rs4.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.rs(d, a, rs=-5.0, sw=True)
pipe.write("rs5.glue", d, a, overwrite=True)
