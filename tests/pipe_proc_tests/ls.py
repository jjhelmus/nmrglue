#! /usr/bin/env python
""" Create files for ls unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.ls(d, a, ls=2.0, sw=True)
pipe.write("ls1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.ls(d, a, ls=-3.0, sw=True)
pipe.write("ls2.glue", d, a, overwrite=True)

# freq domain
d, a = pipe.read("freq_real.ft2")
d, a = p.ls(d, a, ls=2.0, sw=True)
pipe.write("ls3.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.ls(d, a, ls=17.0, sw=True)
pipe.write("ls4.glue", d, a, overwrite=True)

d, a = pipe.read("freq_real.ft2")
d, a = p.ls(d, a, ls=-5.0, sw=True)
pipe.write("ls5.glue", d, a, overwrite=True)
