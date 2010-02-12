#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.rev(d,a,sw=True)
pipe.write("rev1.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.rev(d,a,sw=True)
pipe.write("rev2.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.rev(d,a,sw=False)
pipe.write("rev3.glue",d,a,overwrite=True)
