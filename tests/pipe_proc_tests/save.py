#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.save(d,a,"save.glue",overwrite=True)
pipe.write("save2.glue",d,a,overwrite=True)
