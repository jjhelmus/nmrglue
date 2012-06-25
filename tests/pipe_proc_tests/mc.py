#! /usr/bin/env python
""" tests for MC function """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.mc(d, a, mode="pow")
pipe.write("mc1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.mc(d, a, mode="mod")
pipe.write("mc2.glue", d, a, overwrite=True)
