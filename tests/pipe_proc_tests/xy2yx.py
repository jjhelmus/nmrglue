#! /usr/bin/env python
""" Create files for xy2yx unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.xy2yx(d, a, auto=True)
pipe.write("xy2yx1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.xy2yx(d, a, hyper=True)
pipe.write("xy2yx2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.xy2yx(d, a, nohyper=True)
pipe.write("xy2yx3.glue", d, a, overwrite=True)
