#! /usr/bin/env python
""" Create files for shuf unit test """

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="ri2c")
pipe.write("shuf1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="c2ri")
pipe.write("shuf2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="ri2rr")
pipe.write("shuf3.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="exlr")
pipe.write("shuf4.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="rolr")
pipe.write("shuf5.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="swap")
pipe.write("shuf6.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.shuf(d, a, mode="inv")
pipe.write("shuf7.glue", d, a, overwrite=True)
