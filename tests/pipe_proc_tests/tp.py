#! /usr/bin/env python
""" Create files for tp unit test """

from subprocess import check_output

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, auto=True)
pipe.write("tp1.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, hyper=True)
pipe.write("tp2.glue", d, a, overwrite=True)

d, a = pipe.read("time_complex.fid")
d, a = p.tp(d, a, nohyper=True)
pipe.write("tp3.glue", d, a, overwrite=True)


pipe_command = """\
nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT  -auto \
-ov -out time-freq.c.ft1"""
check_output(pipe_command, shell=True)

d, a = pipe.read("time-freq.c.ft1")
d, a = p.tp(d, a, hyper=True)
pipe.write("tp4.glue", d, a, overwrite=True)

d, a = pipe.read("time-freq.c.ft1")
d, a = p.tp(d, a, auto=True)
pipe.write("tp5.glue", d, a, overwrite=True)

os.remove("time-freq.c.ft1")

pipe_command = """\
nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT  -auto -di \
-ov -out time-freq.r.ft1"""
check_output(pipe_command, shell=True)

d, a = pipe.read("time-freq.r.ft1")
d, a = p.tp(d, a)
pipe.write("tp6.glue", d, a, overwrite=True)

d, a = pipe.read("time-freq.r.ft1")
d, a = p.tp(d, a, auto=True)
pipe.write("tp7.glue", d, a, overwrite=True)

os.remove("time-freq.r.ft1")
