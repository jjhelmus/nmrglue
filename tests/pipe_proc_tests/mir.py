#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="left",sw=True)
pipe.write("mir1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="right",sw=True)
pipe.write("mir2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="right",invr=True,sw=True)
pipe.write("mir3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="left",invr=True,sw=True)
pipe.write("mir4.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="center",sw=True)
pipe.write("mir5.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="ps90-180",sw=True)
pipe.write("mir6.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.mir(d,a,mode="ps0-0",sw=True)
pipe.write("mir7.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="left",sw=True)
pipe.write("mir8.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="right",sw=True)
pipe.write("mir9.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="right",invr=True,sw=True)
pipe.write("mir10.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="left",invr=True,sw=True)
pipe.write("mir11.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="center",sw=True)
pipe.write("mir12.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="ps90-180",sw=True)
pipe.write("mir13.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.mir(d,a,mode="ps0-0",sw=True)
pipe.write("mir14.glue",d,a,overwrite=True)
