#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("1D_freq_real.dat")
d,a = p.base(d,a,nl=[100,200,300])
pipe.write("base1.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_complex.dat")
d,a = p.base(d,a,nl=[100,200])
pipe.write("base2.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.base(d,a,nl=[100,200])
pipe.write("base3.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.base(d,a,nl=[200,300,800])
pipe.write("base4.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.base(d,a,nl=[200,300,650],first=True)
pipe.write("base5.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.base(d,a,nl=[200,560],last=True)
pipe.write("base6.glue",d,a,overwrite=True)

d,a = pipe.read("1D_freq_real.dat")
d,a = p.base(d,a,nl=[150,250],nw=3)
pipe.write("base7.glue",d,a,overwrite=True)
