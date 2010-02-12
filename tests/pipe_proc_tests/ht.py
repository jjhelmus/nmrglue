#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a)
pipe.write("ht1.glue",d,a,overwrite=True)

d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a,td=True)
pipe.write("ht2.glue",d,a,overwrite=True)

d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a,mode="ps0-0")
pipe.write("ht3.glue",d,a,overwrite=True)

d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a,zf=True)
pipe.write("ht5.glue",d,a,overwrite=True)

d,a = pipe.read("1D_time_real.fid")
d,a = p.ht(d,a,auto=True)
pipe.write("ht6.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.ht(d,a)
pipe.write("ht7.glue",d,a,overwrite=True)

d,a = pipe.read("freq_real.ft2")
d,a = p.ht(d,a,zf=True,td=True)
pipe.write("ht8.glue",d,a,overwrite=True)
