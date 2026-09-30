#! /usr/bin/env python
"""Create the nmrglue file for the HT ps90-180 comparison."""

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p


d, a = pipe.read("1D_time_real.fid")
d, a = p.ht(d, a, mode="ps90-180")
pipe.write("ht41.glue", d, a, overwrite=True)
