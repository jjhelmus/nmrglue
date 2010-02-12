#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

d,a = pipe.read("time_complex.fid")
d,a = p.qmix(d,a,ic=4,oc=2,cList=[0.1,0.2,0.3,0.4,
                                  0.5,0.6,0.7,0.8],time=False)
pipe.write("qmix.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.qmix(d,a,ic=2,oc=4,cList=[1.0,0.5,
                                  0.7,0.8,
                                  0.2,0.6,
                                  0.1,0.9],time=True)
pipe.write("qmix2.glue",d,a,overwrite=True)
