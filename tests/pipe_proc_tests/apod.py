#! /usr/bin/env python

import nmrglue.fileio.pipe as pipe
import nmrglue.process.pipe_proc as p

# Apod Tests
d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="SP",q1=0.35,q2=0.98,q3=2.0)
pipe.write("apod1.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="EM",q1=10.0)
pipe.write("apod2.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="GM",q1=2.35,q2=1.25,q3=1.2)
pipe.write("apod3.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="GMB",q1=1.25,q2=3.0)
pipe.write("apod4.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="TM",q1=100,q2=200)
pipe.write("apod5.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="TRI",q1=500,q2=0.50,q3=0.8)
pipe.write("apod6.glue",d,a,overwrite=True)

d,a = pipe.read("time_complex.fid")
d,a = p.apod(d,a,qName="JMOD",q1=5.0,q2=2.5,q3=1.2)
pipe.write("apod7.glue",d,a,overwrite=True)
