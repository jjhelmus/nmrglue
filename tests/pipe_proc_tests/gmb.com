#!/bin/csh

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn GMB -lb 2.0 -gb 0.5 -c 1.0 -inv\
-ov -out gmb1.dat

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn GMB -lb 10.0 -gb 0.2 -c 0.5 -start 20 -size 800\
-ov -out gmb2.dat
