#!/bin/csh

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn JMOD -cos -j 10.0 -lb 5.0 -c 1.0 \
-ov -out jmod1.dat

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn JMOD -sin -j 18.0 -lb 1.4 -c 0.5 -start 100 -size 800 -one \
-ov -out jmod2.dat
