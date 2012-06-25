#!/bin/csh

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn EM -lb 40 -c 0.5 -start 100 -size 750 -one\
-ov -out em1.dat

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn EM -lb 20 -c 1.5 -size 900 \
-ov -out em2.dat
