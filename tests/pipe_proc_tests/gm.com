#!/bin/csh

nmrPipe -in ./time_complex.fid          \
| nmrPipe  -fn GM -g1 5.0 -g2 2.0 -g3 0.0 -c 1.0 \
-ov -out gm1.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe -fn GM -g1 2.0 -g2 5.0 -g3 0.5 -c 1.5 -start 100 -size 654 -one\
 -ov -out gm2.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe -fn GM -g1 2.0 -g2 5.0 -g3 0.5 -c 1.5 -start 100 \
 -ov -out gm3.dat
