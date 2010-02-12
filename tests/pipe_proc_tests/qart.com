#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn QART -a 1.0 -f 0.5  \
-ov -out qart.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn QART -a 0.8 -f 1.2  \
-ov -out qart2.dat
