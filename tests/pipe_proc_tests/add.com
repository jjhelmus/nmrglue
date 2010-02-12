#!/bin/csh

nmrPipe -in ./time_complex.fid     \
| nmrPipe  -fn ADD -r 2.0 -i -1.0  \
-ov -out add1.dat

nmrPipe -in ./time_complex.fid     \
| nmrPipe  -fn ADD -c 10.0         \
-ov -out add2.dat

nmrPipe -in ./time_complex.fid            \
| nmrPipe  -fn ADD -c 10.0 -x1 10 -xn 400 \
-ov -out add3.dat

nmrPipe -in ./time_complex.fid        \
| nmrPipe  -fn ADD -ri -x1 50 -xn 300 \
-ov -out add4.dat

