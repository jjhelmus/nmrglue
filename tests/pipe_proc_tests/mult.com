#!/bin/csh

nmrPipe -in ./time_complex.fid    \
| nmrPipe  -fn MULT -c 1.5 -inv   \
-ov -out mult1.dat

nmrPipe -in ./time_complex.fid      \
| nmrPipe  -fn MULT -r 2.0 -i 1.2   \
-ov -out mult2.dat

nmrPipe -in ./time_complex.fid      \
| nmrPipe  -fn MULT -r 1.5 -i 0.9 -x1 90 -xn 600   \
-ov -out mult3.dat

