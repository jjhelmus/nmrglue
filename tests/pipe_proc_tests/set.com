#!/bin/csh

nmrPipe -in ./time_complex.fid    \
| nmrPipe  -fn SET  -r 1.0 -i 2.0 \
-ov -out set1.dat

nmrPipe -in ./time_complex.fid    \
| nmrPipe  -fn SET  -c 3.0        \
-ov -out set2.dat

nmrPipe -in ./time_complex.fid             \
| nmrPipe  -fn SET  -r 8.0 -i -2.0 -c 10.0 \
-ov -out set3.dat

nmrPipe -in ./time_complex.fid                   \
| nmrPipe  -fn SET  -r 1.0 -i 2.0 -x1 50 -xn 700 \
-ov -out set4.dat
