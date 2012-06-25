#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZD -wide 1 -x0 20 -slope 2 -func 0 \
-ov -out zd1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZD -wide 5 -x0 10 -slope 5 -func 1 \
-ov -out zd2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZD -wide 3 -x0 100 -slope 1 -func 2 \
-ov -out zd3.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZD -wide 8 -x0 15 -slope 3 -func 3 -g 20 \
-ov -out zd4.dat
