#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX  -auto\
-ov -out xy.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX -hyper \
-ov -out xy2.dat


nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX  -nohyper \
-ov -out xy3.dat
