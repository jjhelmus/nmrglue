#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX  -auto\
-ov -out xy2yx1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX -hyper \
-ov -out xy2yx2.dat


nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn XY2YX  -nohyper \
-ov -out xy2yx3.dat
