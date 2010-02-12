#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn SMO  \
-ov -out smo.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn SMO -n 5\
-ov -out smo2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn SMO -center  \
-ov -out smo3.dat
