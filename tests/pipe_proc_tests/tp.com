#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn TP  -auto\
-ov -out tp1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn TP -hyper \
-ov -out tp2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn TP  -nohyper \
-ov -out tp3.dat
