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


nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT -auto \
-ov -out time-freq.c.ft1

nmrPipe -in ./time-freq.c.ft1                \
| nmrPipe  -fn TP -hyper\
-ov -out tp4.dat

nmrPipe -in ./time-freq.c.ft1                \
| nmrPipe  -fn TP -auto\
-ov -out tp5.dat

rm -f time-freq.c.ft1

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT -auto -di \
-ov -out time-freq.r.ft1

nmrPipe -in ./time-freq.r.ft1                \
| nmrPipe -fn TP \
-ov -out tp6.dat

nmrPipe -in ./time-freq.r.ft1                \
| nmrPipe  -fn TP -auto \
-ov -out tp7.dat

rm -f time-freq.r.ft1
