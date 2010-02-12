#!/bin/csh

nmrPipe -in ./time_complex.fid                              \
| nmrPipe  -fn TRI -loc 500 -lHi 0.5 -rHi 0.8  -inv \
-ov -out tri.dat

nmrPipe -in ./time_complex.fid                              \
| nmrPipe  -fn TRI -loc 750 -lHi 0.1 -rHi 0.5  \
-ov -out tri2.dat
