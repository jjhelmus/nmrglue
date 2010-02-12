#!/bin/csh

nmrPipe -in ./1D_time.fid  \
| nmrPipe  -fn HA          \
-ov -out ha1.dat

nmrPipe -in ./1D_time.fid  \
| nmrPipe  -fn HA -inv     \
-ov -out ha2.dat
