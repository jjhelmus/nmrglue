#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn DX  \
-ov -out dx.dat
