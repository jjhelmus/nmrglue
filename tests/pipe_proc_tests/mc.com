#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn MC -pow \
-ov -out mc1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn MC -mod \
-ov -out mc2.dat
