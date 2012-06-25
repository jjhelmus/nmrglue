#!/bin/csh

nmrPipe -in ./time_complex.fid               \
| nmrPipe  -fn INTEG  \
-ov -out integ1.dat
