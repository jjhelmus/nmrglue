#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn CBF  \
-ov -out cbf1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn CBF -last 30 \
-ov -out cbf2.dat

# reg and slice not implemented
nmrPipe -in ./time_complex.fid            \
| nmrPipe  -fn CBF  -reg 2 100  \
-ov -out cbf3.dat

nmrPipe -in ./time_complex.fid      \
| nmrPipe  -fn CBF  -slice 10 20    \
-ov -out cbf4.dat
