#!/bin/csh
nmrPipe -in test.fid \
| nmrPipe  -fn SP -off 0.35 -end 0.98 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                                \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 151.0 -p1 0.0 -di               \
     -verb -ov -out test.ft
