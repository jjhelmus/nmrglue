#!/bin/csh

nmrPipe -in ./nmrpipe_1d/test.fid \
| nmrPipe  -fn SP -off 0.35 -end 0.98 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                                \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 -17.7 -p1 -36.0 -di              \
   -verb -ov -out ./nmrpipe_1d/test.ft
