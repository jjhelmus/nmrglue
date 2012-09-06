#!/bin/csh

nmrPipe -in ./nmrpipe_2d/test.fid                              \
| nmrPipe  -fn SP -off 0.35 -end 0.98 -pow 1 -c 1.0 \
| nmrPipe  -fn ZF -auto                             \
| nmrPipe  -fn FT -auto                             \
| nmrPipe  -fn PS -p0 -29.0 -p1  0.0 -di            \
| nmrPipe  -fn TP                                   \
| nmrPipe  -fn SP -off 0.35 -end 0.9 -pow 1 -c 0.5 \
| nmrPipe  -fn ZF -size 2048                        \
| nmrPipe  -fn FT -auto                             \
| nmrPipe  -fn PS -p0 0 -p1 0   -di           	    \
| nmrPipe  -fn TP                                   \
    -verb -ov -out ./nmrpipe_2d/test.ft2
