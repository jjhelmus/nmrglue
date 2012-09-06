#!/bin/csh

nmrPipe -in ./nmrpipe_2d_tppi/test.fid \
| nmrPipe  -fn SP -off 0.35 -end 0.98 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                                \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 151.0 -p1 0.0 -di                \
| nmrPipe  -fn TP                                      \
| nmrPipe  -fn SP -off 0.35 -end 0.98 -pow 2 -c 0.5    \
| nmrPipe  -fn ZF -auto                                \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 0 -p1 0 -di                      \
| nmrPipe  -fn REV -sw                                 \
| nmrPipe  -fn TP                                      \
     -verb -ov -out ./nmrpipe_2d_tppi/test.ft2
