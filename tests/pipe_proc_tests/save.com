#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SAVE -name save.dat \
-ov -out save2.dat
