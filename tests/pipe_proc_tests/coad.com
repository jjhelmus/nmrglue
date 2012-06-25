#!/bin/csh

nmrPipe -in ./time_complex.fid                 \
| nmrPipe  -fn COAD -cList 1 1 -axis X -time  \
-ov -out coad1.dat

nmrPipe -in ./time_complex.fid                 \
| nmrPipe  -fn COAD -cList 1 0 -5 8 -axis Y -time  \
-ov -out coad2.dat
