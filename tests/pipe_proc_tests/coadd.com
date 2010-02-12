#!/bin/csh

nmrPipe -in ./time_complex.fid                 \
| nmrPipe  -fn COADD -cList 1 1 -axis X -time  \
-ov -out ca1.dat

nmrPipe -in ./time_complex.fid                 \
| nmrPipe  -fn COADD -cList 1 0 -5 8 -axis Y -time  \
-ov -out ca2.dat
