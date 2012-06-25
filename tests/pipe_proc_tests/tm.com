#!/bin/csh

nmrPipe -in ./time_complex.fid                              \
| nmrPipe  -fn TM -t1 100 -t2 200 -c 1.0 -inv       \
-ov -out tm1.dat

nmrPipe -in ./time_complex.fid                              \
| nmrPipe  -fn TM -t1 10 -t2 1400 -c 1.5 -start 5 -size 1490       \
-ov -out tm2.dat
