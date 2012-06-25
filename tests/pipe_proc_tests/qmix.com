#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn QMIX -ic 4  -oc 2 -cList 0.1 0.2 0.3 0.4\
                                             0.5 0.6 0.7 0.8\
-ov -out qmix1.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn QMIX -ic 2 -oc 4 -time -cList 1.0 0.5\
                                             0.7 0.8\
                                             0.2 0.6\
                                             0.1 0.9\
-ov -out qmix2.dat
