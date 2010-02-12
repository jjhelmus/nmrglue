#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn REV -sw  \
-ov -out rev1.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn REV -sw  \
-ov -out rev2.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn REV \
-ov -out rev3.dat
