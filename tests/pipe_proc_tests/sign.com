#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -ri\
-ov -out s1.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -r\
-ov -out s2.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -i\
-ov -out s3.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -left\
-ov -out s4.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -right\
-ov -out s5.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -alt\
-ov -out s6.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -abs\
-ov -out s7.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -sign\
-ov -out s8.dat
