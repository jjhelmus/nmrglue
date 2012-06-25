#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -ri\
-ov -out sign1.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -r\
-ov -out sign2.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -i\
-ov -out sign3.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -left\
-ov -out sign4.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -right\
-ov -out sign5.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -alt\
-ov -out sign6.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -abs\
-ov -out sign7.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SIGN -sign\
-ov -out sign8.dat
