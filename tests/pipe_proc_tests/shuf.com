#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -ri2c       \
-ov -out s1.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -c2ri       \
-ov -out s2.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -ri2rr      \
-ov -out s3.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -exlr       \
-ov -out s4.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -rolr       \
-ov -out s5.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -swap       \
-ov -out s6.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -inv        \
-ov -out s7.dat
