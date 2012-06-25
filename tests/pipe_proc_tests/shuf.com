#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -ri2c       \
-ov -out shuf1.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -c2ri       \
-ov -out shuf2.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -ri2rr      \
-ov -out shuf3.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -exlr       \
-ov -out shuf4.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -rolr       \
-ov -out shuf5.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -swap       \
-ov -out shuf6.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -inv        \
-ov -out shuf7.dat
