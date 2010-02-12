#!/bin/csh

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT         \
-ov -out ft1.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -real        \
-ov -out ft2.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -inv        \
-ov -out ft3.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -alt        \
-ov -out ft4.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -neg        \
-ov -out ft5.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -null        \
-ov -out ft6.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -bruk        \
-ov -out ft7.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn FT -auto        \
-ov -out ft8.dat
