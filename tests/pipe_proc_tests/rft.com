#!/bin/csh

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn RFT             \
-ov -out rft1.dat

nmrPipe -in ./1D_time.fid      \
| nmrPipe  -fn RFT             \
-ov -out rft2.dat

nmrPipe -in ./time_real.fid    \
| nmrPipe  -fn RFT             \
-ov -out rft3.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn RFT             \
-ov -out rft4.dat

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn RFT -inv        \
-ov -out rft5.dat

nmrPipe -in ./1D_time.fid      \
| nmrPipe  -fn RFT -inv        \
-ov -out rft6.dat

nmrPipe -in ./time_real.fid    \
| nmrPipe  -fn RFT -inv        \
-ov -out rft7.dat

nmrPipe -in ./time_complex.fid \
| nmrPipe  -fn RFT -inv        \
-ov -out rft8.dat

nmrPipe -in ./1D_freq_complex.dat \
| nmrPipe  -fn RFT -inv           \
-ov -out rft12.dat

nmrPipe -in ./1D_freq_real.dat \
| nmrPipe  -fn RFT -inv        \
-ov -out rft13.dat

nmrPipe -in ./freq_real.ft2 \
| nmrPipe  -fn RFT -inv     \
-ov -out rft14.dat


