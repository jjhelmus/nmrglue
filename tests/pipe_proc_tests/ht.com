#!/bin/csh

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn HT  \
-ov -out ht1.dat

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn HT -td \
-ov -out ht2.dat

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn HT -ps0-0 \
-ov -out ht3.dat


nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn HT -zf  \
-ov -out ht5.dat

nmrPipe -in ./1D_time_real.fid \
| nmrPipe  -fn HT -auto  \
-ov -out ht6.dat

# 2D data
nmrPipe -in ./freq_real.ft2 \
| nmrPipe  -fn HT \
-ov -out ht7.dat

nmrPipe -in ./freq_real.ft2 \
| nmrPipe  -fn HT -zf -td   \
-ov -out ht8.dat
