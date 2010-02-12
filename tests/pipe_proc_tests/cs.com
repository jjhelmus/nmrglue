#!/bin/csh

nmrPipe -in ./time_complex.fid      \
| nmrPipe  -fn CS -rs 5.0 -neg  -sw \
-ov -out cs1.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe  -fn CS -rs -2.0  -sw  \
-ov -out cs2.dat

nmrPipe -in ./time_complex.fid     \
| nmrPipe  -fn CS -ls 3.0 -neg -sw \
-ov -out cs3.dat

nmrPipe -in ./time_complex.fid         \
| nmrPipe  -fn CS -ls -8.0 -neg  -sw   \
-ov -out cs4.dat

# freq domain
nmrPipe -in ./freq_real.ft2           \
| nmrPipe  -fn CS -ls 7.0 -neg  -sw   \
-ov -out cs5.dat

nmrPipe -in ./freq_real.ft2       \
| nmrPipe  -fn CS -ls -3.0  -sw   \
-ov -out cs6.dat

nmrPipe -in ./freq_real.ft2           \
| nmrPipe  -fn CS -rs 9.0 -neg  -sw   \
-ov -out cs7.dat

nmrPipe -in ./freq_real.ft2           \
| nmrPipe  -fn CS -rs 3.0 -sw   \
-ov -out cs8.dat
