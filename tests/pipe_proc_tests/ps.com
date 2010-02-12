#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn PS -p0 11.1 -p1 16.0           \
-ov -out ps1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn PS -p0 11.1 -p1 16.0  -noup    \
-ov -out ps2.dat

nmrPipe -in ./time_complex.fid                  \
| nmrPipe  -fn PS -exp -p0 0.1 -p1 0.5 -tc 1.1  \
-ov -out ps3.dat

nmrPipe -in ./1D_freq_complex.dat                \
| nmrPipe  -fn PS -p0 -27.0 -p1 0.0      \
-ov -out ps4.dat

nmrPipe -in ./1D_freq_real.dat                \
| nmrPipe  -fn PS -ht -p0 -27.0 -p1 0.0      \
-ov -out ps5.dat

nmrPipe -in ./1D_freq_real.dat                \
| nmrPipe  -fn PS -ht -zf -p0 -27.0 -p1 0.0      \
-ov -out ps6.dat
