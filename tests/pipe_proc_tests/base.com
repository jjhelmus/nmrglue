#!/bin/csh

nmrPipe -in ./1D_freq_real.dat   \
| nmrPipe  -fn BASE -nl 100 200 300  \
-ov -out base1.dat

nmrPipe -in ./1D_freq_complex.dat   \
| nmrPipe  -fn BASE -nl 100 200  \
-ov -out base2.dat

nmrPipe -in ./freq_real.ft2      \
| nmrPipe  -fn BASE -nl 100 200  \
-ov -out base3.dat

nmrPipe -in ./1D_freq_real.dat      \
| nmrPipe  -fn BASE -nl 200 300 800 \
-ov -out base4.dat

nmrPipe -in ./1D_freq_real.dat             \
| nmrPipe  -fn BASE -nl 200 300 650 -first \
-ov -out base5.dat

nmrPipe -in ./1D_freq_real.dat        \
| nmrPipe  -fn BASE -nl 200 560 -last \
-ov -out base6.dat

nmrPipe -in ./1D_freq_real.dat        \
| nmrPipe  -fn BASE -nl 150 250 -nw 3 \
-ov -out base7.dat
