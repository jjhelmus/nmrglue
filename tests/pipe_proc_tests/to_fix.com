#!/bin/csh

# A collection of tests which fail for one reason or another

# FSH
nmrPipe -in ./1D_freq_real.dat  \
| nmrPipe  -fn FSH  -ls 1 -sw      \
-ov -out fsh1.dat

nmrPipe -in ./1D_freq_real.dat  \
| nmrPipe  -fn FSH  -rs 8 -sw      \
-ov -out fsh2.dat

nmrPipe -in ./1D_freq_real.dat  \
| nmrPipe  -fn FSH  -rs 6.7 -sw    \
-ov -out fsh3.dat

nmrPipe -in ./1D_freq_complex.dat  \
| nmrPipe  -fn FSH  -ls 9.5  -sw   \
-ov -out fsh4.dat

# FSH
nmrPipe -in ./1D_freq_complex.dat \
| nmrPipe  -fn RFT                \
-ov -out rft9.dat

nmrPipe -in ./1D_freq_real.dat \
| nmrPipe  -fn RFT             \
-ov -out rft10.dat

nmrPipe -in ./freq_real.ft2 \
| nmrPipe  -fn RFT          \
-ov -out rft11.dat


# HT
nmrPipe -in ./1D_freq_real.dat \
| nmrPipe  -fn HT -ps90-180 \
-ov -out ht4.dat

# Integration tests
nmrPipe -in time_real.fid                   \
|   nmrPipe -fn GMB -gb 0.1 -lb -8 -c 0.5             \
|   nmrPipe -fn ZF -auto                \
|   nmrPipe -fn FT -alt                 \
|   nmrPipe -fn PS -p0 0 -p1 0       \
|   nmrPipe -fn TP -hyper                    \
|   nmrPipe -fn SP -off 0.5 -pow 2.0 -c 0.5  \
|   nmrPipe -fn ZF -auto            \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn PS -p0 0 -p1 0 -di     \
-ov -out 2d_mixed_processing1.dat

nmrPipe -in time_real.fid                   \
|   nmrPipe -fn EM -lb 8           \
|   nmrPipe -fn ZF -auto                \
|   nmrPipe -fn FT -auto -di                    \
|   nmrPipe -fn TP                       \
|   nmrPipe -fn SP -off 0.5 -pow 1.0 -c 0.5  \
|   nmrPipe -fn ZF -auto            \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn MC    \
-ov -out 2d_mixed_processing2.dat