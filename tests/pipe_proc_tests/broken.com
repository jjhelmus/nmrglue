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

# ADD
nmrPipe -in ./time_complex.fid               \
| nmrPipe  -fn ADD -c 10.0  -i 5.0 -r 3.0    \
-ov -out add5.dat

nmrPipe -in ./time_complex.fid             \
| nmrPipe  -fn MULT -c 1.8 -r 2.0 -i 1.2   \
-ov -out mult4.dat


# SHUF
nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -rr2ri      \
-ov -out s8.dat

nmrPipe -in ./time_real.fid \
| nmrPipe -fn SHUF -rr2ri   \
-ov -out s9.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -bswap      \
-ov -out s10.dat

#nmrPipe -in ./time_complex.fid  \
#| nmrPipe  -fn SHUF -r2i        \
#-ov -out s11.dat
#nmrPipe -in ./time_complex.fid  \
#| nmrPipe  -fn SHUF -i2r        \
#-ov -out s12.dat

# DEV
#nmrPipe -in ./time_complex.fid \
#| nmrPipe  -fn DEV  \
#-ov -out dev.dat

