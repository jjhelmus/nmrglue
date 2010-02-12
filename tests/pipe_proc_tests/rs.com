#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn RS -rs 2.0 -sw \
-ov -out rs1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn RS -rs -3.0 -sw \
-ov -out rs2.dat

# freq domain 
nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn RS -rs 2.0 -sw \
-ov -out rs3.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn RS -rs 17.0 -sw \
-ov -out rs4.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn RS -rs -5.0 -sw \
-ov -out rs5.dat
