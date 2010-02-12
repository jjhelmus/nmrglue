#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn LS -ls 2.0 -sw \
-ov -out ls1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn LS -ls -3.0 -sw \
-ov -out ls2.dat

# freq domain 
nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn LS -ls 2.0 -sw \
-ov -out ls3.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn LS -ls 17.0 -sw \
-ov -out ls4.dat

nmrPipe -in ./freq_real.ft2                \
| nmrPipe  -fn LS -ls -5.0 -sw \
-ov -out ls5.dat
