#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZF -zf 2  \
-ov -out zf.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZF -pad 200  \
-ov -out zf2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZF -size 8000  \
-ov -out zf3.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZF -size 4096 -mid  \
-ov -out zf4.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn ZF -size 4096 -inter  \
-ov -out zf5.dat
