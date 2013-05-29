#! /bin/csh

nmrPipe -in ./time_complex.fid                              \
| nmrPipe -fn SP -off 0.35 -end 0.98 -pow 2 -c 1.0  \
-ov -out sp1.dat

nmrPipe -in ./time_complex.fid                              \
| nmrPipe -fn SP -off 0.10 -end 0.75 -pow 1 -c 0.5  -size 200 -one \
-ov -out sp2.dat
