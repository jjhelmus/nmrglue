#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn YTP  -auto\
-ov -out ytp1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn YTP -hyper \
-ov -out ytp2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn YTP  -nohyper \
-ov -out ytp3.dat
