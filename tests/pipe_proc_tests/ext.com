#!/bin/csh

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -left -sw \
-ov -out ext1.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -right  -sw\
-ov -out ext2.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -mid -sw \
-ov -out ext3.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -x1 1 -xn 100 -sw \
-ov -out ext4.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -xn 200 -y1 50 -yn 75 -sw \
-ov -out ext5.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -x1 5 -xn 200 -pow2  -sw\
-ov -out ext6.dat

nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn EXT -x1 5 -xn 200 -round 10 -sw \
-ov -out ext7.dat


nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT -auto \
-ov -out time-freq.c.ft1

nmrPipe -in ./time-freq.c.ft1                \
| nmrPipe  -fn EXT -left -sw\
-ov -out ext8.dat

nmrPipe -in ./time-freq.c.ft1                \
| nmrPipe  -fn EXT -right -sw\
-ov -out ext9.dat

rm -f time-freq.c.ft1


nmrPipe -in ./time_complex.fid                \
| nmrPipe  -fn FT -auto -di \
-ov -out time-freq.r.ft1

nmrPipe -in ./time-freq.r.ft1                \
| nmrPipe  -fn EXT -left -sw\
-ov -out ext10.dat

nmrPipe -in ./time-freq.r.ft1                \
| nmrPipe  -fn EXT -right -sw\
-ov -out ext11.dat

rm -f time-freq.r.ft1
