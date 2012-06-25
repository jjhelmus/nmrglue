#!/bin/csh

# A collection of tests which fail due to intended differences between NMRPipe
# and nmrglue.  These should NOT be fixed.

# ADD
nmrPipe -in ./time_complex.fid               \
| nmrPipe  -fn ADD -c 10.0  -i 5.0 -r 3.0    \
-ov -out add5.dat

# MULT
nmrPipe -in ./time_complex.fid             \
| nmrPipe  -fn MULT -c 1.8 -r 2.0 -i 1.2   \
-ov -out mult4.dat

# SHUF
nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -rr2ri      \
-ov -out shuf8.dat

nmrPipe -in ./time_real.fid \
| nmrPipe -fn SHUF -rr2ri   \
-ov -out shuf9.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn SHUF -bswap      \
-ov -out shuf10.dat

#nmrPipe -in ./time_complex.fid  \
#| nmrPipe  -fn SHUF -r2i        \
#-ov -out shuf11.dat

#nmrPipe -in ./time_complex.fid  \
#| nmrPipe  -fn SHUF -i2r        \
#-ov -out shuf12.dat

# DEV
#nmrPipe -in ./time_complex.fid \
#| nmrPipe  -fn DEV  \
#-ov -out dev.dat

