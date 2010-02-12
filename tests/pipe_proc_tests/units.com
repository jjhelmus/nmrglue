#!/bin/csh

# points
nmrPipe -in ./1D_time.fid                \
| nmrPipe  -fn EXT -x1 30 -xn 340.6 -sw           \
-ov -out u1.dat

# percent
nmrPipe -in ./1D_time.fid                \
| nmrPipe  -fn EXT -x1 11% -xn 60.2% -sw           \
-ov -out u2.dat

nmrPipe -in ./1D_time.fid                \
| nmrPipe  -fn EXT -x1 20.5% -xn 78.9% -sw           \
-ov -out u3.dat

# Hz
nmrPipe -in ./1D_freq_real.dat                  \
| nmrPipe  -fn EXT -x1 24000Hz -xn 12000Hz -sw           \
-ov -out u4.dat

nmrPipe -in ./1D_freq_real.dat                  \
| nmrPipe  -fn EXT -x1 15693.7Hz -xn -1634.1Hz -sw           \
-ov -out u5.dat

# PPM
nmrPipe -in ./1D_freq_real.dat                  \
| nmrPipe  -fn EXT -x1 120.2ppm -xn -12.1ppm -sw           \
-ov -out u6.dat

nmrPipe -in ./1D_freq_real.dat                  \
| nmrPipe  -fn EXT -x1 67.23ppm -xn 11.2ppm -sw           \
-ov -out u7.dat

# 2D points
nmrPipe -in ./freq_real.ft2                  \
| nmrPipe  -fn EXT -x1 193 -xn 1812 -sw      \
-ov -out u8.dat

nmrPipe -in ./freq_real.ft2                             \
| nmrPipe  -fn EXT -x1 679 -xn 1456 -y1 45 -yn 80 -sw      \
-ov -out u9.dat

nmrPipe -in ./time_complex.fid                  \
| nmrPipe  -fn EXT -x1 189 -xn 798 -sw          \
-ov -out u10.dat

nmrPipe -in ./time_complex.fid                  \
| nmrPipe  -fn EXT -x1 87 -xn 991 -y1 182 -yn 310 -sw \
-ov -out u11.dat

# 2D percent
nmrPipe -in ./freq_real.ft2                  \
| nmrPipe  -fn EXT -x1 22.4% -xn 67.8% -sw      \
-ov -out u12.dat

nmrPipe -in ./freq_real.ft2                             \
| nmrPipe  -fn EXT -x1 12.6% -xn 20% -y1 89% -yn 99% -sw      \
-ov -out u13.dat

# 2D Hz
nmrPipe -in ./freq_real.ft2                  \
| nmrPipe  -fn EXT -x1 13203hz -xn -1560hz -sw      \
-ov -out u14.dat

nmrPipe -in ./freq_real.ft2                             \
| nmrPipe  -fn EXT -x1 10239hz -xn -19341hz -y1 1333hz -yn -1234hz -sw      \
-ov -out u15.dat

# 2D ppm
nmrPipe -in ./freq_real.ft2                  \
| nmrPipe  -fn EXT -x1 145.2ppm -xn -11.2ppm -sw      \
-ov -out u16.dat

nmrPipe -in ./freq_real.ft2                             \
| nmrPipe  -fn EXT -x1 12.3ppm -xn -101.2ppm -y1 23.0ppm -yn -14.2ppm -sw      \
-ov -out u17.dat


