#!/bin/csh

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn MIR -left  -sw      \
-ov -out mir1.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe  -fn MIR -right -sw       \
-ov -out mir2.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe  -fn MIR -right -invr -sw       \
-ov -out mir3.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe  -fn MIR -left -invr -sw       \
-ov -out mir4.dat

nmrPipe -in ./time_complex.fid   \
| nmrPipe  -fn MIR -center  -sw     \
-ov -out mir5.dat

nmrPipe -in ./time_complex.fid  \
| nmrPipe  -fn MIR -ps90-180  -sw  \
-ov -out mir6.dat

nmrPipe -in ./time_complex.fid        \
| nmrPipe  -fn MIR -ps0-0    -sw         \
-ov -out mir7.dat

nmrPipe -in ./1D_freq_real.dat  \
| nmrPipe  -fn MIR -left  -sw      \
-ov -out mir8.dat

nmrPipe -in ./1D_freq_real.dat   \
| nmrPipe  -fn MIR -right -sw       \
-ov -out mir9.dat

nmrPipe -in ./1D_freq_real.dat   \
| nmrPipe  -fn MIR -right -invr -sw       \
-ov -out mir10.dat

nmrPipe -in ./1D_freq_real.dat   \
| nmrPipe  -fn MIR -left -invr -sw       \
-ov -out mir11.dat

nmrPipe -in ./1D_freq_real.dat   \
| nmrPipe  -fn MIR -center  -sw     \
-ov -out mir12.dat

nmrPipe -in ./1D_freq_real.dat  \
| nmrPipe  -fn MIR -ps90-180  -sw  \
-ov -out mir13.dat

nmrPipe -in ./1D_freq_real.dat        \
| nmrPipe  -fn MIR -ps0-0    -sw         \
-ov -out mir14.dat
