#!/bin/csh


# Script to process S3E filtered data
# By: Jonathan Helmus (helmus.3@osu.edu)

# data must be converted to nmrpipe data with twice the number of points
# in the indirect dimension using var2pipe and the A/B S3E block the intermost
# array

# var2pipe sign alternates the Y vectors for States so we have:
# A_0 -B_0 A_1 -B_1
# For the sum spectra we want to return:
# A_0+B_0 -(A_1+B_1) 
# and for the difference spectra:
# A_0-B_0 -(A_1-B_1)
# this is accomplished by COADD Y vectors and then sign alternating 
# using the QMIX function


# process the sum spectrum
nmrPipe -in test.fid \
| nmrPipe  -fn COADD -axis Y -cList 1 -1 -time         \
| nmrPipe  -fn QMIX -ic 2 -oc 2 -cList 1  0            \
                                       0 -1            \
| nmrPipe  -fn SP -off 0.45 -end 0.95 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -size 8192                           \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 0.0 -p1 0.0                      \
| nmrPipe  -fn FSH -ls 27.5Hz                       \
| nmrPipe  -fn PS -p0 101.0 -p1 0.0 -di                     \
| nmrPipe  -fn TP                                      \
| nmrPipe  -fn SP -off 0.45 -end 0.95 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 2048                           \
| nmrPipe  -fn FT -neg                                \
| nmrPipe  -fn PS -p0 0 -p1 0 -di                      \
#| nmrPipe  -fn REV 				                   \
| nmrPipe  -fn TP                                      \
     -verb -ov -out sum.ft2

# process the difference spectrum
nmrPipe -in test.fid \
| nmrPipe  -fn COADD -axis Y -cList 1 1 -time          \
| nmrPipe  -fn QMIX -ic 2 -oc 2 -cList 1  0            \
                                       0 -1            \
| nmrPipe  -fn SP -off 0.45 -end 0.95 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -size 8192                           \
| nmrPipe  -fn FT -auto                                \
| nmrPipe  -fn PS -p0 -90.0 -p1 0.0                    \
| nmrPipe  -fn FSH -rs 27.5Hz                       \
| nmrPipe  -fn PS -p0 101.0 -p1 0.0 -di                     \
| nmrPipe  -fn TP                                      \
| nmrPipe  -fn SP -off 0.45 -end 0.95 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 2048                           \
| nmrPipe  -fn FT -neg                                \
| nmrPipe  -fn PS -p0 0 -p1 0 -di                      \
#| nmrPipe  -fn REV 				                   \
| nmrPipe  -fn TP                                      \
     -verb -ov -out dif.ft2

addNMR -in1 sum.ft2 -in2 dif.ft2 -out test.ft2 -add
rm sum.ft2 dif.ft2
