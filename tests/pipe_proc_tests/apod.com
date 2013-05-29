#! /bin/csh

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName SP -q1 0.35 -q2 0.98 -q3 2.0 \
-ov -out apod1.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName EM -q1 10.0 \
-ov -out apod2.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName GM -q1 2.35 -q2 1.25 -q3 1.2 \
-ov -out apod3.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName GMB -q1 1.25 -q2 3.00 \
-ov -out apod4.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName TM -q1 100 -q2 200 \
-ov -out apod5.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName TRI -q1 500 -q2 0.50 -q3 0.8 \
-ov -out apod6.dat

nmrPipe -in ./time_complex.fid                         \
| nmrPipe -fn APOD -qName JMOD -q1 5.0 -q2 2.5 -q3 1.2 \
-ov -out apod7.dat


