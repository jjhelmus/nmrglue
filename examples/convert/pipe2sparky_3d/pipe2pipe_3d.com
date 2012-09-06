#!/bin/sh

xyz2pipe -in ./nmrpipe_3d/data/test%03d.fid -x -verb \
| nmrPipe -fn ZF -auto                               \
| nmrPipe -fn FT -auto                               \
| nmrPipe -fn PS -p0 0.0 -p1 0.0  -di                \
| nmrPipe -fn TP                                     \
| nmrPipe -fn ZF -auto                               \
| nmrPipe -fn FT -auto                               \
| nmrPipe -fn PS -p0 0 -p1 0  -di                    \
| pipe2xyz -out ./nmrpipe_3d/ft/test%03d.ft2 -y

xyz2pipe -in ./nmrpipe_3d/ft/test%03d.ft2 -z -verb   \
| nmrPipe -fn ZF -auto                               \
| nmrPipe -fn FT -auto                               \
| nmrPipe -fn PS -p0 -92 -p1 65  -di                 \
| pipe2xyz -out ./nmrpipe_3d/ft/test%03d.ft3 -z		    
