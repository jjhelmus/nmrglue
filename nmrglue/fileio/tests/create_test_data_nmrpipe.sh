#!/bin/csh

# nmrpipe_1d_time.fid
# 1D time containing 1 - 1j, 2 - 2j, 0, 0, ... 
simTimeND                \
  -xN                32  \
  -xT                16  \
  -xMODE        Complex  \
  -xSW        50000.000  \
  -xOBS         500.000  \
  -xCAR          99.000  \
  -xLAB              H1  \
  -ndim               1  \
  -out data/nmrpipe_1d_time.fid -ov

nmrPipe -in data/nmrpipe_1d_time.fid         \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1 \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2 \
   -ov -out data/nmrpipe_1d_time.fid

# nmrpipe_1d_freq.ft2
# 1D NMRPipe freq containing 1, 2
nmrPipe -in data/nmrpipe_1d_time.fid         \
| nmrPipe -fn SET                            \
| nmrPipe -fn FT -di                         \
| nmrPipe -fn SET -r 1.0 -x1 1 -xn 1         \
| nmrPipe -fn SET -r 2.0 -x1 2 -xn 2         \
   -ov -out data/nmrpipe_1d_freq.fid

# nmrpipe_1d_ext.ft2
# 1D NMRPipe freq extracted containing 1, 2
nmrPipe -in data/nmrpipe_1d_freq.fid         \
| nmrPipe -fn SET                            \
| nmrPipe -fn EXT -x1 4 -xn 11 -sw           \
| nmrPipe -fn SET -r 1.0 -x1 1 -xn 1         \
| nmrPipe -fn SET -r 2.0 -x1 2 -xn 2         \
   -ov -out data/nmrpipe_1d_ext.fid

# nmrpipe_2d_time.fid
# 2D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
simTimeND                                      \
  -xN                16  -yN                 4 \
  -xT                 8  -yT                 2 \
  -xMODE        Complex  -yMODE         States \
  -xSW        50000.000  -ySW        20000.000 \
  -xOBS         500.000  -yOBS         125.000 \
  -xCAR           4.700  -yCAR          99.000 \
  -xLAB              H1  -yLAB             C13 \
  -ndim               2  -aq2D          States \
  -out data/nmrpipe_2d_time.fid -ov

nmrPipe -in data/nmrpipe_2d_time.fid         \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1 \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2 \
   -ov -out data/nmrpipe_2d_time.fid

# nmrpipe_2d_freq.ft2
# 2D NMRPipe freq containing 1, 2
nmrPipe -in data/nmrpipe_2d_time.fid         \
| nmrPipe -fn SET                            \
| nmrPipe -fn FT -di                         \
| nmrPipe -fn TP                             \
| nmrPipe -fn FT -di                         \
| nmrPipe -fn TP                             \
| nmrPipe -fn SET -r 1.0 -x1 1 -xn 1         \
| nmrPipe -fn SET -r 2.0 -x1 2 -xn 2         \
   -ov -out data/nmrpipe_2d_freq.ft2

# nmrpipe_2d_time_tp.fid
# 2D NMRPipe time transposed
nmrPipe -in data/nmrpipe_2d_time.fid         \
| nmrPipe -fn TP -hyper                      \
   -ov -out data/nmrpipe_2d_time_tp.fid

# nmrpipe_2d_freq_tp.ft2
# 2D NMRPipe freq transposed
nmrPipe -in data/nmrpipe_2d_freq.ft2         \
| nmrPipe -fn TP                             \
   -ov -out data/nmrpipe_2d_freq_tp.ft2


# nmrpipe_3d_time.fid
# 3D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
simTimeND \
  -xN             16  -yN            6 -zN            4 \
  -xT              8  -yT            3 -zT            2 \
  -xMODE     Complex  -yMODE    States -zMODE    States \
  -xSW     50000.000  -ySW   20000.000 -zSW   10000.000 \
  -xOBS      500.000  -yOBS    125.000 -zOBS     50.000 \
  -xCAR        4.700  -yCAR     99.000 -zCAR    120.000 \
  -xLAB           H1  -yLAB        C13 -zLAB        N15 \
  -ndim            3  -aq2D          States             \
  -out data/nmrpipe_3d_time.dir/nmrpipe_3d_time_%03d.fid -ov

xyz2pipe -in data/nmrpipe_3d_time.dir/nmrpipe_3d_time_%03d.fid \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1 \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2 \
| pipe2xyz -out data/nmrpipe_3d_time.dir/nmrpipe_3d_time_%03d.fid

# nmrpipe_3d_freq.fid
# 3D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
xyz2pipe -in data/nmrpipe_3d_time.dir/nmrpipe_3d_time_%03d.fid -x   \
| nmrPipe -fn SET                                                   \
| nmrPipe -fn FT -di                                                \
| nmrPipe -fn TP                                                    \
| nmrPipe -fn FT -di                                                \
| pipe2xyz -out data/nmrpipe_3d_freq.dir/nmrpipe_3d_freq_%03d.ft2 -y

xyz2pipe -in data/nmrpipe_3d_freq.dir/nmrpipe_3d_freq_%03d.ft2 -z   \
| nmrPipe -fn FT -di                                                \
| nmrPipe -fn SET -r 1.0 -x1 1 -xn 1                                \
| nmrPipe -fn SET -r 2.0 -x1 2 -xn 2                                \
| pipe2xyz -out data/nmrpipe_3d_freq.dir/nmrpipe_3d_freq_%03d.ft3 -z

rm data/nmrpipe_3d_freq.dir/*.ft2

# nmrpipe_3d_time.fid
# 3D NMRPipe stream file (single file 3D file)
xyz2pipe -in data/nmrpipe_3d_time.dir/nmrpipe_3d_time_%03d.fid -x > data/nmrpipe_3d_time.fid

xyz2pipe -in data/nmrpipe_3d_freq.dir/nmrpipe_3d_freq_%03d.ft3 -x > data/nmrpipe_3d_freq.ft3

# nmrpipe_4d_time_2.dir
# 4D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
simTimeND \
  -xN           10  -yN            8 -zN            6 -aN            4 \
  -xT            5  -yT            4 -zT            3 -aT            2 \
  -xMODE   Complex  -yMODE    States -zMODE    States -aMODE    States \
  -xSW   50000.000  -ySW   20000.000 -zSW   10000.000 -aSW   30000.000 \
  -xOBS    500.000  -yOBS    125.000 -zOBS     50.000 -aOBS    150.000 \
  -xCAR      4.700  -yCAR     99.000 -zCAR    120.000 -aCAR     80.000 \
  -xLAB         H1  -yLAB        C13 -zLAB        N15 -aLAB        P31 \
  -ndim          4  -aq2D     States                                   \
| pipe2xyz -out data/nmrpipe_4d_time_2.dir/nmrpipe_4d_time_%03d_%03d.fid -x -ov -to 0

xyz2pipe -in data/nmrpipe_4d_time_2.dir/nmrpipe_4d_time_%03d_%03d.fid -x       \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1 \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2 \
| pipe2xyz -out data/nmrpipe_4d_time_2.dir/nmrpipe_4d_time_%03d_%03d.fid -x -ov

# nmrpipe_4d_time_1.dir
# 4D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
simTimeND \
  -xN           10  -yN            8 -zN            6 -aN            4 \
  -xT            5  -yT            4 -zT            3 -aT            2 \
  -xMODE   Complex  -yMODE    States -zMODE    States -aMODE    States \
  -xSW   50000.000  -ySW   20000.000 -zSW   10000.000 -aSW   30000.000 \
  -xOBS    500.000  -yOBS    125.000 -zOBS     50.000 -aOBS    150.000 \
  -xCAR      4.700  -yCAR     99.000 -zCAR    120.000 -aCAR     80.000 \
  -xLAB         H1  -yLAB        C13 -zLAB        N15 -aLAB        P31 \
  -ndim          4  -aq2D     States                                   \
  -out data/nmrpipe_4d_time_1.dir/nmrpipe_4d_time_%03d.fid -verb -ov

xyz2pipe -in data/nmrpipe_4d_time_1.dir/nmrpipe_4d_time_%03d.fid -x       \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1 \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2 \
| pipe2xyz -out data/nmrpipe_4d_time_1.dir/nmrpipe_4d_time_%03d.fid -x -ov

# nmrpipe_4d_time.fid
# 4D Stream file
xyz2pipe -in data/nmrpipe_4d_time_2.dir/nmrpipe_4d_time_%03d_%03d.fid -x > data/nmrpipe_4d_time.fid


# nmrpipe_4d_freq_1.dir
# 4D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
xyz2pipe -in data/nmrpipe_4d_time_1.dir/nmrpipe_4d_time_%03d.fid -x -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft1 -x -ov

xyz2pipe -in data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft1 -y -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft2 -y -ov

xyz2pipe -in data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft2 -z -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft3 -z -ov

xyz2pipe -in data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft3 -a -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1                                \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2                                \
| pipe2xyz -out data/nmrpipe_4d_freq_1.dir/nmrpipe_4d_freq_%03d.ft4 -a -ov

rm data/nmrpipe_4d_freq_1.dir/*.ft1
rm data/nmrpipe_4d_freq_1.dir/*.ft2
rm data/nmrpipe_4d_freq_1.dir/*.ft3

# nmrpipe_4d_freq_2.dir
# 4D NMRPipe time containing 1 - 1j, 2 - 2j, 0 + 0j along each vector
xyz2pipe -in data/nmrpipe_4d_time_2.dir/nmrpipe_4d_time_%03d_%03d.fid -x -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft1 -x -ov

xyz2pipe -in data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft1 -y -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft2 -y -ov

xyz2pipe -in data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft2 -z -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| pipe2xyz -out data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft3 -z -ov

xyz2pipe -in data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft3 -a -verb   \
| nmrPipe -fn SET                                                           \
| nmrPipe -fn FT -di                                                        \
| nmrPipe -fn SET -r 1.0 -i -1.0 -x1 1 -xn 1                                \
| nmrPipe -fn SET -r 2.0 -i -2.0 -x1 2 -xn 2                                \
| pipe2xyz -out data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft4 -a -ov

rm data/nmrpipe_4d_freq_2.dir/*.ft1
rm data/nmrpipe_4d_freq_2.dir/*.ft2
rm data/nmrpipe_4d_freq_2.dir/*.ft3

# nmrpipe_4d_freq.ft4
# 4D Stream file
xyz2pipe -in data/nmrpipe_4d_freq_2.dir/nmrpipe_4d_freq_%03d_%03d.ft4 -x > data/nmrpipe_4d_freq.ft4
