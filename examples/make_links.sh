#! /bin/sh

# FULL path to data
set DATADIR = /home/jhelmus/sandbox/nmrglue_test/nmrglue/data

ln -sf $DATADIR/agilent_1d      convert/agilent2pipe_1d/agilent_1d
ln -sf $DATADIR/agilent_1d      convert/agilent2pipe_1d/agilent_1d
ln -sf $DATADIR/agilent_2d      convert/agilent2pipe_2d/agilent_2d
ln -sf $DATADIR/agilent_2d_tppi convert/agilent2pipe_2d_tppi/agilent_2d_tppi
ln -sf $DATADIR/agilent_3d      convert/agilent2pipe_3d/agilent_3d
ln -sf $DATADIR/bruker_1d       convert/bruker2pipe_1d/bruker_1d
ln -sf $DATADIR/bruker_2d       convert/bruker2pipe_2d/bruker_2d
ln -sf $DATADIR/bruker_3d       convert/bruker2pipe_3d/bruker_3d
ln -sf $DATADIR/nmrpipe_2d      convert/pipe2sparky_2d/nmrpipe_2d
ln -sf $DATADIR/nmrpipe_3d      convert/pipe2sparky_3d/nmrpipe_3d
ln -sf $DATADIR/nmrpipe_2d      integration/integrate_2d/nmrpipe_2d
ln -sf $DATADIR/nmrpipe_1d      interactive/interactive_1d/nmrpipe_1d
ln -sf $DATADIR/nmrpipe_2d      interactive/interactive_2d/nmrpipe_2d
ln -sf $DATADIR/nmrpipe_1d      plotting/plot_1d/nmrpipe_1d
ln -sf $DATADIR/nmrpipe_2d      plotting/plot_2d/nmrpipe_2d
ln -sf $DATADIR/nmrpipe_1d      processing/nmrpipe_1d
ln -sf $DATADIR/nmrpipe_2d      processing/nmrpipe_2d
ln -sf $DATADIR/nmrpipe_2d_tppi processing/nmrpipe_2d_tppi
ln -sf $DATADIR/nmrpipe_3d      processing/nmrpipe_3d
