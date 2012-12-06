#! /bin/csh

# this file creates example zip files which should be uploaded to the
# nmrglue google code page and linked to the examples

mkdir zip_files

#####################
# coadding examples #
#####################
cd coadd

# coadd_1d_pipe
zip example_coadd_1d_pipe.zip \
coadd_1d_pipe/coadd_1d.py     \
coadd_1d_pipe/README          \
coadd_1d_pipe/test-run01.fid  \
coadd_1d_pipe/test-run02.fid  \
coadd_1d_pipe/test-run03.fid  \
coadd_1d_pipe/test-run04.fid

mv example_coadd_1d_pipe.zip ../zip_files/.

# coadd_2d_pipe
zip example_coadd_2d_pipe.zip \
coadd_2d_pipe/coadd_2d.py     \
coadd_2d_pipe/README          \
coadd_2d_pipe/test-run01.fid  \
coadd_2d_pipe/test-run02.fid  \
coadd_2d_pipe/test-run03.fid

mv example_coadd_2d_pipe.zip ../zip_files/.

# coadd_pseudo3d_pipe
zip example_coadd_pseudo3d_pipe.zip                     \
coadd_pseudo3d_pipe/coadd_pseudo3d.py                   \
coadd_pseudo3d_pipe/README                              \
coadd_pseudo3d_pipe/run02.dir/nredor_20.fid/test.fid    \
coadd_pseudo3d_pipe/run02.dir/nredor_50.fid/test.fid    \
coadd_pseudo3d_pipe/run02.dir/nredor_60.fid/test.fid    \
coadd_pseudo3d_pipe/run02.dir/nredor_40.fid/test.fid    \
coadd_pseudo3d_pipe/run01.dir/nredor_20.fid/test.fid    \
coadd_pseudo3d_pipe/run01.dir/nredor_50.fid/test.fid    \
coadd_pseudo3d_pipe/run01.dir/nredor_60.fid/test.fid    \
coadd_pseudo3d_pipe/run01.dir/nredor_40.fid/test.fid    \

mv example_coadd_pseudo3d_pipe.zip ../zip_files/.

cd ..

############################
# data conversion examples #
############################
cd convert

# agilent2pipe_1d
zip example_agilent2pipe_1d.zip      \
agilent2pipe_1d/agilent2pipe_1d.py   \
agilent2pipe_1d/agilent2pipe_1d.com  \
agilent2pipe_1d/compare.py           \
agilent2pipe_1d/README               \
agilent2pipe_1d/agilent_1d/fid       \
agilent2pipe_1d/agilent_1d/procpar

mv example_agilent2pipe_1d.zip ../zip_files/.

# agilent2pipe_2d
zip example_agilent2pipe_2d.zip      \
agilent2pipe_2d/agilent2pipe_2d.py   \
agilent2pipe_2d/agilent2pipe_2d.com  \
agilent2pipe_2d/compare.py           \
agilent2pipe_2d/README               \
agilent2pipe_2d/agilent_2d/fid       \
agilent2pipe_2d/agilent_2d/procpar

mv example_agilent2pipe_2d.zip ../zip_files/.

# agilent2pipe_2d_tppi
zip example_agilent2pipe_2d_tppi.zip            \
agilent2pipe_2d_tppi/agilent2pipe_2d_tppi.py    \
agilent2pipe_2d_tppi/agilent2pipe_2d_tppi.com   \
agilent2pipe_2d_tppi/compare.py                 \
agilent2pipe_2d_tppi/README                     \
agilent2pipe_2d_tppi/agilent_2d_tppi/fid        \
agilent2pipe_2d_tppi/agilent_2d_tppi/procpar

mv example_agilent2pipe_2d_tppi.zip ../zip_files/.

# agilent2pipe_3d
zip example_agilent2pipe_3d.zip      \
agilent2pipe_3d/agilent2pipe_3d.py   \
agilent2pipe_3d/agilent2pipe_3d.com  \
agilent2pipe_3d/compare.py           \
agilent2pipe_3d/README               \
agilent2pipe_3d/agilent_3d/fid       \
agilent2pipe_3d/agilent_3d/procpar

mv example_agilent2pipe_3d.zip ../zip_files/.

# bruker2pipe_1d
zip example_bruker2pipe_1d.zip      \
bruker2pipe_1d/bruker2pipe_1d.py    \
bruker2pipe_1d/bruker2pipe_1d.com   \
bruker2pipe_1d/compare.py           \
bruker2pipe_1d/README               \
bruker2pipe_1d/bruker_1d/acqu       \
bruker2pipe_1d/bruker_1d/acqus      \
bruker2pipe_1d/bruker_1d/fid        \
bruker2pipe_1d/bruker_1d/pulseprogram

mv example_bruker2pipe_1d.zip ../zip_files/.

# bruker2pipe_2d
zip example_bruker2pipe_2d.zip      \
bruker2pipe_2d/bruker2pipe_2d.py    \
bruker2pipe_2d/bruker2pipe_2d.com   \
bruker2pipe_2d/compare.py           \
bruker2pipe_2d/README               \
bruker2pipe_2d/bruker_2d/acqu       \
bruker2pipe_2d/bruker_2d/acqus      \
bruker2pipe_2d/bruker_2d/acqu2      \
bruker2pipe_2d/bruker_2d/acqu2s     \
bruker2pipe_2d/bruker_2d/ser        \
bruker2pipe_2d/bruker_2d/pulseprogram

mv example_bruker2pipe_2d.zip ../zip_files/.

# bruker2pipe_3d
zip example_bruker2pipe_3d.zip      \
bruker2pipe_3d/bruker2pipe_3d.py    \
bruker2pipe_3d/bruker2pipe_3d.com   \
bruker2pipe_3d/compare.py           \
bruker2pipe_3d/README               \
bruker2pipe_3d/bruker_3d/acqu       \
bruker2pipe_3d/bruker_3d/acqus      \
bruker2pipe_3d/bruker_3d/acqu2      \
bruker2pipe_3d/bruker_3d/acqu2s     \
bruker2pipe_3d/bruker_3d/ser        \
bruker2pipe_3d/bruker_3d/pulseprogram

mv example_bruker2pipe_3d.zip ../zip_files/.

# pipe2sparky_2d
zip example_pipe2sparky_2d.zip      \
pipe2sparky_2d/pipe2sparky_2d.py    \
pipe2sparky_2d/pipe2sparky_2d.com   \
pipe2sparky_2d/README               \
pipe2sparky_2d/compare.py           \
pipe2sparky_2d/nmrpipe_2d/test.ft2

mv example_pipe2sparky_2d.zip ../zip_files/.

# pipe2sparky_3d
zip example_pipe2sparky_3d.zip      \
pipe2sparky_3d/pipe2sparky_3d.py    \
pipe2sparky_3d/pipe2sparky_3d.com   \
pipe2sparky_3d/README               \
pipe2sparky_3d/compare.py           \
pipe2sparky_3d/nmrpipe_3d/data/test*.fid

mv example_pipe2sparky_3d.zip ../zip_files/.

cd ..

#####################
# plotting examples #
#####################

cd plotting

# plot_1d
zip example_plot_1d.zip         \
plot_1d/plot_1d_pipe_freq.py    \
plot_1d/plot_1d_pipe_time.py    \
plot_1d/README                  \
plot_1d/nmrpipe_1d/test.fid     \
plot_1d/nmrpipe_1d/test.ft      \

mv example_plot_1d.zip ../zip_files/.

# plot_2d
zip example_plot_2d.zip         \
plot_2d/plot_spectrum.py        \
plot_2d/plot_spectrum_pts.py    \
plot_2d/plot_assignments.py     \
plot_2d/plot_boxes.py           \
plot_2d/limits.in               \
plot_2d/README                  \
plot_2d/nmrpipe_2d/test.ft2     \

mv example_plot_2d.zip ../zip_files/.

cd ..

#######################
# processing examples #
#######################

# process_pipe_1d
zip example_process_pipe_1d.zip     \
processing/README                   \
processing/process_pipe_1d.py       \
processing/compare_1d.py            \
processing/pipe2pipe_1d.com         \
processing/nmrpipe_1d/test.fid      

mv example_process_pipe_1d.zip zip_files/.

# process_pipe_2d
zip example_process_pipe_2d.zip     \
processing/README                   \
processing/process_pipe_2d.py       \
processing/compare_2d.py            \
processing/pipe2pipe_2d.com         \
processing/nmrpipe_2d/test.fid

mv example_process_pipe_2d.zip zip_files/.

# process_pipe_2d_tppi
zip example_process_pipe_2d_tppi.zip    \
processing/README                       \
processing/process_pipe_2d_tppi.py      \
processing/compare_2d_tppi.py           \
processing/pipe2pipe_2d_tppi.com        \
processing/nmrpipe_2d_tppi/test.fid

mv example_process_pipe_2d_tppi.zip zip_files/.

# process_pipe_3d
zip example_process_pipe_3d.zip     \
processing/README                   \
processing/process_pipe_3d.py       \
processing/compare_3d.py            \
processing/pipe2pipe_3d.com         \
processing/nmrpipe_3d/data/*.fid

mv example_process_pipe_3d.zip zip_files/.

################################
# Interactive plotting example #
################################

cd interactive

# interactive_1d
zip example_interactive_1d.zip      \
interactive_1d/1d_interactive.py    \
interactive_1d/README               \
interactive_1d/nmrpipe_1d/test.ft

mv example_interactive_1d.zip ../zip_files/.

# interactive_2d
zip example_interactive_2d.zip      \
interactive_2d/2d_interactive.py    \
interactive_2d/README               \
interactive_2d/nmrpipe_2d/test.ft2

mv example_interactive_2d.zip ../zip_files/.

cd ..

########################
# Integration examples #
########################

cd integration

# intergrate_1d
zip example_integrate_1d.zip    \
integrate_1d/README             \
integrate_1d/integrate_1d.py    \
integrate_1d/1d_data.ft         \
integrate_1d/limits.in

mv example_integrate_1d.zip ../zip_files/.

# intergrate_2d
zip example_integrate_2d.zip        \
integrate_2d/README                 \
integrate_2d/integrate_2d.py        \
integrate_2d/nmrpipe_2d/test.ft2    \
integrate_2d/limits.in

mv example_integrate_2d.zip ../zip_files/.

cd ..

#####################
# Separate examples #
#####################

cd separate

# separate_1d_varian
zip example_separate_1d_varian.zip      \
separate_1d_varian/README               \
separate_1d_varian/separate.py          \
separate_1d_varian/arrayed_data.dir/*

mv example_separate_1d_varian.zip ../zip_files/.

# separate_2d_varian
zip example_separate_2d_varian.zip      \
separate_2d_varian/README               \
separate_2d_varian/separate.py          \
separate_2d_varian/arrayed_data.dir/*

mv example_separate_2d_varian.zip ../zip_files/.

# separate_2d_bruker
zip example_separate_2d_bruker.zip      \
separate_2d_bruker/README               \
separate_2d_bruker/separate.py          \
separate_2d_bruker/arrayed_data.dir/*

mv example_separate_2d_bruker.zip ../zip_files/.

cd ..

####################
# Fitting examples #
####################

# XXX this does not include the data as it is too large.

cd fitting_data

# fitting_t1_data
zip example_fitting_t1_data.zip     \
t1_measurements/README              \
t1_measurements/boxes.in            \
t1_measurements/extract_trajs.py    \
t1_measurements/fit_exp_leastsq.py  \
t1_measurements/pt.py               \
t1_measurements/time.dat            \
t1_measurements/spectra.in          \
t1_measurements/proc.com            \
t1_measurements/xy_s3e.com          \
t1_measurements/data/Ytau_1000000.fid/test.fid \
t1_measurements/data/Ytau_100000.fid/test.fid  \
t1_measurements/data/Ytau_100.fid/test.fid     \
t1_measurements/data/Ytau_1500000.fid/test.fid \
t1_measurements/data/Ytau_2000000.fid/test.fid \
t1_measurements/data/Ytau_250000.fid/test.fid  \
t1_measurements/data/Ytau_3000000.fid/test.fid \
t1_measurements/data/Ytau_4000000.fid/test.fid \
t1_measurements/data/Ytau_500000.fid/test.fid  \
t1_measurements/data/Ytau_750000.fid/test.fid

mv example_fitting_t1_data.zip ../zip_files/.

cd ..

#####################
# All non-test data #
#####################

# all non-test example data
zip all_none_test_example_data.zip \
coadd/coadd_2d_pipe/test-run01.fid  \
coadd/coadd_2d_pipe/test-run02.fid  \
coadd/coadd_2d_pipe/test-run03.fid  \
coadd/coadd_pseudo3d_pipe/run02.dir/nredor_20.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run02.dir/nredor_50.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run02.dir/nredor_60.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run02.dir/nredor_40.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run01.dir/nredor_20.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run01.dir/nredor_50.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run01.dir/nredor_60.fid/test.fid    \
coadd/coadd_pseudo3d_pipe/run01.dir/nredor_40.fid/test.fid    \
fitting_data/t1_measurements/data/Ytau_1000000.fid/test.fid \
fitting_data/t1_measurements/data/Ytau_100000.fid/test.fid  \
fitting_data/t1_measurements/data/Ytau_100.fid/test.fid     \
fitting_data/t1_measurements/data/Ytau_1500000.fid/test.fid \
fitting_data/t1_measurements/data/Ytau_2000000.fid/test.fid \
fitting_data/t1_measurements/data/Ytau_250000.fid/test.fid  \
fitting_data/t1_measurements/data/Ytau_3000000.fid/test.fid \
fitting_data/t1_measurements/data/Ytau_4000000.fid/test.fid \
fitting_data/t1_measurements/data/Ytau_500000.fid/test.fid  \
fitting_data/t1_measurements/data/Ytau_750000.fid/test.fid  \
separate/separate_1d_varian/arrayed_data.dir/* \
separate/separate_2d_varian/arrayed_data.dir/* \
separate/separate_2d_bruker/arrayed_data.dir/*

mv all_none_test_example_data.zip zip_files/.

