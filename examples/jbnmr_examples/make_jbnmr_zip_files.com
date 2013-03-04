#! /bin/csh

# this file creates example zip files which should be uploaded to the
# nmrglue google code page and linked to the examples

mkdir zip_files

# Listing S1 sparky_to_nmrpipe
zip jbnmr_s1_sparky_to_nmrpipe.zip             \
s1_sparky_to_nmrpipe/README.rst                \
s1_sparky_to_nmrpipe/sparky2pipe_convert.py    \
s1_sparky_to_nmrpipe/data.ucsf

mv jbnmr_s1_sparky_to_nmrpipe.zip ./zip_files/.

# Listing S2 and S3 1D plotting
zip jbnmr_s2_s3_1d_plotting.zip                 \
s2_s3_1d_plotting/plot_1d_pipe_spectrum.py      \
s2_s3_1d_plotting/plot_1d_pipe_time.py          \
s2_s3_1d_plotting/README.rst                    \
s2_s3_1d_plotting/test.fid                      \
s2_s3_1d_plotting/test.ft                       

mv jbnmr_s2_s3_1d_plotting.zip ./zip_files/.

# Listing S4 2d plotting
zip jbnmr_s4_2d_plotting.zip                \
s4_2d_plotting/plot_2d_pipe_spectrum.py     \
s4_2d_plotting/README.rst                   \
s4_2d_plotting/test.ft2

mv jbnmr_s4_2d_plotting.zip ./zip_files/.


# Listing S5 seperate interleaved
zip jbnmr_s5_seperate_interleaved.zip               \
s5_seperate_interleaved/seperate.py                 \
s5_seperate_interleaved/README.rst                  \
s5_seperate_interleaved/arrayed_data.fid            \
s5_seperate_interleaved/arrayed_data.fid/procpar    \
s5_seperate_interleaved/arrayed_data.fid/text       \
s5_seperate_interleaved/arrayed_data.fid/log        \
s5_seperate_interleaved/arrayed_data.fid/fid

mv jbnmr_s5_seperate_interleaved.zip ./zip_files/.

# Listing S6 seperate inner_phase
zip jbnmr_s6_seperate_inner_phase.zip               \
s6_seperate_inner_phase/seperate.py                 \
s6_seperate_inner_phase/README.rst                  \
s6_seperate_inner_phase/arrayed_data.fid            \
s6_seperate_inner_phase/arrayed_data.fid/procpar    \
s6_seperate_inner_phase/arrayed_data.fid/text       \
s6_seperate_inner_phase/arrayed_data.fid/log        \
s6_seperate_inner_phase/arrayed_data.fid/fid

mv jbnmr_s6_seperate_inner_phase.zip ./zip_files/.


# Listing S7, S8 and S9 processing S3E filtered data
zip jbnmr_s7-s9_s3e_processing.zip      \
s7-s9_s3e_processing/convert.py         \
s7-s9_s3e_processing/fid                \
s7-s9_s3e_processing/procpar            \
s7-s9_s3e_processing/README.rst         \
s7-s9_s3e_processing/seperate_s3e.py    \
s7-s9_s3e_processing/xy_s3e.py

mv jbnmr_s7-s9_s3e_processing.zip ./zip_files/.


# Listing S10 covariance processing
zip jbnmr_s10_covariance_processing.zip     \
s10_covariance_processing/cov_process.py    \
s10_covariance_processing/test.fid          \
s10_covariance_processing/x.com             \
s10_covariance_processing/cov_plot.py       \
s10_covariance_processing/README.rst        \
s10_covariance_processing/test.ft

mv jbnmr_s10_covariance_processing.zip ./zip_files/.

# Listing S11 strip plots
zip jbnmr_s11_strip_plots_part1.zip                 \
s11_strip_plots/README.rst                          \
s11_strip_plots/make_strip_plots.py                 \
s11_strip_plots/ass.tab                             \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/fid.com   \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/procpar   \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/text      \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/xyz.com   \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/log       \
s11_strip_plots/GB3-CN-3DCONCA-041007.fid/fid       \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/fid.com   \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/procpar   \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/text      \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/xyz.com   \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/log       \
s11_strip_plots/GB3-CN-3DNCOCX-040807.fid/fid       

zip jbnmr_s11_strip_plots_part2.zip                 \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/fid.com   \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/procpar   \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/text      \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/xyz.com   \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/log       \
s11_strip_plots/GB3-CN-3DNCACX-040507.fid/fid

mv jbnmr_s11_strip_plots_part1.zip ./zip_files/.
mv jbnmr_s11_strip_plots_part2.zip ./zip_files/.


# Listing S12-S15 strip plots
zip jbnmr_s12-s15_relaxation_analysis_part1.zip             \
s12-s15_relaxation_analysis/README.rst                      \
s12-s15_relaxation_analysis/fit_exp.py                      \
s12-s15_relaxation_analysis/extract_traj.py                 \
s12-s15_relaxation_analysis/plot_boxes.py                   \
s12-s15_relaxation_analysis/plot_trajectories.py            \
s12-s15_relaxation_analysis/spectra.in                      \
s12-s15_relaxation_analysis/boxes.in                        \
s12-s15_relaxation_analysis/relaxation_times.in             \
s12-s15_relaxation_analysis/data/Ytau_100.fid/test.ft2      \
s12-s15_relaxation_analysis/data/Ytau_100000.fid/test.ft2   \
s12-s15_relaxation_analysis/data/Ytau_250000.fid/test.ft2   \

zip jbnmr_s12-s15_relaxation_analysis_part2.zip             \
s12-s15_relaxation_analysis/data/Ytau_500000.fid/test.ft2   \
s12-s15_relaxation_analysis/data/Ytau_750000.fid/test.ft2   \
s12-s15_relaxation_analysis/data/Ytau_1000000.fid/test.ft2

zip jbnmr_s12-s15_relaxation_analysis_part3.zip             \
s12-s15_relaxation_analysis/data/Ytau_1500000.fid/test.ft2  \
s12-s15_relaxation_analysis/data/Ytau_2000000.fid/test.ft2  \
s12-s15_relaxation_analysis/data/Ytau_3000000.fid/test.ft2  

zip jbnmr_s12-s15_relaxation_analysis_part4.zip             \
s12-s15_relaxation_analysis/data/Ytau_4000000.fid/test.ft2

mv jbnmr_s12-s15_relaxation_analysis_part1.zip ./zip_files/.
mv jbnmr_s12-s15_relaxation_analysis_part2.zip ./zip_files/.
mv jbnmr_s12-s15_relaxation_analysis_part3.zip ./zip_files/.
mv jbnmr_s12-s15_relaxation_analysis_part4.zip ./zip_files/.
