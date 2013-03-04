Relaxation trajectory analysis example
======================================

Introduction
------------

This example is taken from Listing S12 - S15 in the 2013 JBNMR nmrglue paper.
In this example a series of 3D NMRPipe files containing relaxation trajectories
for a solid state NMR experment and analyzed.  



Instructions
------------

Execute `python extract_traj.py` to extract the relaxation trajectories from
the data set.  'XXX.dat' files are created for the peaks defined in the
`boxes.in` file.  The `spectra.in` file defines which spectra the trajectories
will be extracted from.

Execute `python plot_boxes.py` to create plots showing the peak and 
integration limits for all peaks defined in `boxes.in`. 
`peak_XXX_spectrum_X.png` files are created for all peaks and spectra.  A
example plot is provided as Figure 7 of the paper, which corresponds to 
`peak_D40_spectrum_0.png`

Execute `python fit_exp.py` to fit all relaxation trajectories.  The fitting
results are provided in the `fits.txt` file.  The `relaxation_times.in` file
defines the relaxation times for each spectra.

Execute `python plot_trajectories` to create plots of all experimental and fit
relaxation tranjectories.  This script creates a series of `XXX_plot.png` 
files. An example plot is provided as Figure 8 of the paper, which corresponds
to `D40_plot.png`.

