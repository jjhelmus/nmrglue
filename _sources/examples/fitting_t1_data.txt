.. _fitting_t1_data:

fitting example: fitting_t1_data
================================

This example shows how to use nmrglue and `scipy <http://www.scipy.org/>`_ 
optimize module to fit T1 relaxation trajectories.  Three scripts are used in
the process.

First the ``extract_trajs.py`` script reads in box limits from ``boxes.in`` and
a list of spectra from ``spectra.in``.  The script integrates each peak in each
spectrum and writies the trajectory for each peak to disk as ``traj.npy`` in 
`numpy <http://numpy.scipy.org/>`_ ``.npy`` format.

[`extract_trajs.py <el/fitting_data/t1_measurements/extract_trajs.py>`_]

.. literalinclude:: el/fitting_data/t1_measurements/extract_trajs.py

[`boxes.in <el/fitting_data/t1_measurements/boxes.in>`_]

.. literalinclude:: el/fitting_data/t1_measurements/boxes.in

[`spectra.in <el/fitting_data/t1_measurements/spectra.in>`_]

.. literalinclude:: el/fitting_data/t1_measurements/spectra.in


The second script ``fit_exp_leastsq.py`` reads in this ``traj.npy`` file and the
T1 relaxation times associated with the spectra collected from ``time.dat``.  
Each trajectory is fit using the least squares approach. Other optimization
algorithms can be substituted with small changes to the code, see the 
`scipy.optimize <http://docs.scipy.org/doc/scipy/reference/optimize.html>_`
documentation).  The resulting fits are saved to a `fits.pickle` file for 
easy reading into python as well as the human readable ``fits.txt`` file. 


[`fit_exp_leastsq.py <el/fitting_data/t1_measurements/fit_exp_leastsq.py>`_]

.. literalinclude:: el/fitting_data/t1_measurements/fit_exp_leastsq.py

[`time.dat <el/fitting_data/t1_measurements/time.dat>`_]

.. literalinclude:: el/fitting_data/t1_measurements/time.dat

Results: 

[`fits.txt <el/fitting_data/t1_measurements/fits.txt>`_]

.. literalinclude:: el/fitting_data/t1_measurements/fits.txt

The last script ``pt.py`` reads in the fits, trajectories and T1 
relaxation times and plots the experimental points and best fit to a series 
of ``*_plot.png`` files.


[`pt.py <el/fitting_data/t1_measurements/pt.py>`_]

.. literalinclude:: el/fitting_data/t1_measurements/pt.py

Results:

[`A24_plot.png <el/fitting_data/t1_measurements/A24_plot.png>`_]

.. image:: el/fitting_data/t1_measurements/A24_plot.png
