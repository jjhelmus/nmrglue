#! /usr/bin/env python
# Plot trajectories and fitting results

import pickle
import matplotlib.pyplot as plt
import numpy as np

# the same fit_func as in fit_exp_leastsq.py
def fit_func(p, x):

    A, R2 = p

    # bound A between 0.98 and 1.02 (although fits do not reflect this)
    if A > 1.02:
        A = 1.02
    if A < 0.98:
        A = 0.98

    return A * np.exp(-1.0 * np.array(x) * R2 / 1.0e6)

# read in the trajectories, fitting results, and times
fits = pickle.load(open("fits.pickle"))
trajs = np.load("traj.npy")
times = np.recfromtxt("time.dat")

sim_times = np.linspace(times[0], times[-1], 2000)

# loop over the peaks
for peak, params in fits.iteritems():

    print "Plotting:", peak
    exp_traj = trajs[peak]
    sim_traj = fit_func(params[0], sim_times)

    # create the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(times, exp_traj, 'or')
    ax.plot(sim_times, sim_traj, '-k')
    ax.set_title(peak)

    # save the figure
    fig.savefig(peak + "_plot.png")
