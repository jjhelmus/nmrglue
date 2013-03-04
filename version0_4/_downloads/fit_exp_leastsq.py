#! /usr/bin/env python
# fit a collection to T1 trajectories to a decaying exponential

import scipy.optimize
import numpy as np
import pickle

# read in the trajectories and times
trajs = np.load("traj.npy")
t1 = np.recfromtxt("time.dat")

# fitting function and residual calculation
def fit_func(p, x):
    A, R2 = p
    # bound A between 0.98 and 1.02 (although fits do not reflect this)
    if A > 1.02:
        A = 1.02
    if A < 0.98:
        A = 0.98

    return A * np.exp(-1.0 * np.array(x) * R2 / 1.0e6)

def residuals(p, y, x):
    err = y - fit_func(p, x)
    return err

p0 = [1.0, 0.05] # initial guess

fits = {}
# loop over the peak trajectories
for peak in trajs.dtype.names:

    print "Fitting Peak:", peak

    # get the trajectory to fit
    traj = trajs[peak]

    # fit the trajectory using leastsq (fmin, etc can also be used)
    results = scipy.optimize.leastsq(residuals, p0, args=(traj, t1))
    fits[peak] = results

# pickle the fits 
f = open("fits.pickle", 'w')
pickle.dump(fits, f)
f.close()

# output the fits nicely to file
f = open("fits.txt", 'w')
f.write("#Peak\tA\t\tR2\t\tier\n")
for k, v in fits.iteritems():
    f.write(k + "\t" + str(v[0][0]) + "\t" + str(v[0][1]) + "\t" + str(v[1]) + 
            "\n")
f.close()
