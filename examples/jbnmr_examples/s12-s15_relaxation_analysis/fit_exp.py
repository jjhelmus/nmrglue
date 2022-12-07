import glob

import numpy as np
from nmrglue.analysis.leastsqbound import leastsqbound

# exponential function to fit data to.
def fit_func(p, x):
    A, R2 = p
    return A * np.exp(-1.0 * np.array(x) * R2 / 1.0e6)

# residuals between fit and experimental data.
def residuals(p, y, x):
    err = y - fit_func(p, x)
    return err

# prepare fitting parameters
relaxation_times = np.loadtxt("relaxation_times.in")
x0 = [1.0, 0.10]  # initial fitting parameter
bounds = [(0.98, 1.02), (None, None)] # fitting constraints

# create an output file to record the fitting results
output = open('fits.txt', 'w')
output.write("#Peak\tA\t\tR2\t\tier\n")

# loop over the trajecory files
for filename in glob.glob('*.dat'):

    peak = filename[:3]
    print "Fitting Peak:", peak

    # fit the trajectory using constrained least squares optimization
    trajectory = np.loadtxt(filename)
    x, ier = leastsqbound(residuals, x0, bounds=bounds,
                            args=(trajectory, relaxation_times))

    # write fitting results to output file
    output.write('%s\t%.6f\t%.6f\t%i\n' % (peak, x[0], x[1], ier))

output.close()  # close the output file
