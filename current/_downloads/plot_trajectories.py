import numpy as np
import matplotlib.pyplot as plt

# exponential function used to fit the data
def fit_func(p,x):
    A, R2 = p
    return A * np.exp(-1.0 * np.array(x) * R2 / 1.0e6)

fitting_results = np.recfromtxt('fits.txt')
experimental_relaxation_times = np.loadtxt("relaxation_times.in")
simulated_relaxation_times = np.linspace(0,4000000,2000)

# loop over the fitting results
for peak, A, R2, ier in fitting_results:

    print "Plotting:", peak
    
    # load the experimental and simulated relaxation trajectories
    experimental_trajectory = np.loadtxt(peak + '.dat')
    simulated_trajectory = fit_func((A, R2), simulated_relaxation_times)
    
    # create the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(experimental_relaxation_times, experimental_trajectory, 'or')
    ax.plot(simulated_relaxation_times, simulated_trajectory, '-k')
    ax.set_title(peak)
    fig.savefig(peak+"_plot.png")
