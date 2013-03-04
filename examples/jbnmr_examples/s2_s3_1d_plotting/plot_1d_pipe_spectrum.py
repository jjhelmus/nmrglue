import nmrglue as ng
import matplotlib.pyplot as plt

# read in the data from a NMRPipe file
dic, data = ng.pipe.read("test.ft")

# create a unit conversion object for the axis
uc = ng.pipe.make_uc(dic, data)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(uc.ppm_scale(), data, 'k-')

# decorate axes
ax.set_yticklabels([])
ax.set_xlabel("13C ppm")
ax.set_xlim(200, 0)
ax.set_ylim(-80000, 2500000)

# save the figure
fig.savefig("spectrum.png")
