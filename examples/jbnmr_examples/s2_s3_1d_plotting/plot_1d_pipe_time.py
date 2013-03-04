import nmrglue as ng
import matplotlib.pyplot as plt

# read in the data from a NMRPipe file
dic, data = ng.pipe.read("test.fid")

# make a unit conversion object for the axis
uc = ng.pipe.make_uc(dic, data)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(uc.ms_scale(), data.real, 'k-')

# decorate axes
ax.set_yticklabels([])
ax.set_xlabel("Time (ms)")
ax.set_ylim(-100000, 100000)

# save the figure
fig.savefig("fid.png")
