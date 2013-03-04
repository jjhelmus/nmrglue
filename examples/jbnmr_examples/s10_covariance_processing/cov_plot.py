import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

# open the data
dic, data = ng.pipe.read("test.ft")
uc = ng.pipe.make_uc(dic, data, 1)
x0, x1 = uc.ppm_limits()

# compute the covariance
C = np.cov(data.T)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
cl = [1e8 * 1.30 ** x for x in range(20)]
ax.contour(C, cl, colors='blue', extent=(x0, x1, x0, x1), origin=None)
ax.set_xlabel("13C ppm")
ax.set_xlim(62.5, 24.5)
ax.set_ylabel("13C ppm")
ax.set_ylim(62.5, 24.5)
plt.savefig('covariance_figure.png')
