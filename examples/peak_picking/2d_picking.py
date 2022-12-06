import nmrglue as ng
import matplotlib.pyplot as plt

# read data
dic, data = ng.pipe.read("test.ft2")
uc = {
    "n": ng.pipe.make_uc(dic, data, dim=0),
    "c": ng.pipe.make_uc(dic, data, dim=1),
}

threshold = 8.5e4


# detect all peaks with a threshold
peaks = ng.peakpick.pick(data, pthres=threshold, algorithm="thres", msep=[1, 1])


# contour levels and axes
ext = [*uc["c"].ppm_limits(), *uc["n"].ppm_limits()]
clevs = [threshold * 1.4 ** i for i in range(10)]


# plot and indicate all peaks
fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
ax.contour(data.real, levels=clevs, extent=ext, cmap="bone", linewidths=0.5)
ax.set_xlim(183, 168)
ax.set_ylim(135, 100)


# add markers for peak positions
x = uc["c"].ppm(peaks["X_AXIS"])
y = uc["n"].ppm(peaks["Y_AXIS"])
ax.scatter(x, y, marker=".", s=10, color="r")


ax.set_xlabel("$^{13}$C (ppm)")
ax.set_ylabel("$^{15}$N (ppm)")

plt.show()
