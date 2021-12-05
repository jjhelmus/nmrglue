import nmrglue as ng
import matplotlib.pyplot as plt

# read data
dic, data = ng.pipe.read("../integration/integrate_1d/1d_data.ft")
uc = ng.pipe.make_uc(dic, data)

threshold = 1e6

# detect all peaks with a threshold
peaks = ng.peakpick.pick(data, pthres=threshold, algorithm="downward")

# plot and indicate all peaks
fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
ax.plot(uc.ppm_scale(), data.real)

# add markers for peak positions
for n, peak in enumerate(peaks):
    height = data[int(peak["X_AXIS"])]
    ppm = uc.ppm(peak["X_AXIS"])
    ax.scatter(ppm, height, marker="o", color="r", s=100, alpha=0.5)
    ax.text(ppm, height + 5e5, n + 1, ha="center", va="center")

ax.hlines(threshold, *uc.ppm_limits(), linestyle="--", color="k")
ax.set_ylim(top=1.2e7)
ax.set_xlim(200, 0)

plt.show()
