import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.offsetbox as offsetbox

# Load data
data = np.loadtxt("results/energies/energy_values.csv", delimiter=",", skiprows=1)
vs = data[:, 0]
m0_15 = data[:, 1]
m0_05 = data[:, 2]
m0_ADE = data[:, 3]
m0_ADE_a1 = data[:, 4]
m0_05_a16 = data[:, 5]

NAMES = [r"$l_0$ = $1.5$ nm", r"$l_0$ = $0.5$ nm", r"ADE ($\alpha=4$)", r"ADE ($\alpha=1$)", r"$l_0$ = $0.5$ nm ($\alpha=16$)"]

fig, ax = plt.subplots(1, 1, figsize=(6,4))

ax.scatter(vs, m0_15, s=60, color="blue", zorder=3)
ax.scatter(vs, m0_05, s=60, color="black", zorder=2)

# cubic spline fit 
f = interp1d(vs, m0_15, kind="cubic")
xnew = np.linspace(0.76, 0.97, num=100, endpoint=True)
ax.plot(xnew, f(xnew), color="blue", zorder=1, lw=2.0)

f = interp1d(vs, m0_05, kind="cubic")
xnew = np.linspace(0.76, 0.97, num=100, endpoint=True)
ax.plot(xnew, f(xnew), color="black", zorder=1, lw=2.0)

# ADE_a1
f = interp1d(vs, m0_ADE_a1, kind="cubic")
xnew = np.linspace(0.76, 0.97, num=100, endpoint=True)
ax.plot(xnew, f(xnew), color="red", zorder=0, lw=2.0, linestyle="--")


# add text near each plot, tilt the text
#ax.text(0.92, 2.03, NAMES[0], fontsize=14, color="blue", rotation=56)
#ax.text(0.905, 1.752, NAMES[1], fontsize=14, color="black", rotation=30)
#ax.text(0.88, 1.225, NAMES[2], fontsize=14, color="gray", alpha=0.25, rotation=3)
#ax.text(0.825, 1.61, NAMES[3], fontsize=14, color="red", rotation=9)

# load gifs
ab = offsetbox.AnnotationBbox(offsetbox.OffsetImage(plt.imread("results/energies/prolate.gif"), zoom=0.045), (0.94, 1.5), frameon=False, zorder=-1)
ax.add_artist(ab)
ab = offsetbox.AnnotationBbox(offsetbox.OffsetImage(plt.imread("results/energies/bud.gif"), zoom=0.045), (0.81, 2.2), frameon=False, zorder=-1)
ax.add_artist(ab)
ab = offsetbox.AnnotationBbox(offsetbox.OffsetImage(plt.imread("results/energies/dumbbell.gif"), zoom=0.045), (0.8, 1.33), frameon=False, zorder=-1)
ax.add_artist(ab)

ax.set_xlabel(r"$\nu$", fontsize=20)
ax.set_ylabel(r"$\frac{m_0}{4 \pi}$", fontsize=20)
ax.set_xlim(0.76, 0.97)
ax.set_ylim(1.0, 2.6)
ax.tick_params(axis="both", labelsize=15)
[ax.spines[spine].set_linewidth(2.5) for spine in ax.spines]
[ax.spines[spine].set_zorder(10) for spine in ax.spines]
[ax.tick_params(axis=axis, width=2.5) for axis in ["x", "y"]]
plt.savefig("results/energies/phase_diagram.svg", bbox_inches="tight", dpi=1000)
plt.savefig("results/energies/phase_diagram.png", bbox_inches="tight", dpi=300)
plt.clf()


fig, ax = plt.subplots(1, 1, figsize=(6,3))

ax.scatter(vs, m0_05_a16, s=60, color="green", zorder=2)

# ADE
f = interp1d(vs, m0_ADE, kind="cubic")
xnew = np.linspace(0.76, 0.97, num=100, endpoint=True)
ax.plot(xnew, f(xnew), color="gray", zorder=1, lw=2.0, alpha=0.25, linestyle="--")

# m0_05_a16
f = interp1d(vs, m0_05_a16, kind="cubic")
xnew = np.linspace(0.76, 0.97, num=100, endpoint=True)
ax.plot(xnew, f(xnew), color="green", zorder=0, lw=2.0)

# add text near each plot, tilt the text
#ax.text(0.92, 2.03, NAMES[0], fontsize=14, color="blue", rotation=56)
#ax.text(0.905, 1.752, NAMES[1], fontsize=14, color="black", rotation=30)
ax.text(0.88, 1.23, NAMES[2], fontsize=14, color="gray", alpha=0.5, rotation=3)
ax.text(0.825, 1.33, NAMES[4], fontsize=14, color="green", rotation=0)
#ax.text(0.825, 1.61, NAMES[3], fontsize=14, color="red", rotation=9)

ax.set_xlabel(r"$\nu$", fontsize=20)
ax.set_ylabel(r"$\frac{m_0}{4 \pi}$", fontsize=20)
ax.set_xlim(0.76, 0.97)
ax.set_ylim(1.0, 1.5)
ax.tick_params(axis="both", labelsize=15)
[ax.spines[spine].set_linewidth(2.5) for spine in ax.spines]
[ax.spines[spine].set_zorder(10) for spine in ax.spines]
[ax.tick_params(axis=axis, width=2.5) for axis in ["x", "y"]]
plt.savefig("results/energies/phase_diagram_si.pdf", bbox_inches="tight", dpi=1000)
plt.savefig("results/energies/phase_diagram_si.png", bbox_inches="tight", dpi=300)
plt.clf()
