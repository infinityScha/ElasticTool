using PyCall
using CSV
using DataFrames

@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.image as mpimg
@pyimport matplotlib.offsetbox as offsetbox

# plots energy separation, taken from the mixed system

SYSTEMS = ["Bending", "Tilting", "Stretching", "Φing"]
NAMES = ["Bending", "Tilting", "Stretching", raw"$f^{\mathrm{mix}}$"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

fig = plt.figure(figsize=(8, 5))
ax = plt.gca()

data = CSV.read("for_en_sep.csv", DataFrame)

ξ_arr = LinRange(0.0, 1.0, 10)

for (SYSTEM, NAME, COLOR) in zip(SYSTEMS, NAMES, COLORS)
    property_arr = (data[!, SYSTEM] .- data[1, SYSTEM])/4.114
    ax.plot(ξ_arr, property_arr, "-", color=COLOR, linewidth=2.0, alpha=1.0, zorder=0)
    
    # plot the energy as scatter
    plt.scatter(ξ_arr, property_arr, s=60, color=COLOR, label=NAME, alpha=1.0, zorder=10)
end

ax.set_xlim(0.0, 1.0)
[i.set_linewidth(1.5) for i in ax.spines.values()]
# make axis load last
[spine.set_zorder(10) for spine in ax.spines.values()]
ax.xaxis.set_tick_params(width=1.5)
ax.yaxis.set_tick_params(width=1.5)
ax.set_xlabel("ξ", fontsize=20)
ax.set_ylabel(raw"ΔF [k$_\mathrm{B}$T]", fontsize=20)
ax.tick_params(labelsize=20)
plt.legend(fontsize=13, loc="lower left", ncol=1, frameon=false, borderaxespad=0.0, handletextpad=0.0, columnspacing=0.0, title_fontsize=13)
plt.savefig("results/figS8.pdf", bbox_inches="tight")
plt.savefig("results/figS8.png", bbox_inches="tight", dpi=300)
plt.close("all")