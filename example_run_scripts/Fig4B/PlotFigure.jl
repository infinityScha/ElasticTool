using PyCall
using CSV
using DataFrames

@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.image as mpimg
@pyimport matplotlib.offsetbox as offsetbox

SYSTEMS = ["0.0", "0.4", "0.8", "1.2", "1.6", "2.0"]
NAMES = ["χ=0", "χ=0.4", "χ=0.8", "χ=1.2", "χ=1.6", "χ=2"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

fig = plt.figure(figsize=(8, 5))
ax = plt.gca()

i=0
for (SYSTEM, NAME, COLOR) in zip(SYSTEMS, NAMES, COLORS)
    global i
    # find the last energy path in the directory
    dir = "results/" * SYSTEM * "/figs/"
    files = filter(x -> occursin("energy_path", x), readdir(dir))
    files = filter(x -> !occursin("interp", x), files)
    # extract the number of the last image
    numbers = map(x -> parse(Int, split(x, "energy_path")[2][1:end-4]), files)
    counter_image = maximum(numbers)
    # load the data
    d = CSV.read("results/" * SYSTEM * "/figs/energy_path$(counter_image).csv", DataFrame)
    ξ_arr, energy_arr = d[:, 1], d[:, 2]
    d = CSV.read("results/" * SYSTEM * "/figs/energy_path_interpolated$(counter_image).csv", DataFrame)
    ξ_interp_arr, en_interp_arr = d[:, 1], d[:, 2]

    # plot the energy with the cubic interpolation
    ax.plot(ξ_interp_arr, en_interp_arr, "-", color=COLOR, linewidth=2.0, alpha=1.0, zorder=0)
    
    # plot the energy as scatter
    plt.scatter(ξ_arr, energy_arr, s=60, color=COLOR, label=NAME, alpha=1.0, zorder=10-i)
    i += 1
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
plt.legend(fontsize=13, loc="lower left", ncol=1, frameon=false, borderaxespad=0.0, handletextpad=0.0, columnspacing=0.0)
plt.savefig("results/fig4.pdf", bbox_inches="tight")
plt.savefig("results/fig4.png", bbox_inches="tight", dpi=300)
plt.close("all")