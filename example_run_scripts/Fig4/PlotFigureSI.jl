using PyCall
using CSV
using DataFrames

@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.image as mpimg
@pyimport matplotlib.offsetbox as offsetbox

SYSTEMS = ["mixed/5", "mixed", "mixed/15"]
NAMES = ["M=5", "M=10", "M=15"]
LINES = ["-", "--", ":"]

fig = plt.figure(figsize=(8, 5))
ax = plt.gca()


for (SYSTEM, NAME, LINE) in zip(SYSTEMS, NAMES, LINES)
    # find the last energy path in the directory
    dir = "results/" * SYSTEM * "/path/figs/"
    files = filter(x -> occursin("energy_path", x), readdir(dir))
    files = filter(x -> !occursin("interp", x), files)
    # extract the number of the last image
    numbers = map(x -> parse(Int, split(x, "energy_path")[2][1:end-4]), files)
    counter_image = maximum(numbers)
    # load the data
    d = CSV.read("results/" * SYSTEM * "/path/figs/energy_path$(counter_image).csv", DataFrame)
    ξ_arr, energy_arr = d[:, 1], d[:, 2]
    d = CSV.read("results/" * SYSTEM * "/path/figs/energy_path_interpolated$(counter_image).csv", DataFrame)
    ξ_interp_arr, en_interp_arr = d[:, 1], d[:, 2]

    # plot the energy with the cubic interpolation
    ax.plot(ξ_interp_arr, en_interp_arr, LINE, label=NAME, color="black", linewidth=2.0, alpha=1.0)
    
    # plot the energy as scatter
    #plt.scatter(ξ_arr, energy_arr, s=60, color="black", linewidths=0.75, edgecolors="black", alpha=1.0)
end

ax.set_xlim(0.0, 1.0)
[i.set_linewidth(1.5) for i in ax.spines.values()]
ax.xaxis.set_tick_params(width=1.5)
ax.yaxis.set_tick_params(width=1.5)
ax.set_xlabel("ξ", fontsize=20)
ax.set_ylabel(raw"ΔF [k$_\mathrm{b}$T]", fontsize=20)
ax.tick_params(labelsize=20)
plt.legend(fontsize=13, loc="lower left", ncol=1, frameon=false, borderaxespad=0.0, handletextpad=0.5, columnspacing=0.0)
plt.savefig("results/figS4.pdf", bbox_inches="tight")
plt.savefig("results/figS4.png", bbox_inches="tight", dpi=300)
plt.close("all")