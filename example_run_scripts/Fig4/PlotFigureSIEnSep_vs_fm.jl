using PyCall
using CSV
using DataFrames

@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.image as mpimg
@pyimport matplotlib.offsetbox as offsetbox

# plots energy separation, taken from the mixed system

COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "black"]

fig = plt.figure(figsize=(8, 5))
ax = plt.gca()

SYSTEMS = ["5", "mixed", "15", "20", "40", "upper"]
NAMES = [raw"$\phi_0 = 0.05, \tilde{J}_s=0.8$ nm⁻¹", raw"$\phi_0 = 0.1, \tilde{J}_s=0.4$ nm⁻¹", raw"$\phi_0 = 0.15, \tilde{J}_s=0.267$ nm⁻¹", raw"$\phi_0 = 0.2, \tilde{J}_s=0.2$ nm⁻¹", raw"$\phi_0 = 0.4, \tilde{J}_s=0.1$ nm⁻¹", raw"no mix, $J^u_s=0.04$ nm⁻¹"]

for (SYSTEM, NAME, COLOR) in zip(SYSTEMS, NAMES, COLORS)
    data = CSV.read("results/$SYSTEM/path/for_en_sep.csv", DataFrame)

    ξ_arr = LinRange(0.0, 1.0, 10)


    f_el_arr = (data[!, "Bending"] .+ data[!, "Tilting"] .+ data[!, "Stretching"])./4.114
    f_el_arr = f_el_arr .- f_el_arr[1]

    if SYSTEM != "upper"
        f_mix_arr = data[!, "Φing"]./4.114
        f_mix_arr = f_mix_arr .- f_mix_arr[1]    
    end

    plt.plot(ξ_arr, f_el_arr, "-", color=COLOR, linewidth=2.0, alpha=1.0, zorder=0)

    # plot the energy as scatter
    plt.scatter(ξ_arr, f_el_arr, s=60, color=COLOR, label=NAME, alpha=1.0, zorder=10)
    if SYSTEM != "upper"
        plt.scatter(ξ_arr, f_mix_arr, s=60, color=COLOR, alpha=1.0, zorder=10)
        plt.plot(ξ_arr, f_mix_arr, "--", color=COLOR, linewidth=2.0, alpha=1.0, zorder=0)
    end
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
plt.legend(fontsize=12, loc="lower left", ncol=1, frameon=false, borderaxespad=0.0, handletextpad=0.0, columnspacing=0.0, title_fontsize=13)
plt.savefig("results/figS9.pdf", bbox_inches="tight")
plt.savefig("results/figS9.png", bbox_inches="tight", dpi=300)
plt.close("all")
