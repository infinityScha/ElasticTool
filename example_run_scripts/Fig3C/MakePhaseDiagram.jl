const MAIN_PATH = "phase_diagram/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 120.0
const PLOTTING_MIN_Z = -70.0
const Φ_MAX = 0.2
const Φ_MIN = 0.0

# is the system unstable? then setting it to yes will improve stability
const IS_UNSTABLE = false # false or true
# is the system PLANAR?
const IS_PLANAR = false
# is the system consists of domains (i.e. Lo/Ld)?
const DOMAINS = false # false or true
const DOMAINS_WIDTH = 4.0 #  in nm
const DOMAINS_ENF_STR = 1.0e3 # in pN
const ADD_LINE_TENSION = false # false or true
const LINE_TENSION = 1.0 # in pN
# interaction of membrane in the Z-axis
const IS_Z_SELF_INTERACTION = false # false or true
const ϵ_Z_SELF_INTERACTION = 0.0 #39 * 1.65954 * 0.5 # = 32.36, estimated energy (kJ/mol) per R9 * conversion to pN * nm * R9 surface density
const σ_Z_SELF_INTERACTION = 1.5 # in nm
# mixing interaction paramter
const CHI = 0.0
# smooth line parameter
const Kₗᵢₙₑ = 0.2 * 0.639 # in kT / nm, this gives a penalty of 0.2 kT / nm for a naive linear interface with width of 0.8 nm
# for energy report
const ENERGY_FACTOR = 1.0

include("../../ElasticTool.jl")

result_names = ["results/", "results_reversed/", "results_domain2/"]

JSs = [-0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
CHIs = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

function get_final_energy(result_name, js, chi)
    if chi == 0.0
        chi = 0
    end
    if js == 0.0
        js = 0
    end
    energy_path = result_name * "$(js)/$(chi)/energy.csv"

    # check if the file exists
    if !isfile(energy_path)
        return 100000.0
    end

    energies_file = readlines(energy_path)
    energy = split(energies_file[end], ",")[2]
    if energy == "energy"
        return 100000.0
    end
    return parse(Float64, energy)
end

energy_dict = Dict()
for result_name in result_names
    energy_dict[result_name] = Dict()
    for js in JSs
        energy_dict[result_name][js] = Dict()
        for chi in CHIs
            energy_dict[result_name][js][chi] = get_final_energy(result_name, js, chi)
        end
    end
end

# make array with minimal energy for each js and chi of the possible results
min_energy = zeros(length(JSs), length(CHIs))
which_has_min_energy = zeros(length(JSs), length(CHIs))
for (i, js) in enumerate(JSs)
    for (j, chi) in enumerate(CHIs)
        min_energy[i, j] = minimum([energy_dict[result_name][js][chi] for result_name in result_names])
        which_has_min_energy[i, j] = argmin([energy_dict[result_name][js][chi] for result_name in result_names])
    end
end

function get_enrichment(which, js, chi)
    if chi == 0.0
        chi = 0
    end
    if js == 0.0
        js = 0
    end
    # load the final state
    coords, indep_coords, vol_patches, rem_patches = load_state(result_names[which]*"$(js)/$(chi)/states/final.coords", result_names[which]*"$(js)/$(chi)/states/final.vol_patches", result_names[which]*"$(js)/$(chi)/states/final.rem_patches")
    # get the enrichment
    rem_patch = rem_patches[1]
    upper_rs = [coord.r for coord in rem_patch.upper_coords]
    upper_zs = [coord.z for coord in rem_patch.upper_coords]
    upper_Φ₊s = [coord.Φ₊ for coord in rem_patch.upper_coords]
    upper_Φs = [Φ(Φ₊) for Φ₊ in upper_Φ₊s]
    # get curve length representation
    ds = [sqrt((r0 - r1)^2 + (z0 - z1)^2) for (r0, r1, z0, z1) in zip(upper_rs[1:end-1], upper_rs[2:end], upper_zs[1:end-1], upper_zs[2:end])]
    dA = [0.5*(r0+r1)*ds_ for (r0, r1, ds_) in zip(upper_rs[1:end-1], upper_rs[2:end], ds)]
    dn = [0.5*(Φ0+Φ1)*dA_ for (Φ0, Φ1, dA_) in zip(upper_Φs[1:end-1], upper_Φs[2:end], dA)]
    pushfirst!(ds, 0.0)
    pushfirst!(dA, 0.0)
    pushfirst!(dn, 0.0)
    s = cumsum(ds)
    A_ = cumsum(dA)
    n = cumsum(dn)
    # find the last maximum of r
    k = length(upper_rs) - 1
    while upper_rs[k] < upper_rs[k-1]
        k -= 1
    end
    # find the next minimum
    while upper_rs[k] > upper_rs[k-1]
        k -= 1
    end
    # get the total area up to the first minimum
    A₀ = A_[k]
    # get the total n up to the first minimum
    n₀ = n[k]
    Φ₀ = n₀/A₀
    # get the total area from the first minimum
    A₁ = A_[end] - A₀
    # get the total n from the first minimum
    n₁ = n[end] - n₀
    Φ₁ = n₁/A₁
    #plt.plot(s[1:k], upper_rs[1:k], color="black")
    #plt.plot(s[k:end], upper_rs[k:end], color="red")
    #plt.xlabel("s (nm)")
    #plt.ylabel("r (nm)")
    #plt.savefig("phase_diagram/test_r.png", bbox_inches="tight")
    #plt.clf()
    #plt.plot(s[1:k], upper_Φs[1:k], color="black")
    #plt.plot(s[k:end], upper_Φs[k:end], color="red")
    #plt.xlabel("s (nm)")
    #plt.ylabel("Φ")
    #plt.savefig("phase_diagram/test_Φ.png", bbox_inches="tight")
    #plt.clf()
    #println("js: $(js), chi: $(chi), Φ₀: $(Φ₀), Φ₁: $(Φ₁)")
    #println("A₀: $(A₀), A₁: $(A₁)")
    #println("n₀: $(n₀), n₁: $(n₁)")
    return Φ₀, Φ₁
end

function how_many_interfaces(which, js, chi)
    if chi == 0.0
        chi = 0
    end
    if js == 0.0
        js = 0
    end
    # load the final state
    coords, indep_coords, vol_patches, rem_patches = load_state(result_names[which]*"$(js)/$(chi)/states/final.coords", result_names[which]*"$(js)/$(chi)/states/final.vol_patches", result_names[which]*"$(js)/$(chi)/states/final.rem_patches")
    rem_patch = rem_patches[1]
    upper_Φ₊s = [coord.Φ₊ for coord in rem_patch.upper_coords]
    upper_Φs = [Φ(Φ₊) for Φ₊ in upper_Φ₊s]
    count_change = 0
    which_now = upper_Φs[1] > 0.5 ? 1 : 0
    for Φ in upper_Φs
        if abs(which_now - Φ) > 0.5
            count_change += 1
            which_now = 1 - which_now
        end
    end
    return count_change
end


#get_enrichment(1, 0.0, 3.0)
#get_enrichment(1, 0.2, 3.0)


Φ₀_matrix = zeros(length(JSs), length(CHIs))
Φ₁_matrix = zeros(length(JSs), length(CHIs))
JSs_matrix = zeros(length(JSs), length(CHIs))
CHIs_matrix = zeros(length(JSs), length(CHIs))
interfaces_matrix = zeros(length(JSs), length(CHIs))
for (i, js) in enumerate(JSs)
    for (j, chi) in enumerate(CHIs)
        #if js == 0.0 && chi == 3.0
        #    which = -1
        #    Φ₀_matrix[i, j] = 0.0909
        #    Φ₁_matrix[i, j] = 0.0909
        #    JSs_matrix[i, j] = js
        #    CHIs_matrix[i, j] = chi
        #    interfaces_matrix[i, j] = -1
        #    continue
        #end
        which = Int(which_has_min_energy[i, j])
        Φ₀, Φ₁ = get_enrichment(which, js, chi)
        Φ₀_matrix[i, j] = Φ₀
        Φ₁_matrix[i, j] = Φ₁
        JSs_matrix[i, j] = js
        CHIs_matrix[i, j] = chi
        interfaces_matrix[i, j] = how_many_interfaces(which, js, chi)
        if js < 0.05 && which == 3 && interfaces_matrix[i, j] == 1
            interfaces_matrix[i, j] = -2
        end
        #if js < 0.05 && chi == 3.0
        #    interfaces_matrix[i, j] = -2
        #    println(Φ₀, " ", Φ₁, " ", which)
        #end
        #if js == -0.2 && chi == 2.8
        #    interfaces_matrix[i, j] = -2
        #end
    end
end

# plot the enrichment heatmap, without smoothing
plt.figure(figsize=(8, 6))
#plt.imshow(Φ₀_matrix .- Φ₁_matrix, cmap="nipy_spectral", extent=[CHIs[1], CHIs[end], JSs[1], JSs[end]], aspect="auto", alpha=0.3)
plt.pcolormesh(CHIs_matrix, JSs_matrix, Φ₀_matrix .- Φ₁_matrix, cmap="nipy_spectral", shading="auto", alpha=0.5)
plt.colorbar(label="Enrichment")
plt.xlabel("χ")
plt.ylabel("Jₛ")
plt.savefig("phase_diagram/enrichment_unsmooth.png", bbox_inches="tight")
plt.close()



# interpolate 
@pyimport scipy.interpolate as spi
@pyimport numpy as np

x = range(-0.05, stop=0.3, length=100)
y = range(0.0, stop=3.0, length=100)
grid_x, grid_y = np.meshgrid(x, y)

Φ₀_matrix = spi.griddata((JSs_matrix[:], CHIs_matrix[:]), Φ₀_matrix[:], (grid_x, grid_y), method="cubic")
Φ₁_matrix = spi.griddata((JSs_matrix[:], CHIs_matrix[:]), Φ₁_matrix[:], (grid_x, grid_y), method="cubic")
which_has_min_energy = spi.griddata((JSs_matrix[:], CHIs_matrix[:]), which_has_min_energy[:], (grid_x, grid_y), method="nearest")
#interfaces_matrix = spi.griddata((JSs_matrix[:], CHIs_matrix[:]), interfaces_matrix[:], (grid_x, grid_y), method="nearest")

enrichment_matrix = Φ₀_matrix .- Φ₁_matrix

orig_JSs_matrix = JSs_matrix
orig_CHIs_matrix = CHIs_matrix

JSs_matrix = grid_x
CHIs_matrix = grid_y

# the matrices y-axis
Φ₀_matrix = Φ₀_matrix[end:-1:1, :]
Φ₁_matrix = Φ₁_matrix[end:-1:1, :]
enrichment_matrix = enrichment_matrix[end:-1:1, :]
JSs_matrix = JSs_matrix[end:-1:1, :]
CHIs_matrix = CHIs_matrix[end:-1:1, :]
#interfaces_matrix = interfaces_matrix[end:-1:1, :]
which_has_min_energy = which_has_min_energy[end:-1:1, :]

# plot the enrichment heatmap, use larger fontsize

# set large border widths
fig = plt.figure(figsize=(8, 5))
#plt.pcolormesh(CHIs_matrix, JSs_matrix, enrichment_matrix, cmap="nipy_spectral", shading="gouraud", alpha=0.5, zorder=0)
re_enrichment_matrix = np.flip(np.flip(enrichment_matrix, axis=1), axis=0)'
plt.imshow(re_enrichment_matrix, cmap="nipy_spectral", extent=[CHIs[1], CHIs[end], JSs[1], JSs[end]], aspect="auto", alpha=0.5, zorder=0)
ax = plt.gca()
[i.set_linewidth(2.5) for i in ax.spines.values()]
# make axis load last
[spine.set_zorder(10) for spine in ax.spines.values()]
ax.xaxis.set_tick_params(width=1.5)
ax.yaxis.set_tick_params(width=1.5)
# add a color bar based on the imshow and not on the contour
cbar = plt.colorbar()
[i.set_linewidth(2.0) for i in cbar.ax.spines.values()]
cbar.ax.tick_params(labelsize=15)
cbar.ax.xaxis.set_tick_params(width=1.5)
cbar.ax.yaxis.set_tick_params(width=1.5)
# set the label on the top of the color bar
cbar.ax.set_title(raw"$\Delta \Phi$", fontsize=20)
# add a countor line for the 0.05 enrichment, write it next to the line
plt.contour(CHIs_matrix, JSs_matrix, enrichment_matrix, levels=[0.065], colors="black", linewidths=3, linestyles="--", zorder=2)
plt.text(1.14, 0.17, "ΔΦ = 0.065", fontsize=18, color="black")
# add dots in different colors for which has the minimum energy
colors_dict = Dict(0 => "yellow", 1 => "blue", 2 => "magenta", 3 => "cyan", 4 => "red", -1 => "black", -2 => "lime")
for (js, chi, interfaces) in zip(orig_JSs_matrix[:], orig_CHIs_matrix[:], interfaces_matrix[:])
    plt.scatter([chi], [js], color=colors_dict[Int(interfaces)], s=60, edgecolors="black", linewidths=0.75, zorder=3)
end
plt.xlabel(raw"$\chi$", fontsize=20)
plt.ylabel(raw"$\tilde{J}_s$  [nm$^{-1}$]", fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig("phase_diagram/enrichment.pdf", bbox_inches="tight")
plt.savefig("phase_diagram/enrichment.png", bbox_inches="tight", dpi=300)
plt.close()

# plot ϕ₀ heatmap
plt.figure(figsize=(8, 6))
plt.imshow(Φ₀_matrix, cmap="viridis", extent=[CHIs[1], CHIs[end], JSs[1], JSs[end]], aspect="auto", alpha=0.3)
plt.colorbar(label="ϕ₀")
plt.xlabel("χ")
plt.ylabel("Jₛ")
plt.savefig("phase_diagram/Φ₀.png", bbox_inches="tight")
plt.close()

# plot ϕ₁ heatmap
plt.figure(figsize=(8, 6))
plt.imshow(Φ₁_matrix, cmap="viridis", extent=[CHIs[1], CHIs[end], JSs[1], JSs[end]], aspect="auto", alpha=0.3)
plt.colorbar(label="ϕ₁")
plt.xlabel("χ")
plt.ylabel("Jₛ")
plt.savefig("phase_diagram/Φ₁.png", bbox_inches="tight")
plt.close()

# plot the which heatmap
plt.figure(figsize=(8, 6))
plt.imshow(which_has_min_energy, cmap="viridis", extent=[CHIs[1], CHIs[end], JSs[1], JSs[end]], aspect="auto", alpha=0.3)
plt.colorbar(label="Which")
plt.xlabel("χ")
plt.ylabel("Jₛ")
plt.savefig("phase_diagram/which.png", bbox_inches="tight")
plt.close()
