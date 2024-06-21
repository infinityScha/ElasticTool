using PyCall
using DataFrames
import Plots:Animation, buildanimation  

@pyimport matplotlib
@pyimport matplotlib.pyplot as plt
@pyimport numpy as np
@pyimport mpl_toolkits.axes_grid1 as axes_grid1
@pyimport matplotlib.collections as mcoll
LineCollection = mcoll.LineCollection
make_axes_locatable = axes_grid1.make_axes_locatable
plt.rcParams["svg.fonttype"] = "none"
matplotlib.use("Agg")

counter_image = 0

function plot_coords(vol_patches::Array{VolumePatch,1})
    global counter_image
    counter_image += 1
    # dont plot the water patch 
    vol_patches_not_water = []
    for vol_patch in vol_patches
        if length(vol_patch.vol_elems[1].coords) != 2
            push!(vol_patches_not_water, vol_patch)
        end
    end

    # get virdis colormap
    autumn_cmap = plt.cm.get_cmap("viridis")
    norm = matplotlib.colors.Normalize(Φ_MIN, Φ_MAX)
    new_cmap = plt.cm.ScalarMappable(cmap=autumn_cmap, norm=norm)

    fig, ax = plt.subplots()
    ANY_Φ = false
    for vol_patch in vol_patches_not_water
        Φ₀ = vol_patch.Φ₀
        for vol_elem in vol_patch.vol_elems
            x = [coord_.r for coord_ in vol_elem.coords[1:4]]
            y = [coord_.z for coord_ in vol_elem.coords[1:4]]
            field = vol_elem.vol_patch.use_Φ ? mean([Φ(coord_.Φ₊) for coord_ in vol_elem.coords[1:2]]) : 0.0
            ax.fill(x, y, color=new_cmap.to_rgba(field))
            if vol_elem.vol_patch.use_Φ
                ANY_Φ = true
                # write the value of the field in the middle of the element
                ax.text(mean(x), mean(y), string(round(field*100, digits=2)), horizontalalignment="center", verticalalignment="center", fontsize=5)
            end
            push!(x, x[1])
            push!(y, y[1])
            ax.plot(x, y, "indigo", linewidth=4)
            ax.plot(x[1:2], y[1:2], "indigo", linewidth=8)
            ax.plot(x[3:4], y[3:4], "indigo", linewidth=8)
        end
    end
    plt.xlim([PLOTTING_MIN_R, PLOTTING_MAX_R])
    plt.ylim([PLOTTING_MIN_Z, PLOTTING_MAX_Z])
    ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    ax.set_aspect("equal")
    for spine in ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.xlabel("r [nm]", fontsize=150)
    plt.ylabel("z [nm]", fontsize=150)
    zoom = 20
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)
    if ANY_Φ
        cbar = plt.colorbar(new_cmap, cax=cax)
        cbar.ax.tick_params(labelsize=120, width=12, length=50, pad=20)
        cbar.ax.set_ylabel("Φ", fontsize=120, labelpad=20)
    end
    plt.savefig(FIG_DIR * "/conf" * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=1.0)
    plt.close("all")
end

function plot_energy_density(vol_patches::Array{VolumePatch,1}, which::String)
    global counter_image
    counter_image += 1
    which_dict = Dict("total" => [ENERGY_MAP["bending"], ENERGY_MAP["tilting"], ENERGY_MAP["stretching"], ENERGY_MAP["Φing"], ENERGY_MAP["selfing"]], "bending" => [ENERGY_MAP["bending"]], "tilting" => [ENERGY_MAP["tilting"]], "stretching" => [ENERGY_MAP["stretching"]], "Φing" => [ENERGY_MAP["Φing"]], "selfing" => [ENERGY_MAP["selfing"]])

    # dont plot the water patch 
    vol_patches_not_water = []
    for vol_patch in vol_patches
        if length(vol_patch.vol_elems[1].coords) != 2
            push!(vol_patches_not_water, vol_patch)
        end
    end
    en_arr = Float64[]
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            en = sum(vol_elem.energy[which_dict[which]])
            push!(en_arr, en)
        end
    end
    blim = minimum(en_arr)
    ulim = maximum(en_arr)
    autumn_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_colormap", ["thistle", "rebeccapurple"])
    norm = matplotlib.colors.Normalize(blim, ulim)
    new_cmap = plt.cm.ScalarMappable(cmap=autumn_cmap, norm=norm)

    fig, ax = plt.subplots()
    count = 0
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            en = sum(vol_elem.energy[which_dict[which]])
            x = [coord_.r for coord_ in vol_elem.coords[1:4]]
            y = [coord_.z for coord_ in vol_elem.coords[1:4]]
            ax.fill(x, y, color=new_cmap.to_rgba(en))
            ax.plot(x, y, "indigo", linewidth=4)
            ax.plot(x[1:2], y[1:2], "indigo", linewidth=8)
            ax.plot(x[3:4], y[3:4], "indigo", linewidth=8)
            count += 1
        end
    end
    plt.xlim([PLOTTING_MIN_R, PLOTTING_MAX_R])
    plt.ylim([PLOTTING_MIN_Z, PLOTTING_MAX_Z])
    ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    ax.set_aspect("equal")
    for spine in ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.xlabel("r [nm]", fontsize=150)
    plt.ylabel("z [nm]", fontsize=150)
    zoom = 20
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)
    cbar = plt.colorbar(new_cmap, cax=cax)
    cbar.ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    cbar.ax.set_ylabel("Energy density [pN/nm]", fontsize=120, labelpad=20)
    for spine in cbar.ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.savefig(FIG_DIR * "/en_dens_" * which * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=1.0)
    plt.close("all")
end

function plot_thickness(vol_patches::Array{VolumePatch,1})
    global counter_image
    counter_image += 1
    # dont plot the water patch 
    vol_patches_not_water = []
    for vol_patch in vol_patches
        if length(vol_patch.vol_elems[1].coords) != 2
            push!(vol_patches_not_water, vol_patch)
        end
    end
    thickness_arr = Float64[]
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            r1, r2, r3 , r4 = [coord_.r for coord_ in vol_elem.coords[1:4]]
            z1, z2, z3 , z4 = [coord_.z for coord_ in vol_elem.coords[1:4]]
            thickness1 = sqrt((r1-r4)^2 + (z1-z4)^2)
            thickness2 = sqrt((r2-r3)^2 + (z2-z3)^2)
            push!(thickness_arr, (thickness1+thickness2)/2)
        end
    end
    blim = minimum(thickness_arr)
    ulim = maximum(thickness_arr)
    autumn_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_colormap", ["thistle", "rebeccapurple"]) #plt.cm.get_cmap("autumn")
    norm = matplotlib.colors.Normalize(blim, ulim)
    new_cmap = plt.cm.ScalarMappable(cmap=autumn_cmap, norm=norm)

    fig, ax = plt.subplots()
    count = 0
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            x = [coord_.r for coord_ in vol_elem.coords[1:4]]
            y = [coord_.z for coord_ in vol_elem.coords[1:4]]
            thickness = 0.5*(sqrt((x[1]-x[4])^2 + (y[1]-y[4])^2) + sqrt((x[2]-x[3])^2 + (y[2]-y[3])^2))
            ax.fill(x, y, color=new_cmap.to_rgba(thickness))
            ax.plot(x, y, "indigo", linewidth=4)
            ax.plot(x[1:2], y[1:2], "indigo", linewidth=8)
            ax.plot(x[3:4], y[3:4], "indigo", linewidth=8)
            # write the value of the field in the middle of the element
            ax.text(mean(x), mean(y), string(round(thickness, digits=3)), horizontalalignment="center", verticalalignment="center", fontsize=5)
            count += 1
        end
    end
    plt.xlim([PLOTTING_MIN_R, PLOTTING_MAX_R])
    plt.ylim([PLOTTING_MIN_Z, PLOTTING_MAX_Z])
    ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    ax.set_aspect("equal")
    for spine in ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.xlabel("r [nm]", fontsize=150)
    plt.ylabel("z [nm]", fontsize=150)
    zoom = 20
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)
    cbar = plt.colorbar(new_cmap, cax=cax)
    cbar.ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    cbar.ax.set_ylabel("Length [nm]", fontsize=120, labelpad=20)
    for spine in cbar.ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.savefig(FIG_DIR * "/length_" * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=1.0)
    plt.close("all")
end

function plot_energy(vol_patches::Array{VolumePatch,1}, which::String)
    global counter_image
    counter_image += 1
    which_dict = Dict("total" => [2, 3, 4, 5, 7], "bending" => [2], "g_curving" => [3], "tilting" => [4], "stretching" => [5], "Φing" => [6], "selfing" => [7], "dsing" => [8])
    which_dict = Dict("total" => [ENERGY_MAP["bending"], ENERGY_MAP["tilting"], ENERGY_MAP["stretching"], ENERGY_MAP["Φing"], ENERGY_MAP["selfing"]], "bending" => [ENERGY_MAP["bending"]], "g_curving" => [ENERGY_MAP["g_curving"]], "tilting" => [ENERGY_MAP["tilting"]], "stretching" => [ENERGY_MAP["stretching"]], "Φing" => [ENERGY_MAP["Φing"]], "selfing" => [ENERGY_MAP["selfing"]], "dsing" => [ENERGY_MAP["dsing"]])
    # dont plot the water patch 
    vol_patches_not_water = []
    for vol_patch in vol_patches
        if length(vol_patch.vol_elems[1].coords) != 2
            push!(vol_patches_not_water, vol_patch)
        end
    end
    en_arr = Float64[]
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            en = sum(vol_elem.energy[which_dict[which]]) * vol_elem.energy[ENERGY_MAP["dS"]]
            push!(en_arr, en)
        end
    end
    blim = minimum(en_arr)
    ulim = maximum(en_arr)
    autumn_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_colormap", ["thistle", "rebeccapurple"]) #plt.cm.get_cmap("autumn")
    norm = matplotlib.colors.Normalize(blim, ulim)
    new_cmap = plt.cm.ScalarMappable(cmap=autumn_cmap, norm=norm)

    fig, ax = plt.subplots()
    count = 0
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            en = sum(vol_elem.energy[which_dict[which]]) * vol_elem.energy[ENERGY_MAP["dS"]]
            x = [coord_.r for coord_ in vol_elem.coords[1:4]]
            y = [coord_.z for coord_ in vol_elem.coords[1:4]]
            ax.fill(x, y, color=new_cmap.to_rgba(en))
            ax.plot(x, y, "indigo", linewidth=4)
            ax.plot(x[1:2], y[1:2], "indigo", linewidth=8)
            ax.plot(x[3:4], y[3:4], "indigo", linewidth=8)
            count += 1
        end
    end
    plt.xlim([PLOTTING_MIN_R, PLOTTING_MAX_R])
    plt.ylim([PLOTTING_MIN_Z, PLOTTING_MAX_Z])
    ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    ax.set_aspect("equal")
    for spine in ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.xlabel("r [nm]", fontsize=150)
    plt.ylabel("z [nm]", fontsize=150)
    zoom = 20
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * zoom, h * zoom)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)
    cbar = plt.colorbar(new_cmap, cax=cax)
    cbar.ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    cbar.ax.set_ylabel("Energy [pN nm]", fontsize=120, labelpad=20)
    for spine in cbar.ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.savefig(FIG_DIR * "/en_" * which * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=1.0)
    plt.close("all")
end

function plot_energy_dsing(vol_patches::Array{VolumePatch,1})
    global counter_image
    counter_image += 1
    # dont plot the water patch 
    vol_patches_not_water = []
    for vol_patch in vol_patches
        if length(vol_patch.vol_elems[1].coords) != 2
            push!(vol_patches_not_water, vol_patch)
        end
    end
    en_arr = Float64[]
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            push!(en_arr, vol_elem.energy[ENERGY_MAP["dsing"]])
        end
    end
    blim = minimum(en_arr)
    ulim = maximum(en_arr)
    autumn_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_colormap", ["thistle", "rebeccapurple"])
    norm = matplotlib.colors.Normalize(blim, ulim)
    new_cmap = plt.cm.ScalarMappable(cmap=autumn_cmap, norm=norm)

    fig, ax = plt.subplots()
    count = 0
    for vol_patch in vol_patches_not_water
        for vol_elem in vol_patch.vol_elems
            en = vol_elem.energy[ENERGY_MAP["dsing"]]
            x = [coord_.r for coord_ in vol_elem.coords[1:4]]
            y = [coord_.z for coord_ in vol_elem.coords[1:4]]
            ax.fill(x, y, color=new_cmap.to_rgba(en))
            ax.plot(x, y, "indigo", linewidth=4)
            ax.plot(x[1:2], y[1:2], "indigo", linewidth=8)
            ax.plot(x[3:4], y[3:4], "indigo", linewidth=8)
            count += 1
        end
    end
    ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    ax.set_aspect("equal")
    for spine in ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.xlabel("r [nm]", fontsize=150)
    plt.ylabel("z [nm]", fontsize=150)

    w, h = fig.get_size_inches()
    zoom = 20
    fig.set_size_inches(w * zoom, h * zoom)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)
    cbar = plt.colorbar(new_cmap, cax=cax)
    cbar.ax.tick_params(labelsize=120, width=12, length=50, pad=20)
    cbar.ax.set_ylabel("Energy [pN nm]", fontsize=120, labelpad=20)
    for spine in cbar.ax.spines.values()
        spine.set_linewidth(12)
    end
    plt.savefig(FIG_DIR * "/dsing" * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=1.0)
    plt.close("all")
end

function plot_reaction(rem_patches_arr::Vector{Vector{RemeshingPatch}})
    global counter_image
    counter_image += 1

    norm = matplotlib.colors.Normalize(Φ_MIN, Φ_MAX)

    images = String[]
    for (i, rem_patches) in enumerate(rem_patches_arr)
        fig, ax = plt.subplots(figsize=(8, 8*(PLOTTING_MAX_Z-PLOTTING_MIN_Z)/(2*PLOTTING_MAX_R)))
        for rem_patch in rem_patches
            ax.plot([coord.r for coord in rem_patch.mid_coords], [coord.z for coord in rem_patch.mid_coords], "grey", linestyle="dotted")
            ax.plot([-coord.r for coord in rem_patch.mid_coords], [coord.z for coord in rem_patch.mid_coords], "grey", linestyle="dotted")
            for (coords, vol_patch) in zip([rem_patch.upper_coords, rem_patch.lower_coords], [rem_patch.upper_vol_patch, rem_patch.lower_vol_patch])
                if coords ≠ nothing
                    r = [coord.r for coord in coords]
                    z = [coord.z for coord in coords]
                    Φ_ = [vol_patch.use_Φ ? Φ(coord.Φ₊) : 0.0 for coord in coords]
                    for sign_ in [-1.0,1.0]
                        points = np.array([sign_*r, z])'
                        points = np.reshape(points, (-1, 1, 2))
                        segments = np.concatenate([points[1:end-1, :, :], points[2:end, :, :]], axis=1)
                        lc = LineCollection(segments, cmap="viridis", norm=norm)
                        lc.set_array(Φ_)
                        ax.add_collection(lc)
                    end
                end
            end
        end
        ax.set_xlim([-PLOTTING_MAX_R, PLOTTING_MAX_R])
        ax.set_ylim([PLOTTING_MIN_Z, PLOTTING_MAX_Z]) 
        ax.set_aspect("equal")
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.spines["bottom"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        temp_file = @sprintf("%06d.png", i)
        plt.savefig(FIG_DIR * "/" * temp_file, dpi = 300, transparent = true, facecolor = "white", bbox_inches="tight")
        plt.close("all")
        push!(images, temp_file)
    end
    anim = Animation(FIG_DIR,images)
    buildanimation(anim, FIG_DIR * "/path" * string(counter_image) * ".gif", fps=5, show_msg=false, loop=0) 
    # delete the temporary files
    for temp_file in images
        rm(FIG_DIR * "/" * temp_file)
    end
    rm(FIG_DIR * "/" * "palette.bmp")
end

function plot_energy_path(energy_arr::Vector{Float64}, coords_arr, indep_coords_arr, neb_coords_placement_in_grad)
    global counter_image
    counter_image += 1
    # shift the energy array to have zero energy at the beginning
    energy_arr = copy(energy_arr) .- energy_arr[1]
    # get the forces on the NEB coordinates
    temp = [gradient_parallel(indep_coords) for indep_coords in indep_coords_arr]
    grad_arr = [temp_[neb_coords_placement_in_grad_] for (temp_, neb_coords_placement_in_grad_) in zip(temp, neb_coords_placement_in_grad)]
    # keep only the part of the gradient which is tangential to the path and get the coefficients for the cubic spline interpolation
    pos_arr = [pos(coords)[neb_coords_placement_in_grad_] for (coords, neb_coords_placement_in_grad_) in zip(coords_arr, neb_coords_placement_in_grad)]
    r, a, b, c, d = Float64[], Float64[], Float64[], Float64[], Float64[]
    for (pos_i, pos_i_1, grad_i, grad_i_1, energy_i, energy_i_1) in zip(pos_arr[1:end-1], pos_arr[2:end], grad_arr[1:end-1], grad_arr[2:end], energy_arr[1:end-1], energy_arr[2:end])
        τ = pos_i_1 .- pos_i
        R = norm(τ)
        τ = τ ./ R
        force_i, force_i_1 = -dot(grad_i, τ), -dot(grad_i_1, τ)
        push!(r, R)
        push!(a, -2.0*(energy_i_1-energy_i)/R/R/R-(force_i+force_i_1)/R/R)
        push!(b, 3.0*(energy_i_1-energy_i)/R/R+(2.0*force_i+force_i_1)/R)
        push!(c, -force_i)
        push!(d, energy_i)
    end

    pushfirst!(r, 0.0)
    r = cumsum(r)

    function cubic_interpolation(x)
        i = findlast(r .< x)
        x -= r[i]
        return a[i]*x*x*x + b[i]*x*x + c[i]*x + d[i]
    end

    # plot the energy with the cubic interpolation
    fig, ax = plt.subplots()
    x = range(0.001, stop=r[end]-0.001, length=1000)
    y = cubic_interpolation.(x)./4.114
    energy_arr /= 4.114
    ax.plot(x./r[end], y, color="black")
    ax.set_xlim(0.0, 1.0)
    
    # plot the energy as scatter
    ax.plot(r./r[end], energy_arr, "o", color="black")
    ax.set_xlabel("ξ", fontsize=20)
    ax.set_ylabel("F [kT]", fontsize=20)
    ax.tick_params(labelsize=20)
    ax.set_title("Energy path", fontsize=20)
    #add subtitle with the barrier energy
    ax.text(0.5, 0.02, "Barrier: "*string(round(maximum(y)-y[1], digits=2))*" kT", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes, fontsize=10)
    plt.savefig(FIG_DIR * "/energy_path" * string(counter_image) * ".svg", bbox_inches="tight", pad_inches=0.5)
    plt.close("all")
    # Print the initial energy, the final energy the maximum energy and the barrier up to 2 decimals
    println("Initial energy: ", round(y[1], digits=2), " kT Final energy: ", round(y[end], digits=2), " kT Maximum energy: ", round(maximum(y), digits=2), " kT Barrier: ", round(maximum(y)-y[1], digits=2), " kT")
    # save the information to a csv
    CSV.write(FIG_DIR * "/energy_path" * string(counter_image) * ".csv", DataFrame(ξ=r./r[end], F=energy_arr))
    # also the interpolated energy
    CSV.write(FIG_DIR * "/energy_path_interpolated" * string(counter_image) * ".csv", DataFrame(ξ=x./r[end], F=y))
end