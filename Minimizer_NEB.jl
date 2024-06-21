const K = 5.0e2

struct NEBCoordinatePlace
    rem_patch_idx::Int64
    at_start::Bool
end

function _ds_over_coords(rem_patches)
    ds_arr = Float64[]
    for rem_patch in rem_patches
        ds = [sqrt((coord0.r - coord1.r) * (coord0.r - coord1.r) + (coord0.z - coord1.z) * (coord0.z - coord1.z)) for (coord0, coord1) in zip(rem_patch.mid_coords[1:end-1], rem_patch.mid_coords[2:end])]
        ds = cumsum(ds)
        push!(ds_arr, ds[end])
    end
    return ds_arr
end


function ds_arr_over_coords(rem_patch)
    ds = [sqrt((coord0.r - coord1.r) * (coord0.r - coord1.r) + (coord0.z - coord1.z) * (coord0.z - coord1.z)) for (coord0, coord1) in zip(rem_patch.mid_coords[1:end-1], rem_patch.mid_coords[2:end])]
    ds = cumsum(ds)
    pushfirst!(ds, 0.0)
    return ds
end

function _get_placement_in_gradient(coords_arr, rem_patches_arr, neb_coords_place::Vector{NEBCoordinatePlace})
    neb_indices_arr = Vector{Int64}[]
    for (coords, rem_patches) in zip(coords_arr, rem_patches_arr)
        neb_coords_indices = [rem_patches[place.rem_patch_idx].upper_coords[place.at_start ? 1 : end].idx for place in neb_coords_place]
        
        temp = Dict{Int64, Vector{Int64}}()
        counter_ = 1
        for coord in coords
            to_add = coord.idx in neb_coords_indices
            if to_add
                idx = findfirst(x -> x == coord.idx, neb_coords_indices)
                temp[idx] = Int64[]
            end
            if !coord.fixed_z
                if to_add
                    push!(temp[idx], counter_)
                end
                counter_ += 1
            end
            if !coord.fixed_r
                if to_add
                    push!(temp[idx], counter_)
                end
                counter_ += 1
            end
            if !coord.fixed_Φ₊
                counter_ += 1
            end
        end
        # sort the neb indices array
        temp2 = vcat([temp[i] for i in 1:length(neb_coords_indices)]...)
        push!(neb_indices_arr, temp2)
    end
    return neb_indices_arr
end

function choose_climbing_image(vol_patches_arr)
    en_arr = [energy(vol_patches, whole=true) for vol_patches in vol_patches_arr]
    return argmax(en_arr)
end

# TODO: implement energy factor!

function neb_method(coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr, neb_coords_place::Vector{NEBCoordinatePlace}, ; μᵥ::Float64=floatmax(), μᵩ::Float64=floatmax(), initial_step_num::Int64=0, climbing_image::Int64=-1, initial_ΔEₘₐₓ::Float64=1.0, final_ΔEₘₐₓ::Float64=2.5e-3, remesh_spacing::Int64=1, mult_val::Float64=1.0, max_steps::Int64=-1)
    if climbing_image ≠ -1
        error("CI-NEB is not properly implemented at the moment")
    end
    
    M = length(rem_patches_arr)

    constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches_arr)
    neb_indices_arr = _get_placement_in_gradient(coords_arr, rem_patches_arr, neb_coords_place)
    pos_arr = [pos(coords) for coords in coords_arr]

    constrain_v = length(constrain_v_vol_patches[1]) > 0
    constrain_Φ = length(constrain_Φ_vol_patches[1]) > 0

    function force(x, i)
        x_to_coords!(x, coords_arr[i], vol_patches_arr[i], update_en=false)
        grad, grad_Δsing = gradient_neb_parallel(indep_coords_arr[i])
        # don't mingel with endpoints
        if (i == 1) || (i == M)
            return grad .+ grad_Δsing
        end
        # get the tangent
        Rᵢ = pos_arr[i][neb_indices_arr[i]]
        Rᵢ₊₁ = pos_arr[i+1][neb_indices_arr[i+1]]
        Rᵢ₋₁ = pos_arr[i-1][neb_indices_arr[i-1]]
        τ₊ = Rᵢ₊₁ - Rᵢ
        τ₋ = Rᵢ - Rᵢ₋₁
        if (energy_arr[i+1] - energy_arr[i])*(energy_arr[i] - energy_arr[i-1]) < 0.0
            ΔVₘₐₓ = max(abs(energy_arr[i+1]-energy_arr[i]), abs(energy_arr[i-1]-energy_arr[i]))
            ΔVₘᵢₙ = min(abs(energy_arr[i+1]-energy_arr[i]), abs(energy_arr[i-1]-energy_arr[i]))
            τ = energy_arr[i+1]>energy_arr[i-1] ? τ₊*ΔVₘₐₓ + τ₋*ΔVₘᵢₙ : τ₊*ΔVₘᵢₙ + τ₋*ΔVₘₐₓ
        else
            τ = energy_arr[i+1]>energy_arr[i-1] ? τ₊ : τ₋
        end
        τ /= norm(τ)
        # change the part of the grad that relates to the NEB
        Fₜ = dot(grad[neb_indices_arr[i]], τ)*τ
        if i ≠ climbing_image
            # for non-climbing image
            # 1) keep only the normal direction
            grad[neb_indices_arr[i]] -= Fₜ
            # 2) add the spring force in the tangent direction
            # following the "Modified NEB method" 
            # Maras, Trushin, Stukowski, Ala-Nissila, Jonsson, Comp Phys Comm, 205, 13-21 (2016)
            cosΦᵢ = dot(τ₊, τ₋) / (norm(τ₊)*norm(τ₋))
            Δ = norm(τ₊)-norm(τ₋)
            Fₘ = 0.5*(1.0 + cos(π*cosΦᵢ))*(τ₊-τ₋)
            Fₘ -= dot(Fₘ, τ)*τ
            grad[neb_indices_arr[i]] -= K*(Δ*τ + Fₘ)
        else
            # TODO : FIX CI-NEB, need to invert only the original helfrich
            # for a climbing image invert the forces in the tangent direction
            grad[neb_indices_arr[i]] -= 2.0*Fₜ
        end
        return grad .+ grad_Δsing
    end

    function force_Φ(x, i)
        x_to_coords!(x, coords_arr[i], vol_patches_arr[i], update_en=false)
        grad = gradient_Φ_parallel(indep_coords_arr[i])
        return grad
    end

    function set(x, i)
        pos_arr[i] = x
        x_to_coords!(x, coords_arr[i], vol_patches_arr[i], update_en=true)
        # update the energy array for the next NEB_force calculation
        energy_arr[i] = energy(vol_patches_arr[i], whole=true)
    end

    force_arr = [force_(x) = force(x, i) for i in 1:M]
    force_Φ_arr = [force_(x) = force_Φ(x, i) for i in 1:M]
    set_arr = [set_(x) = set(x, i) for i in 1:M]

    function _energy()
        min_domain_en_sub = -DOMAINS_WIDTH*DOMAINS_ENF_STR*0.25*DOMAINS
        return [sum([vol_patch.total_energy for vol_patch in vol_patches_arr[i]]) + min_domain_en_sub for i in 1:M]
    end

    # update constraints only if a value is given for the constraint
    if μᵥ ≠ floatmax() || μᵩ ≠ floatmax()
        # set to zero if the value of the constraint wasn't defined
        if μᵥ == floatmax()
            μᵥ = 0.0
        end 
        if μᵩ == floatmax()
            μᵩ = 0.0
        end
        c_tot = _update_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, μᵥ, μᵩ)
    else
        # otherwise don't change the constraints, but to ensure that everything works well just multiplie the constraints by 1.0
        μᵥ, μᵩ = 0.0, 0.0
        c_tot = _multiplie_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, 1.0)
    end

    energy_arr = [energy(vol_patches_arr[i], whole=true) for i in 1:M]
    prev_energy_arr = copy(energy_arr)
    
    energy_not_whole_arr = [energy(vol_patches_arr[i], whole=false, include_constr=false) for i in 1:M]
    prev_energy_not_whole_arr = copy(energy_not_whole_arr)

    F = 1.0
    v0, dtcurr, Np, αcurr = nothing, nothing, nothing, nothing

    curr_ΔEₘₐₓ = initial_ΔEₘₐₓ
    counter___ = 0
    counter__ = initial_step_num

    if climbing_image ≠ -1
        c_tot = 10.0
    end

    function force_neb_(x, i)
        x_to_coords!(x, coords_arr[i], vol_patches_arr[i], update_en=false)
        # don't mingel with endpoints
        if (i == 1) || (i == M)
            return 0.0
        end
        # get the tangent
        Rᵢ = pos_arr[i][neb_indices_arr[i]]
        Rᵢ₊₁ = pos_arr[i+1][neb_indices_arr[i+1]]
        Rᵢ₋₁ = pos_arr[i-1][neb_indices_arr[i-1]]
        τ₊ = Rᵢ₊₁ - Rᵢ
        τ₋ = Rᵢ - Rᵢ₋₁
        if (energy_arr[i+1] - energy_arr[i])*(energy_arr[i] - energy_arr[i-1]) < 0.0
            ΔVₘₐₓ = max(abs(energy_arr[i+1]-energy_arr[i]), abs(energy_arr[i-1]-energy_arr[i]))
            ΔVₘᵢₙ = min(abs(energy_arr[i+1]-energy_arr[i]), abs(energy_arr[i-1]-energy_arr[i]))
            τ = energy_arr[i+1]>energy_arr[i-1] ? τ₊*ΔVₘₐₓ + τ₋*ΔVₘᵢₙ : τ₊*ΔVₘᵢₙ + τ₋*ΔVₘₐₓ
        else
            τ = energy_arr[i+1]>energy_arr[i-1] ? τ₊ : τ₋
        end
        τ /= norm(τ)
        
        Δ = norm(τ₊)-norm(τ₋)
        return -K*Δ
    end

    while true & ((counter___ < max_steps) || (max_steps == -1))
        c_tot = _multiplie_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, mult_val, update_λs=true, min_value=1.0e4)

        plt.plot([force_neb_(pos_arr[i], i) for i in 1:M], marker="o")
        plt.savefig("curr_neb.png")
        plt.clf()

        plt.plot(_energy(), label="total", marker="o")
        plt.xlabel("ξ")
        plt.ylabel("F [pN nm]")
        plt.savefig("curr_energy.png")
        plt.clf()
        plt.plot([calc_energy(vol_patches_arr[i], "Φing") - 0.25*DOMAINS_ENF_STR*DOMAINS_WIDTH for i in 1:M], label="Φing", marker="o")
        plt.xlabel("ξ")
        plt.ylabel("F [pN nm]")
        plt.savefig("curr_Φing.png")
        plt.clf()
        plt.plot(prev_energy_arr.-energy_arr, marker="o")
        plt.plot(prev_energy_not_whole_arr.-energy_not_whole_arr, marker="x")
        plt.xlabel("ξ")
        plt.ylabel("ΔE [pN nm]")
        plt.savefig("curr_diff.png")
        plt.clf()
        
        if counter__ % 500 == 0
            plot_energy_path([energy(vol_patches, whole=false, include_constr=false) for vol_patches in vol_patches_arr], coords_arr, indep_coords_arr, neb_indices_arr)
            en_arr = [energy(vol_patches_arr[i], whole=false, include_constr=false) for i in 1:M]
            plt.plot(en_arr, label="w/o constr", marker="o")
            en_arr_wo_Φ = _energy()
            plt.plot(en_arr_wo_Φ, label="w constr w/o Φing", marker="o")
            plt.legend()
            plt.xlabel("ξ")
            plt.ylabel("F [pN nm]")
            plt.savefig(EN_DIR*"/energy_"*string(counter___)*".svg")
            plt.clf()

            en_arr_parts = Dict()
            for which in ["bending", "tilting", "stretching", "dsing", "selfing", "Φing"]
                en_arr_parts[which] = [calc_energy(vol_patches_arr[i], which) for i in 1:M]
                plt.plot(en_arr_parts[which].-en_arr_parts[which][1], label=which, marker="o")
            end
            plt.legend()
            plt.xlabel("ξ")
            plt.ylabel("ΔF [pN nm]")
            plt.savefig(EN_SEP_DIR*"/energy_sep_"*string(counter___)*".svg")
            plt.clf()

            dist_arr = [norm(pos_arr[i][neb_indices_arr[i]].-pos_arr[i+1][neb_indices_arr[i+1]]) for i in 1:(M-1)]
            plt.plot(dist_arr, marker="o")
            plt.xlabel("ξ")
            plt.savefig(EN_SEP_DIR*"/dists_"*string(counter___)*".svg")
            plt.clf()

            for j in 1:length(neb_indices_arr[1])
                path_ = [pos_arr[i][neb_indices_arr[i]][j] for i in 1:M]
                plt.plot(path_, marker="o")
            end
            plt.xlabel("ξ")
            plt.ylabel("x [nm]")
            plt.savefig(EN_SEP_DIR*"/path_"*string(counter___)*".svg")
            plt.clf()

            plt.plot(prev_energy_arr.-energy_arr, marker="o")
            plt.xlabel("ξ")
            plt.ylabel("ΔE [pN nm]")
            plt.savefig(EN_SEP_DIR*"/energy_diff_"*string(counter___)*".svg")
            plt.clf()

            # save all the data from above as a csv
            CSV.write(EN_DIR*"/energy"*string(counter___)*".csv", (Energy=en_arr, Energy_with_constr=en_arr_wo_Φ, En_diff=prev_energy_arr.-energy_arr, Bending=en_arr_parts["bending"], Tilting=en_arr_parts["tilting"], Stretching=en_arr_parts["stretching"], dsing=en_arr_parts["dsing"], Selfing=en_arr_parts["selfing"], Φing=en_arr_parts["Φing"]))

            if constrain_v
                idxs = [vol_patch.idx for vol_patch in constrain_v_vol_patches[1]]
                c0 = Vector{Float64}[]
                for vol_patches in constrain_v_vol_patches
                    push!(c0, Float64[])
                    for vol_patch in vol_patches
                        c = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
                        push!(c0[end], c)
                    end
                end

                c0 = reduce(hcat, c0)
                

                for (idx, c_) in zip(idxs, eachrow(c0))
                    plt.plot(c_, label=string(idx), marker="o")
                end
                plt.legend()
                plt.savefig(CONSTR_DIR*"/c0_"*string(counter___)*".svg")
                plt.clf()
            end

            if constrain_Φ
                c1 = Vector{Float64}[]
                idxs = [vol_patch.idx for vol_patch in constrain_Φ_vol_patches[1]]
                for vol_patches in constrain_Φ_vol_patches
                    push!(c1, Float64[])
                    for vol_patch in vol_patches
                        c = (vol_patch.total_Φ/vol_patch.total_area - vol_patch.Φ₀)
                        push!(c1[end], c)
                    end
                end
                c1 = reduce(hcat, c0)

                for (idx, c_) in zip(idxs, eachrow(c1))
                    plt.plot(c_, label=string(idx), marker="o")
                end
                plt.legend()
                plt.savefig(CONSTR_DIR*"/c1_"*string(counter___)*".svg")
                plt.clf()
            end

            plot_reaction(rem_patches_arr)
            for i in 1:M
                plot_coords(vol_patches_arr[i])
            end
            save_string(coords_arr, vol_patches_arr, rem_patches_arr, STATES_DIR*string(counter___))
        end

        # skip steps in high likelihood if the constraints are pretty well defined.
        avoid_remesh = (counter___ > 1) && ((c_tot < 1.0e-2) || !(constrain_v || constrain_Φ)) && (rand() < 0.5) && (!isnan(F))
        
        # remeshing if previous step was threw a NaN or by the schedule
        if (!avoid_remesh) & ((((counter__ == 0) || (counter__ > 9)) && (counter__ % (remesh_spacing*(10^minimum([3,floor(log10(counter__ + 1))]))) == 0)) || isnan(F))
            for i in 1:M
                coords_arr[i], indep_coords_arr[i], vol_patches_arr[i], rem_patches_arr[i] = remesh(rem_patches_arr[i])
                @static if DOMAINS
                    _minimize_Φ_for_domains(rem_patches_arr[i])
                    x_to_coords!(pos(coords_arr[i]), coords_arr[i], vol_patches_arr[i], update_en=true)
                end
            end
            
            constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches_arr)
            
            neb_indices_arr = _get_placement_in_gradient(coords_arr, rem_patches_arr, neb_coords_place)
            energy_arr = [energy(vol_patches_arr[i], whole=true) for i in 1:M]
            if any([vol_patch.use_Φ for vol_patch in vol_patches_arr[1]])
                println("Remeshed... Minimizing Φs")
                _, _, _, _, _, _ = fire2(force_Φ_arr, set_arr, [pos(coords) for coords in coords_arr], max_steps=100, v0=nothing, dtcurr = nothing, Np = nothing, αcurr = nothing)
                v0, dtcurr, Np, αcurr = nothing, nothing, nothing, nothing
                #x, _, _ = quickmin(force_Φ_arr, set_arr, [pos(coords) for coords in coords_arr], max_steps=100, dt=1.0e-3)
                #v0, F = nothing, 1.0  
            end

            v0 = nothing

            plt.plot(_energy(), label="total", marker="o")
            plt.xlabel("ξ")
            plt.ylabel("F [pN nm]")
            plt.savefig("curr_energy.png")
            plt.clf()
            plt.plot([calc_energy(vol_patches_arr[i], "Φing") - 0.25*DOMAINS_ENF_STR*DOMAINS_WIDTH*DOMAINS for i in 1:M], label="Φing", marker="o")
            plt.xlabel("ξ")
            plt.ylabel("F [pN nm]")
            plt.savefig("curr_Φing.png")
            plt.clf()  
        end

        pos_arr = [pos(coords) for coords in coords_arr]
        energy_arr = [energy(vol_patches_arr[i], whole=true) for i in 1:M]
        prev_energy_arr = copy(energy_arr)
        energy_not_whole_arr = [energy(vol_patches_arr[i], whole=false, include_constr=false) for i in 1:M]
        prev_energy_not_whole_arr = copy(energy_not_whole_arr)
        x, v0, dtcurr, Np, αcurr, F = fire2(force_arr, set_arr, [pos(coords) for coords in coords_arr], max_steps=100, v0=v0, dtcurr = dtcurr, Np = Np, αcurr = αcurr)
        #x, v0, F = quickmin(force_arr, set_arr, [pos(coords) for coords in coords_arr], max_steps=100, v0=v0, dt=1.0e-3)
        energy_arr = [energy(vol_patches_arr[i], whole=true) for i in 1:M]
        energy_not_whole_arr = [energy(vol_patches_arr[i], whole=false, include_constr=false) for i in 1:M]
                    
        counter__ += 1
        counter___ += 1

        println("counter__=", counter__, " counter___=", counter___, " μᵥ=", round(μᵥ, digits=2), " μᵩ=", round(μᵩ, digits=2), " F=", round(minimum(F), digits=5), " ΔE=", round(maximum(abs.(prev_energy_arr.-energy_arr)), digits=5))
        flush(stdout)
        
        dE = maximum(abs.(prev_energy_arr.-energy_arr))
        not_done = (c_tot > 1.0e-4) || (dE > final_ΔEₘₐₓ)
        reset = (dE < curr_ΔEₘₐₓ && counter__ >= 5 && curr_ΔEₘₐₓ > final_ΔEₘₐₓ) || (counter__ >= 10000 && curr_ΔEₘₐₓ == final_ΔEₘₐₓ)

        if !not_done || (!(constrain_v || constrain_Φ) && curr_ΔEₘₐₓ == final_ΔEₘₐₓ)
            break
        end
        
        if !reset
            continue
        end

        curr_ΔEₘₐₓ = max(curr_ΔEₘₐₓ * 0.75, final_ΔEₘₐₓ)
        counter__ = 0
    end

    plot_energy_path([energy(vol_patches, whole=false, include_constr=false) for vol_patches in vol_patches_arr], coords_arr, indep_coords_arr, neb_indices_arr)

    # final plots of all the configurations
    for i in 1:M
        plot_coords(vol_patches_arr[i])
    end
    return coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr
end

function interpolate_wrt_neb(coords_arr, rem_patches_arr, neb_coords_place::Vector{NEBCoordinatePlace}, ; path_distribution=nothing)
    # TODO: fix to also interpolate λs
    error("interpolate_wrt_neb is not properly implemented at the moment")
    
    neb_indices_arr = _get_placement_in_gradient(coords_arr, rem_patches_arr, neb_coords_place)
    pos_arr = [pos(coords)[neb_indices] for (coords, neb_indices) in zip(coords_arr, neb_indices_arr)]
    
    M = length(rem_patches_arr)
    # find distance between subsequent points
    ds_arr = [0.0]
    for (pos_a, pos_b) in zip(pos_arr[1:end-1], pos_arr[2:end])
        ds_ = 0.0
        # distance is measure by the midplane only for simplicity
        for (c_a, c_b) in zip(pos_a, pos_b)
            ds_ += (c_a-c_b)*(c_a-c_b)
        end
        ds_ = sqrt(ds_)
        push!(ds_arr, ds_)
    end
    ds_arr = cumsum(ds_arr)
    new_ds = path_distribution === nothing ? LinRange(0.0, ds_arr[end], M) : copy(path_distribution).*ds_arr[end]

    N = 501
    K = length(rem_patches_arr[1])
    # remesh all patches to be subdivided to 501 points
    rem_patches_div = [remesh_N(rem_patches, N=N)[end] for rem_patches in rem_patches_arr]
    # the ones for change, easier way to make them
    rem_patches_arr_for_change = [remesh_N(rem_patches, N=N)[end] for rem_patches in rem_patches_arr]

    for i in 1:K
        # interpolate each constraint with respect to the distance
        for j in 1:N
            if rem_patches_div[1][i].upper_vol_patch ≠ nothing
                # upper vol patch
                # volume constraint
                λᵥs = [rem_patches[i].upper_vol_patch.λᵥ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵥs = Spline1D(ds_arr, λᵥs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].upper_vol_patch.λᵥ = spl_λᵥs(ds)
                end
                # composition constraint
                λᵩs = [rem_patches[i].upper_vol_patch.λᵩ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵩs = Spline1D(ds_arr, λᵩs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].upper_vol_patch.λᵩ = spl_λᵩs(ds)
                end
            end
            if rem_patches_div[1][i].lower_vol_patch ≠ nothing
                # lower vol patch
                # volume constraint
                λᵥs = [rem_patches[i].lower_vol_patch.λᵥ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵥs = Spline1D(ds_arr, λᵥs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].lower_vol_patch.λᵥ = spl_λᵥs(ds)
                end
                # composition constraint
                λᵩs = [rem_patches[i].lower_vol_patch.λᵩ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵩs = Spline1D(ds_arr, λᵩs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].lower_vol_patch.λᵩ = spl_λᵩs(ds)
                end
            end
            if rem_patches_div[1][i].midplane_area_patch ≠ nothing
                # midplane area patch
                # area constraint
                λᵥs = [rem_patches[i].midplane_area_patch.λᵥ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵥs = Spline1D(ds_arr, λᵥs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].midplane_area_patch.λᵥ = spl_λᵥs(ds)
                end
                # composition constraint
                λᵩs = [rem_patches[i].midplane_area_patch.λᵩ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵩs = Spline1D(ds_arr, λᵩs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].midplane_area_patch.λᵩ = spl_λᵩs(ds)
                end
            end
            if rem_patches_div[1][i].upper_water_patch ≠ nothing
                # upper water patch
                # volume constraint
                λᵥs = [rem_patches[i].upper_water_patch.λᵥ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵥs = Spline1D(ds_arr, λᵥs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].upper_water_patch.λᵥ = spl_λᵥs(ds)
                end
                # composition constraint
                λᵩs = [rem_patches[i].upper_water_patch.λᵩ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵩs = Spline1D(ds_arr, λᵩs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].upper_water_patch.λᵩ = spl_λᵩs(ds)
                end
            end
            if rem_patches_div[1][i].lower_water_patch ≠ nothing
                # lower water patch
                # volume constraint
                λᵥs = [rem_patches[i].lower_water_patch.λᵥ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵥs = Spline1D(ds_arr, λᵥs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].lower_water_patch.λᵥ = spl_λᵥs(ds)
                end
                # composition constraint
                λᵩs = [rem_patches[i].lower_water_patch.λᵩ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_λᵩs = Spline1D(ds_arr, λᵩs)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    rem_patches[i].lower_water_patch.λᵩ = spl_λᵩs(ds)
                end
            end

            rs = [rem_patches[i].mid_coords[j].r for rem_patches in rem_patches_div]
            zs = [rem_patches[i].mid_coords[j].z for rem_patches in rem_patches_div]
            # make cubic splines
            spl_rs = Spline1D(ds_arr, rs)
            spl_zs = Spline1D(ds_arr, zs)
            for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                if !rem_patches[i].mid_coords[j].fixed_r
                    rem_patches[i].mid_coords[j].r = spl_rs(ds)
                end
                if !rem_patches[i].mid_coords[j].fixed_z
                    rem_patches[i].mid_coords[j].z = spl_zs(ds)
                end
            end
        end
        if rem_patches_div[1][i].upper_coords ≠ nothing
            for j in 1:N
                rs = [rem_patches[i].upper_coords[j].r for rem_patches in rem_patches_div]
                zs = [rem_patches[i].upper_coords[j].z for rem_patches in rem_patches_div]
                Φ₊s = [rem_patches[i].upper_coords[j].Φ₊ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_rs = Spline1D(ds_arr, rs)
                spl_zs = Spline1D(ds_arr, zs)
                spl_Φ₊s = Spline1D(ds_arr, Φ₊s)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    if !rem_patches[i].upper_coords[j].fixed_r
                        rem_patches[i].upper_coords[j].r = spl_rs(ds)
                    end
                    if !rem_patches[i].upper_coords[j].fixed_z
                        rem_patches[i].upper_coords[j].z = spl_zs(ds)
                    end
                    if !rem_patches[i].upper_coords[j].fixed_Φ₊
                        rem_patches[i].upper_coords[j].Φ₊ = spl_Φ₊s(ds)
                    end
                end
            end
        end
        if rem_patches_div[1][i].lower_coords ≠ nothing
            for j in 1:N
                rs = [rem_patches[i].lower_coords[j].r for rem_patches in rem_patches_div]
                zs = [rem_patches[i].lower_coords[j].z for rem_patches in rem_patches_div]
                Φ₊s = [rem_patches[i].lower_coords[j].Φ₊ for rem_patches in rem_patches_div]
                # make cubic splines
                spl_rs = Spline1D(ds_arr, rs)
                spl_zs = Spline1D(ds_arr, zs)
                spl_Φ₊s = Spline1D(ds_arr, Φ₊s)
                for (ds, rem_patches) in zip(new_ds, rem_patches_arr_for_change)
                    if !rem_patches[i].lower_coords[j].fixed_r
                        rem_patches[i].lower_coords[j].r = spl_rs(ds)
                    end
                    if !rem_patches[i].lower_coords[j].fixed_z
                        rem_patches[i].lower_coords[j].z = spl_zs(ds)
                    end
                    if !rem_patches[i].lower_coords[j].fixed_Φ₊
                        rem_patches[i].lower_coords[j].Φ₊ = spl_Φ₊s(ds)
                    end
                end
            end
        end
    end

    # create the arrays of the true remeshing patches
    coords_arr = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
    indep_coords_arr = Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}}[]
    vol_patches_arr = Vector{VolumePatch}[]
    rem_patches_arr = Vector{RemeshingPatch}[]
    # remesh everything now
    for interpolated_rem_patches in rem_patches_arr_for_change
        new_coords, new_indep_coords, new_vol_patches, new_rem_patches = remesh(interpolated_rem_patches)
        # add them to the correct arrays
        push!(coords_arr, new_coords)
        push!(indep_coords_arr, new_indep_coords)
        push!(vol_patches_arr, new_vol_patches)
        push!(rem_patches_arr, new_rem_patches)
    end
    return coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr
end