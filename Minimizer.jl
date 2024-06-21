# saving directory for the figures
const FIG_DIR = MAIN_PATH*"figs"
if !ispath(FIG_DIR)
    mkpath(FIG_DIR)
end

const CONSTR_DIR = MAIN_PATH*"constr_figs"
if !ispath(CONSTR_DIR)
    mkpath(CONSTR_DIR)
end

const EN_DIR = MAIN_PATH*"en/"
if !ispath(EN_DIR)
    mkdir(EN_DIR)
end

const EN_SEP_DIR = MAIN_PATH*"en_sep/"
if !ispath(EN_SEP_DIR)
    mkdir(EN_SEP_DIR)
end

const STATES_DIR = MAIN_PATH*"states/"
if !ispath(STATES_DIR)
    mkdir(STATES_DIR)
end

const ENERGY_FILE = MAIN_PATH*"energy.csv"
# make the file, delete it first if exists
if isfile(ENERGY_FILE)
    rm(ENERGY_FILE)
end
open(ENERGY_FILE, "w") do io
    write(io, "step,energy\n")
end

function save_energy(step::Int64, en::Float64)
    open(ENERGY_FILE, "a") do io
        write(io, "$step,$en\n")
    end
end

function quickmin(f, set, x0::Vector{Vector{Float64}}, ; v0::Union{Vector{Vector{Float64}}, Nothing}=nothing, max_steps::Int64 = 1000, dt::Float64 = 1.0e-3)
    f0 = [-f_(x0_) for (f_, x0_) in zip(f, x0)]
    if v0 === nothing
        v0 = [zeros(length(x0_)) for x0_ in x0]
    end
    f0_normalized = [f_/norm(f_) for f_ in f0]
    x = deepcopy(x0)
    M = length(f)
    f0ₘₐₓ = 0.0
    failed = false
    @showprogress 1 "Evolving Structure" for _ in 1:max_steps
        f0ₘₐₓ = 0.0
        # update positions
        for i in 1:M
            v0[i] = max(dot(v0[i], f0_normalized[i]), 0.0)*f0_normalized[i]
            x[i] = x[i] + v0[i].*dt
            set[i](x[i])
            v0[i] += f0[i].*dt
        end
        # update search direction
        for i in 1:M
            f_new = -f[i](x[i])
            norm_f_new = norm(f_new)
            if isnan(norm_f_new)
                failed = true
                break
            end
            f0[i] = f_new
            f0_normalized[i] = f_new/norm_f_new
            
            norm_f_new /= sqrt(length(f_new))
            if f0ₘₐₓ < norm_f_new
                f0ₘₐₓ = norm_f_new
            end
        end
        if failed
            break
        end
    end
    return (x, v0, f0ₘₐₓ)
end

function quickmin(f, set, x0::Vector{Float64}, ; v0::Union{Vector{Float64}, Nothing} = nothing, max_steps::Int64 = 1000, dt::Float64 = 1.0e-4)
  x, v0, f0_norm = quickmin([f], [set], [x0], v0= v0 === nothing ? nothing : [v0], max_steps=max_steps, dt=dt)
  return (x[1], v0[1], f0_norm)
end

function fire2(f, set, x0::Vector{Vector{Float64}}, ; v0::Union{Vector{Vector{Float64}}, Nothing}=nothing, max_steps::Int64 = 1000, tmin::Float64 = 1.0e-5, tmax::Float64 = 1.0e-1, dtcurr::Union{Nothing,Vector{Float64}}=nothing, Np::Union{Nothing, Vector{Int64}}=nothing, delaystep::Int64=5, dtgrow::Float64=1.1, dtshrink::Float64=0.5, α₀::Float64=0.25, αshrink::Float64=0.99, αcurr::Union{Nothing, Vector{Float64}}=nothing)
    M = length(f)
    αcurr = αcurr === nothing ? fill(α₀, M) : αcurr
    dtcurr = dtcurr === nothing ? fill(tmin, M) : dtcurr
    Np = Np === nothing ? fill(0, M) : Np

    f0 = [-f_(x0_) for (f_, x0_) in zip(f, x0)]
    if v0 === nothing
        v0 = [zeros(length(x0_)) for x0_ in x0]
    end
    f0_normalized = [f_/norm(f_) for f_ in f0]
    x = deepcopy(x0)
    failed = false
    f0ₘₐₓ = 0.0
    @showprogress 1 "Evolving Structure" for _ in 1:max_steps
        f0ₘₐₓ = 0.0
        # update positions
        for i in 1:M
            p = dot(v0[i], f0_normalized[i])
            if p > 0.0
                Np[i] += 1
                if Np[i] > delaystep
                    dtcurr[i] = min(dtcurr[i]*dtgrow, tmax)
                    αcurr[i] *= αshrink
                end
            else
                Np[i] = 0
                dtcurr[i] = max(dtcurr[i]*dtshrink, tmin)
                αcurr[i] = α₀
                # up-hill correction
                x[i] = x[i] - 0.5*dtcurr[i]*v0[i]
                v0[i] *= 0.0
            end
        end
        for i in 1:M
            # semi-implicit euler
            v0[i] += f0[i].*dtcurr[i]
            v0[i] = (1.0 - αcurr[i])*v0[i] + αcurr[i]*norm(v0[i])*f0_normalized[i] # mixing
            x[i] = x[i] + v0[i].*dtcurr[i]
        end
        for i in 1:M
            set[i](x[i])
            f_new = -f[i](x[i])
            norm_f_new = norm(f_new)
            if isnan(norm_f_new)
                failed = true
                break
            end
            f0[i] = f_new
            f0_normalized[i] = f_new/norm_f_new
            norm_f_new /= sqrt(length(f_new))
            if f0ₘₐₓ < norm_f_new
                f0ₘₐₓ = norm_f_new
            end
        end
        if failed
            break
        end
    end

    return (x, v0, dtcurr, Np, αcurr, f0ₘₐₓ)
end

function fire2(f, set, x0::Vector{Float64}, ; v0::Union{Vector{Float64}, Nothing}=nothing, max_steps::Int64 = 1000, tmin::Float64 = 1.0e-5, tmax::Float64 = 1.0e-1, dtcurr::Union{Nothing,Float64}=nothing, Np::Union{Nothing, Int64}=nothing, delaystep::Int64=5, dtgrow::Float64=1.1, dtshrink::Float64=0.5, α₀::Float64=0.25, αshrink::Float64=0.99, αcurr::Union{Nothing, Float64}=nothing)
    x, v0, dtcurr, Np, αcurr, f0ₘₐₓ = fire2([f], [set], [x0], v0=v0 === nothing ? nothing : [v0], max_steps=max_steps, tmin=tmin, tmax=tmax, dtcurr=dtcurr === nothing ? nothing : [dtcurr], Np=Np === nothing ? nothing : [Np], delaystep=delaystep, dtgrow=dtgrow, dtshrink=dtshrink, α₀=α₀, αshrink=αshrink, αcurr=αcurr === nothing ? nothing : [αcurr])
    return (x[1], v0[1], dtcurr[1], Np[1], αcurr[1], f0ₘₐₓ)
end

# todo: REIMPLEMENT CG...

function energy(vol_patches_arr::Vector{Vector{VolumePatch}}, ; whole=true, include_constr=true)
    en = DOMAINS ? -DOMAINS_WIDTH*DOMAINS_ENF_STR*0.25*length(vol_patches_arr)*whole : 0.0 # removes the extra energy from a domain just so it'll be easier to compare the total energy to the elastic one!
    for vol_patches in vol_patches_arr
        for vol_patch in vol_patches
            if include_constr ⊻ whole
                penalty = 0.0
                if vol_patch.constrain_v
                    c = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
                    λᵥc = vol_patch.λᵥ*c
                    if λᵥc < 0.0
                        λᵥc *= exp(-50.0*abs(c))
                    end
                    penalty += λᵥc + 0.5*vol_patch.μᵥ*c*c
                end
                if vol_patch.constrain_Φ
                    c = (vol_patch.total_Φ/vol_patch.total_area - vol_patch.Φ₀)
                    λᵩc = vol_patch.λᵩ*c
                    if λᵩc < 0.0
                        λᵩc *= exp(-50.0*abs(c))
                    end
                    penalty += λᵩc + 0.5*vol_patch.μᵩ*c*c
                end
            end
            if whole
                en += vol_patch.total_energy
                if !include_constr
                    en -= penalty
                end
                continue
            end
            en += calc_energy(vol_patch)
            if include_constr
                en += penalty
            end
        end
    end
    return en / length(vol_patches_arr)
end

energy(vol_patches::Vector{VolumePatch}, ; whole=true, include_constr=true) = energy([vol_patches], whole=whole, include_constr=include_constr)

function _get_constrained_patches(vol_patches_arr::Vector{Vector{VolumePatch}})
    constrain_v_vol_patches = Vector{VolumePatch}[]
    constrain_Φ_vol_patches = Vector{VolumePatch}[]
    
    for vol_patches in vol_patches_arr
        push!(constrain_v_vol_patches, VolumePatch[])
        push!(constrain_Φ_vol_patches, VolumePatch[])
        for vol_patch in vol_patches
            if vol_patch.constrain_v
                push!(constrain_v_vol_patches[end], vol_patch)
            end
            if vol_patch.constrain_Φ
                push!(constrain_Φ_vol_patches[end], vol_patch)
            end
        end
    end
    return constrain_v_vol_patches, constrain_Φ_vol_patches
end

function _get_constrained_patches(vol_patches::Vector{VolumePatch}) 
    temp = _get_constrained_patches([vol_patches])
    return temp[1][1], temp[2][1]
end

function _update_constraints(constrain_v_vol_patches::Vector{Vector{VolumePatch}}, constrain_Φ_vol_patches::Vector{Vector{VolumePatch}}, μᵥ::Float64, μᵩ::Float64, ; update_λs = false, update_μs = true, verbose = false)
    if verbose
        println("")
    end
    c0 = Vector{Float64}[]
    for vol_patches in constrain_v_vol_patches
        push!(c0, Float64[])
        for vol_patch in vol_patches
            c = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
            push!(c0[end], c)
            if !update_μs
                μᵥ = vol_patch.μᵥ
            end
            if update_λs
                λᵥc = vol_patch.λᵥ*c
                if λᵥc < 0.0
                    vol_patch.λᵥ *= exp(-50.0*abs(c))
                end
                vol_patch.λᵥ += vol_patch.μᵥ * c
                if verbose
                    println("μᵥ=", μᵥ, " c=", c, " λᵥ=", vol_patch.λᵥ)
                end
            else
                if verbose
                    println("μᵥ=", μᵥ, " c=", c)
                end
            end
            vol_patch.μᵥ = μᵥ
            update!(vol_patch)
        end
    end
    c1 = Vector{Float64}[]
    for vol_patches in constrain_Φ_vol_patches
        push!(c1, Float64[])
        for vol_patch in vol_patches
            c = (vol_patch.total_Φ/vol_patch.total_area - vol_patch.Φ₀)
            push!(c1[end], c)
            if !update_μs
                μᵩ = vol_patch.μᵩ
            end
            if update_λs                
                λᵩc = vol_patch.λᵩ*c
                if λᵩc < 0.0
                    vol_patch.λᵩ *= exp(-50.0*abs(c))
                end
                vol_patch.λᵩ += vol_patch.μᵩ * c
                if verbose
                    println("μᵩ=", μᵩ, " c=", c, " λᵩ=", vol_patch.λᵩ)
                end
            else
                if verbose
                    println("μᵩ=", μᵩ, " c=", c)
                end
            end
            vol_patch.μᵩ = μᵩ
            update!(vol_patch)
        end
    end
    temp = abs.(vcat(vcat(c0...), vcat(c1...)))
    return length(temp) == 0 ? 1.0 : maximum(temp)
end

function _multiplie_constraints(constrain_v_vol_patches::Vector{Vector{VolumePatch}}, constrain_Φ_vol_patches::Vector{Vector{VolumePatch}}, factor::Float64, ; update_λs = false, verbose = false, max_value = floatmax(), min_value = 0.0)
    if verbose
        println("")
    end
    c0 = Vector{Float64}[]
    for vol_patches in constrain_v_vol_patches
        push!(c0, Float64[])
        for vol_patch in vol_patches
            c = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
            push!(c0[end], c)
            if update_λs
                λᵥc = vol_patch.λᵥ*c
                if λᵥc < 0.0
                    vol_patch.λᵥ *= exp(-50.0*abs(c))
                end
                vol_patch.λᵥ += vol_patch.μᵥ * c
                if verbose
                    println("μᵥ=", vol_patch.μᵥ, " c=", c, " λᵥ=", vol_patch.λᵥ)
                end
            else
                if verbose
                    println("μᵥ=", vol_patch.μᵥ, " c=", c)
                end
            end
            vol_patch.μᵥ = max(min(vol_patch.μᵥ*factor, max_value), min_value)
            update!(vol_patch)
        end
    end
    c1 = Vector{Float64}[]
    for vol_patches in constrain_Φ_vol_patches
        push!(c1, Float64[])
        for vol_patch in vol_patches
            c = (vol_patch.total_Φ/vol_patch.total_area - vol_patch.Φ₀)
            push!(c1[end], c)
            if update_λs
                λᵩc = vol_patch.λᵩ*c
                if λᵩc < 0.0
                    vol_patch.λᵩ *= exp(-50.0*abs(c))
                end
                vol_patch.λᵩ += vol_patch.μᵩ * c
                if verbose
                    println("μᵩ=", vol_patch.μᵩ, " c=", c, " λᵩ=", vol_patch.λᵩ)
                end
            else
                if verbose
                    println("μᵩ=", vol_patch.μᵩ, " c=", c)
                end
            end
            vol_patch.μᵩ = max(min(vol_patch.μᵩ*factor, max_value), min_value)
            update!(vol_patch)
        end
    end
    temp = abs.(vcat(vcat(c0...), vcat(c1...)))
    return length(temp) == 0 ? 1.0 : maximum(temp)
end

_update_constraints(constrain_v_vol_patches::Vector{VolumePatch}, constrain_Φ_vol_patches::Vector{VolumePatch}, μᵥ::Float64, μᵩ::Float64, ; update_λs::Bool = false, verbose::Bool = true) = _update_constraints([constrain_v_vol_patches], [constrain_Φ_vol_patches], μᵥ, μᵩ, update_λs = update_λs, verbose = verbose)
_multiplie_constraints(constrain_v_vol_patches::Vector{VolumePatch}, constrain_Φ_vol_patches::Vector{VolumePatch}, factor::Float64, ; update_λs::Bool = false, max_value::Float64 = floatmax(), min_value::Float64 = 0.0) = _multiplie_constraints([constrain_v_vol_patches], [constrain_Φ_vol_patches], factor, update_λs = update_λs, verbose = true, max_value = max_value, min_value = min_value)

function _print_info(counter__::Int64, counter___::Int64, μᵥ::Float64, μᵩ::Float64, ΔE::Float64, E1::Float64, vol_patches_arr::Vector{Vector{VolumePatch}})
    println("counter=", rpad(counter__, 4, " "), " total=", rpad(counter___, 5, " "), " μᵥ=", lpad(round(μᵥ, digits=2), 2), " μᵩ=", lpad(round(μᵩ, digits=2), 2), " ΔE=", lpad(round(ΔE*ENERGY_FACTOR, digits=5), 5), " E1=", lpad(round(E1*ENERGY_FACTOR, digits=5), 5), " energy=", lpad(round(energy(vol_patches_arr, whole=false, include_constr=false)*ENERGY_FACTOR, digits=5), 5), "(", round(energy(vol_patches_arr, include_constr=false)*ENERGY_FACTOR, digits=5),")")
    flush(stdout)
end

_print_info(counter__::Int64, counter___::Int64, μᵥ::Float64, μᵩ::Float64, ΔE::Float64, E1::Float64, vol_patches::Vector{VolumePatch}) = _print_info(counter__, counter___, μᵥ, μᵩ, ΔE, E1, [vol_patches])

function minimize(coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}, indep_coords::Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}}, vol_patches::Vector{VolumePatch}, rem_patches::Vector{RemeshingPatch}, which_to_adjust::Vector{Int64}=Int64[], when_to_adjust::Vector{Int64}=Int64[], adjust_resolution_to::Vector{Float64}=Float64[], ; μᵥ::Float64=floatmax(), μᵩ::Float64=floatmax(), initial_ΔEₘₐₓ::Float64=1.0, final_ΔEₘₐₓ::Float64=1.0e-4, remesh_spacing::Int64=1, mult_val::Float64=1.0, max_steps::Int64=-1)
    
    constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches)

    constrain_v = length(constrain_v_vol_patches) > 0
    constrain_Φ = length(constrain_Φ_vol_patches) > 0

    use_Φ = any([vol_patch.use_Φ for vol_patch in vol_patches])
    
    curr_pos = pos(coords)
    
    function force(x)
        x_to_coords!(x, coords, vol_patches, update_en=false)
        return gradient_parallel(indep_coords)
    end

    function force_Φ(x)
        x_to_coords!(x, coords, vol_patches, update_en=false)
        return gradient_Φ_parallel(indep_coords)
    end

    function set(x)
        curr_pos = x
        x_to_coords!(x, coords, vol_patches, update_en=true)
        curr_energy = energy(vol_patches, whole=true)
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
        _update_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, μᵥ, μᵩ)
    else
        # otherwise don't change the constraints, but to ensure that everything works well just multiplie the constraints by 1.0
        μᵥ, μᵩ = 0.0, 0.0
        _multiplie_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, 1.0)
    end

    curr_energy = energy(vol_patches, whole=true)
    prev_energy = curr_energy + 20.0

    max_to_adjust = length(when_to_adjust) == 0 ? 0 : maximum(when_to_adjust)

    c_tot = 1.0
    F = 1.0
    v0, dtcurr, Np, αcurr = nothing, nothing, nothing, nothing
    curr_ΔEₘₐₓ = initial_ΔEₘₐₓ
    counter___ = 0
    counter__ = 0
    while true
        c_tot = _multiplie_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, mult_val, update_λs=true, min_value = 1.0e4)
        flush(stdout)
        if counter__ % 1000 == 0
            save_energy(counter___, curr_energy)
            plot_coords(vol_patches)
            #plot_thickness(vol_patches)
            #plot_energy_density(vol_patches, "stretching")
            #save_state(coords, vol_patches, rem_patches, STATES_DIR*string(counter___)*".coords", STATES_DIR*string(counter___)*".vol_patches", STATES_DIR*string(counter___)*".rem_patches")
        end

        REMESHED = false

        # change resolution if needed
        if counter___ in when_to_adjust
            new_rem_patches = RemeshingPatch[]
            indices = findall(isequal(counter___), when_to_adjust)
            which_arr = which_to_adjust[indices]
            resolution_arr = adjust_resolution_to[indices]
            for (i, rem_patch) in enumerate(rem_patches)
                j = findfirst(isequal(i), which_arr)
                if j === nothing
                    push!(new_rem_patches, rem_patch)
                else
                    push!(new_rem_patches, change_rem_patch_resolution(rem_patch, resolution_arr[j]))
                end
            end
            coords, indep_coords, vol_patches, rem_patches = remesh(new_rem_patches)
            if use_Φ
                println("Remeshed... Minimizing Φs")
                _, _, _, _, _, _ = fire2(force_Φ, set, pos(coords), max_steps=100, v0=nothing, dtcurr = nothing, Np = nothing, αcurr = nothing)
            else
                println("Remeshed...")
            end 
            v0, dtcurr, Np, αcurr = nothing, nothing, nothing, nothing
            REMESHED = true
        end

        # randomly avoid remeshes in the constraints are low enough
        avoid_remesh = (counter___ > 1) && ((c_tot < 5.0e-3) || !(constrain_v || constrain_Φ)) && (rand() < 0.5) && (!isnan(F))
        avoid_remesh = avoid_remesh && (counter___ > 10000)

        # remeshing if previous step was threw a NaN or by the schedule
        if (!avoid_remesh) && ((((counter__ == 0) || (counter__ > 9)) && (counter__ % (remesh_spacing*(10^minimum([3,floor(log10(counter__ + 1))]))) == 0)) || isnan(F))
            coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)
            if use_Φ
                println("Remeshed... Minimizing Φs")
                _, _, _, _, _, _ = fire2(force_Φ, set, pos(coords), max_steps=100, v0=nothing, dtcurr = nothing, Np = nothing, αcurr = nothing)
            else
                println("Remeshed...")
            end 
            v0, dtcurr, Np, αcurr = nothing, nothing, nothing, nothing
            curr_energy = energy(vol_patches, whole=true)
            REMESHED = true
        end

        if REMESHED
            constrain_v_vol_patches = VolumePatch[]
            constrain_Φ_vol_patches = VolumePatch[]
            for vol_patch in vol_patches
                if vol_patch.constrain_v
                    push!(constrain_v_vol_patches, vol_patch)
                end
                if vol_patch.constrain_Φ
                    push!(constrain_Φ_vol_patches, vol_patch)
                end
            end
            curr_energy = energy(vol_patches, whole = true)
            prev_energy = curr_energy + 20.0
            v0, dtcurr, Np, αcurr, F = nothing, nothing, nothing, nothing, 1.0
        end

        counter__ += 1
        counter___ += 1

        prev_energy = energy(vol_patches, whole = true)
        x, v0, dtcurr, Np, αcurr, F = fire2(force, set, pos(coords), max_steps=100, v0=v0, dtcurr = dtcurr, Np = Np, αcurr = αcurr)
        curr_energy = energy(vol_patches, whole = true)
  
        println(calc_energy(vol_patches, "dsing"), " dsing")

        _print_info(counter__, counter___, μᵥ, μᵩ, prev_energy-curr_energy, curr_energy, vol_patches)

        dE = abs.(prev_energy.-curr_energy)
        not_done = (c_tot > 1.0e-5) || (dE > final_ΔEₘₐₓ) || (counter___ < (max_to_adjust + 100)) 
        reset = (dE < curr_ΔEₘₐₓ && counter__ >= 5 && curr_ΔEₘₐₓ > final_ΔEₘₐₓ) || (counter__ >= 10000 && curr_ΔEₘₐₓ == final_ΔEₘₐₓ)

        if !not_done || (!(constrain_v || constrain_Φ) && curr_ΔEₘₐₓ == final_ΔEₘₐₓ) || (max_steps > 0 && counter___ >= max_steps)
            break
        end
        
        if !reset
            continue
        end

        counter__ = 0
        curr_ΔEₘₐₓ = max(curr_ΔEₘₐₓ * 0.75, final_ΔEₘₐₓ)
    end

    println("minimal energy ", calc_energy(vol_patches)*ENERGY_FACTOR)
    save_energy(counter___, calc_energy(vol_patches)*ENERGY_FACTOR)

    plot_coords(vol_patches)
    plot_energy_density(vol_patches, "total")
    plot_energy_density(vol_patches, "bending")
    plot_energy_density(vol_patches, "tilting")
    plot_energy_density(vol_patches, "stretching")
    plot_energy(vol_patches, "total")
    plot_energy_dsing(vol_patches)
    save_state(coords, vol_patches, rem_patches, STATES_DIR*"final.coords", STATES_DIR*"final.vol_patches", STATES_DIR*"final.rem_patches")
    return coords, indep_coords, vol_patches, rem_patches
end

function _interpolate_helper(rem_patches_A::Vector{RemeshingPatch}, rem_patches_B::Vector{RemeshingPatch}, x::Float64, ; vol_patches_A = nothing, vol_patches_B = nothing)
    # remesh all patches to be subdivided to 201 points
    rem_patches_A_div = remesh_N(rem_patches_A, N=501)[end]
    rem_patches_B_div = remesh_N(rem_patches_B, N=501)[end]
    # interpolate from 0 (state A) to 1 (state B)
    # easiest way to make the interpolated RemeshingPatch
    interploated_rem_patches = remesh_N(rem_patches_A, N=501)[end]
    for (rem_patch, rem_patch_A, rem_patch_B) in zip(interploated_rem_patches, rem_patches_A_div, rem_patches_B_div)
        # interpolate the coordinates of the remeshing patch
        if rem_patch.upper_coords ≠ nothing
            for (coord, coord_A, coord_B) in zip(rem_patch.upper_coords, rem_patch_A.upper_coords, rem_patch_B.upper_coords)
                coord.r = (1.0-x)*coord_A.r+x*coord_B.r
                coord.z = (1.0-x)*coord_A.z+x*coord_B.z
                coord.Φ₊ = (1.0-x)*coord_A.Φ₊+x*coord_B.Φ₊
            end
        end
        if rem_patch.mid_coords ≠ nothing
            for (coord, coord_A, coord_B) in zip(rem_patch.mid_coords, rem_patch_A.mid_coords, rem_patch_B.mid_coords)
                coord.r = (1.0-x)*coord_A.r+x*coord_B.r
                coord.z = (1.0-x)*coord_A.z+x*coord_B.z
                coord.Φ₊ = (1.0-x)*coord_A.Φ₊+x*coord_B.Φ₊
            end
        end
        if rem_patch.lower_coords ≠ nothing
            for (coord, coord_A, coord_B) in zip(rem_patch.lower_coords, rem_patch_A.lower_coords, rem_patch_B.lower_coords)
                coord.r = (1.0-x)*coord_A.r+x*coord_B.r
                coord.z = (1.0-x)*coord_A.z+x*coord_B.z
                coord.Φ₊ = (1.0-x)*coord_A.Φ₊+x*coord_B.Φ₊
            end
        end
    end

    # remesh the interploated patch to get something runnable
    coords, indep_coords, vol_patches, rem_patches = remesh(interploated_rem_patches)
    
    if vol_patches_A ≠ nothing && vol_patches_B ≠ nothing
        # interpolate the volume patches
        for (vol_patch, vol_patch_A, vol_patch_B) in zip(vol_patches, vol_patches_A, vol_patches_B)
            vol_patch.λᵥ = vol_patch_A.λᵥ
            vol_patch.λᵩ = vol_patch_A.λᵩ
        end
    end

    return coords, indep_coords, vol_patches, rem_patches
end

function linear_interpolate(rem_patches_A::Vector{RemeshingPatch}, rem_patches_B::Vector{RemeshingPatch}, M::Int64=20)
    coords_arr = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
    indep_coords_arr = Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}}[]
    vol_patches_arr = Vector{VolumePatch}[]
    rem_patches_arr = Vector{RemeshingPatch}[]
    for x in LinRange(0.0, 1.0, M)
        coords, indep_coords, vol_patches, rem_patches = _interpolate_helper(rem_patches_A, rem_patches_B, x)
        push!(coords_arr, coords)
        push!(indep_coords_arr, indep_coords)
        push!(vol_patches_arr, vol_patches)
        push!(rem_patches_arr, rem_patches)
    end
    return coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr
end

function interpolate(rem_patches_A::Vector{RemeshingPatch}, rem_patches_B::Vector{RemeshingPatch}, M::Int64=20, step::Int64=2, s::Float64=0.1, ; keep_edges::Bool=false, vol_patches_A::Union{Vector{VolumePatch}, Nothing}=nothing, vol_patches_B::Union{Vector{VolumePatch}, Nothing}=nothing)
    coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches_A)

    function force(x)
        x_to_coords!(x, coords, vol_patches, update_en=false)
        return gradient_parallel(indep_coords)
    end

    function set(x)
        x_to_coords!(x, coords, vol_patches, update_en=true)
    end

    coords_arr = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
    indep_coords_arr = Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}}[]
    vol_patches_arr = Vector{VolumePatch}[]
    rem_patches_arr = Vector{RemeshingPatch}[]

    if keep_edges
        push!(coords_arr, coords)
        push!(indep_coords_arr, indep_coords)
        push!(vol_patches_arr, vol_patches)
        push!(rem_patches_arr, rem_patches)
        coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches_A)
        range_ = 2:M-1
    else
        range_ = 1:M
    end

    μᵥ, μᵩ = 1.0e6, 1.0e6
    for i in range_
        for j in 1:step
            coords, indep_coords, vol_patches, rem_patches = _interpolate_helper(rem_patches, rem_patches_B, s, vol_patches_A = vol_patches_A, vol_patches_B = vol_patches_B)
            constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches)
            for k in 1:20
                x, _, _, _, _, _ = fire2(force, set, pos(coords), max_steps=50, tmin = 1.0e-4, tmax = 1.0e-2)
                _update_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, μᵥ, μᵩ, update_λs=true, verbose=false)
                if μᵥ > 1.0e4
                    μᵥ *= 0.99
                end
                if μᵩ > 1.0e4
                    μᵩ *= 0.99
                end
            end
        end
        println("Energy of current", energy(vol_patches))
        push!(coords_arr, coords)
        push!(indep_coords_arr, indep_coords)
        push!(vol_patches_arr, vol_patches)
        push!(rem_patches_arr, rem_patches)
    end

    if keep_edges
        coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches_B)
        push!(coords_arr, coords)
        push!(indep_coords_arr, indep_coords)
        push!(vol_patches_arr, vol_patches)
        push!(rem_patches_arr, rem_patches)
    end

    plot_reaction(rem_patches_arr)
    return coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr
end

function relaxation(rem_patches::Vector{RemeshingPatch}, ; remesh_steps::Int64 = 10, steps_per_remesh::Int64 = 1000, μᵥ::Float64 = 2.5e4, μᵩ::Float64 = 2.5e4)
    coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)
    constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches)
    _update_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, 0.0, 0.0)
    
    use_Φ = any([vol_patch.use_Φ for vol_patch in vol_patches])

    function force(x)
        x_to_coords!(x, coords, vol_patches, update_en=false)
        return gradient_parallel(indep_coords)
    end

    function force_Φ(x)
        x_to_coords!(x, coords, vol_patches, update_en=false)
        return gradient_Φ_parallel(indep_coords)
    end
    
    function set(x)
        x_to_coords!(x, coords, vol_patches, update_en=true)
    end

    if use_Φ
        x, _, _, _, _, _ = fire2(force_Φ, set, pos(coords), max_steps=100)
    end

    μᵥs, μᵩs = LinRange(0.0, μᵥ, remesh_steps), LinRange(0.0, μᵩ, remesh_steps)

    for i in 1:remesh_steps
        println("relaxation step $i of $remesh_steps")
        x, _, _, _, _, _ = fire2(force, set, pos(coords), max_steps=steps_per_remesh)
        println("Energy of current: ", energy(vol_patches))
        coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)
        constrain_v_vol_patches, constrain_Φ_vol_patches = _get_constrained_patches(vol_patches)
        _update_constraints(constrain_v_vol_patches, constrain_Φ_vol_patches, μᵥs[i], μᵩs[i], update_λs=false)
        if use_Φ
            x, _, _, _, _, _ = fire2(force_Φ, set, pos(coords), max_steps=100)
        end
    end

    return coords, indep_coords, vol_patches, rem_patches
end