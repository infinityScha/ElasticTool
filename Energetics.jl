using QuadGK

# resolution-breaking penalty
const k_min = 5.0e3
const k_min2 = 1.0e3
# wall at R~0 for cylindrical symmetry
const R_min = 0.1
# auxiliary variable gradient penalty coefficiect to give the correct boundary width
const A = 2.0 * (DOMAINS_WIDTH/π) * (DOMAINS_WIDTH/π)
const LINE_TENSION_FIX = LINE_TENSION * 4.0 / DOMAINS_WIDTH

# Define Lagrange interpolation function
function lagrange_interpolation(s, f_values)
    x_values = LinRange(0.0, s, length(f_values))
    function lagrange_poly(x)
        result = 0.0
        for (i, xi) in enumerate(x_values)
            term = f_values[i]
            for (j, xj) in enumerate(x_values)
                if i ≠ j
                    term *= (x - xj) / (xi - xj)
                end
            end
            result += term
        end
        return result
    end
    return lagrange_poly
end

function derivative_of_lagrange_interpolation(s, f_values)
    x_values = LinRange(0.0, s, length(f_values))
    function lagrange_poly_derivative(x)
        result = 0.0
        for (i, xi) in enumerate(x_values)
            for (j, xj) in enumerate(x_values)
                if i == j
                    continue
                end
                term = f_values[i]
                for (k, xk) in enumerate(x_values)
                    if ((i ≠ k) & (j ≠ k))
                        term *= (x - xk) / (xi - xk)
                    end
                end
                term *= 1.0 / (xi - xj)
                result += term
            end
        end
        return result
    end
    return lagrange_poly_derivative
end

function norm_2(v)
    return dot(v, v)
end

function update!(vol_patch::VolumePatch)
    en = 0.0
    vol = 0.0
    area_ = 0.0
    Φ_ = 0.0
    
    @threads for vol_elem in vol_patch.vol_elems
        update!(vol_elem)
    end

    for vol_elem in vol_patch.vol_elems
        vol += vol_elem.volume
        area_ += vol_elem.energy[10]
        Φ_ += vol_elem.energy[11]*vol_elem.energy[10]
        en += vol_elem.energy[1]
    end
    vol_patch.total_energy = en
    vol_patch.total_volume = vol
    vol_patch.total_area = area_
    vol_patch.total_Φ = Φ_
    if vol_patch.constrain_v
        dv_ = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
        λᵥdv_ = vol_patch.λᵥ * dv_
        if λᵥdv_ < 0.0
            λᵥdv_ *= exp(-50.0*abs(dv_))
        end
        vol_patch.total_energy += λᵥdv_ + 0.5 * vol_patch.μᵥ * dv_ * dv_
    end
    if vol_patch.constrain_Φ
        dΦ_ = vol_patch.total_Φ/vol_patch.total_area - vol_patch.Φ₀
        λᵩdΦ_ = vol_patch.λᵩ * dΦ_
        if λᵩdΦ_ < 0.0
            λᵩdΦ_ *= exp(-50.0*abs(dΦ_))
        end
        vol_patch.total_energy += λᵩdΦ_ + 0.5 * vol_patch.μᵩ * dΦ_ * dΦ_
    end
end

function update!(vol_elem::VolumeElement{VolumePatch})
    vol_elem.volume = calc_volume(vol_elem)
    vol_elem.energy = calc_energy(vol_elem)
end

# define a function to calculate the cross product of two vectors
function cross_product(a::Coordinate{VolumePatch,VolumeElement{VolumePatch}}, b::Coordinate{VolumePatch,VolumeElement{VolumePatch}})
    return a.r * b.z - a.z * b.r
end

function cross_product(a::Coordinate{VolumePatch,VolumeElement{VolumePatch}}, b::Coordinate{VolumePatch,VolumeElement{VolumePatch}}, c::Coordinate{VolumePatch,VolumeElement{VolumePatch}})
    return (a.r - c.r) * (b.z - c.z) - (a.z - c.z) * (b.r - c.r)
end

# define the Jarvis March algorithm
function jarvis_march(coordinates::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}})
    # find the leftmost coordinate to start with
    p0 = coordinates[argmin([c.r for c in coordinates])]

    hull = [p0]
    while true
        endpoint = coordinates[1]
        for i in 2:4
            if endpoint.idx == hull[end].idx || cross_product(coordinates[i], endpoint, hull[end]) > 0
                endpoint = coordinates[i]
            end
        end

        if endpoint == p0
            break
        end

        push!(hull, endpoint)
    end

    return -calc_volume2d(hull), calc_volume2d(coordinates)
end

function calc_volume2d(coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}})
    n = length(coords)

    area = 0
    for i in 1:n
        j = mod1(i + 1, n)
        area += cross_product(coords[i], coords[j])
    end

    area /= 2

    return -area
end

function calc_volume(coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}})
    n = length(coords)

    # for water patches (VOLUME CAN BE NEGATIVE, the total volume would be positive)
    if n == 2
        @static if IS_PLANAR
            return -0.5*(coords[2].z - coords[1].z) * (coords[2].r + coords[1].r)
        else
            return -π*(coords[2].z - coords[1].z) * (coords[2].r * coords[2].r + coords[1].r * coords[2].r + coords[1].r * coords[1].r) / 3
        end
    end

    vol = 0.0
    for i in 1:n
        j = mod1(i + 1, n)
        @static if IS_PLANAR
            vol += (coords[j].z - coords[i].z) * (coords[j].r + coords[i].r)
        else
            vol += (coords[j].z - coords[i].z) * (coords[j].r * coords[j].r + coords[i].r * coords[j].r + coords[i].r * coords[i].r)
        end
    end

    @static if IS_PLANAR
        vol *= 0.5
    else
        vol *= π / 3.0
    end

    return -vol
end

function calc_volume(vol_elem::VolumeElement{VolumePatch})
    if vol_elem.vol_patch.constrain_s
        return sqrt((vol_elem.coords[2].r - vol_elem.coords[1].r)*(vol_elem.coords[2].r - vol_elem.coords[1].r) + (vol_elem.coords[2].z - vol_elem.coords[1].z)*(vol_elem.coords[2].z - vol_elem.coords[1].z))
    end
    # if there are extra midpoints...
    return length(vol_elem.coords) > 4 ? calc_volume(vol_elem.coords[1:4]) : calc_volume(vol_elem.coords)
end

function calc_energy(vol_elem::VolumeElement{VolumePatch}, ; only_Φ::Bool=false)
    coords = vol_elem.coords

    # check if the volume element is of a water patch, these most have 3 coordinates only
    if length(coords) < 4
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5
    end

    r1, r2, r3, r4 = (coords[i].r for i in 1:4)
    z1, z2, z3, z4 = (coords[i].z for i in 1:4)
    d0 = SVector{2,Float64}(r4 - r1, z4 - z1)
    d1 = SVector{2,Float64}(r3 - r2, z3 - z2)
    l0, l1 = norm(d0), norm(d1)
    d0 /= l0
    d1 /= l1
    d_ = 0.5 * (d0 + d1)
    dd_ = d1 - d0
    r_ = 0.5 * (r1 + r2)
    tau_ = SVector{2,Float64}(r2 - r1, z2 - z1)
    Δs = norm(tau_)
    tau_ /= Δs

    N_ = SVector{2,Float64}(tau_[2], -tau_[1])

    dd_dτ = dot(dd_, tau_) / Δs

    @static if IS_PLANAR
        divd = dd_dτ
        # TODO: implement the planar case
        det∇d = 0.0 
    else
        dᵣ_r = d_[1] / r_
        divd = dᵣ_r + dd_dτ
        det∇d = dot(dd_dτ, dᵣ_r)
    end

    # Precompute values
    κ, Jₛ, κᵍ, kₜ, Kₐ, h₀  = vol_elem.vol_patch.κ, vol_elem.vol_patch.Jₛ, vol_elem.vol_patch.κᵍ, vol_elem.vol_patch.kₜ, vol_elem.vol_patch.Kₐ, vol_elem.vol_patch.h₀
    
    Φing = 0.0
    Δsing = 0.0
    mean_Φ = 0.5
    if vol_elem.vol_patch.use_Φ
        Φ₊_arr = length(coords) > 4 ? [coord_.Φ₊ for coord_ in coords[5:end]] : Float64[]
        pushfirst!(Φ₊_arr, coords[1].Φ₊)
        push!(Φ₊_arr, coords[2].Φ₊)
        Φ_arr = [Φ(Φ₊_) for Φ₊_ in Φ₊_arr]
        mean_Φ = 0.5 * (Φ_arr[1] + Φ_arr[end])
        κ  = (1.0 - mean_Φ) * κ + mean_Φ * vol_elem.vol_patch.κᴵ
        Jₛ = (1.0 - mean_Φ) * Jₛ + mean_Φ * vol_elem.vol_patch.Jₛᴵ
        κᵍ = (1.0 - mean_Φ) * κᵍ + mean_Φ * vol_elem.vol_patch.κᵍᴵ
        kₜ = (1.0 - mean_Φ) * kₜ + mean_Φ * vol_elem.vol_patch.kₜᴵ
        Kₐ = (1.0 - mean_Φ) * Kₐ + mean_Φ * vol_elem.vol_patch.Kₐᴵ
        h₀ = (1.0 - mean_Φ) * h₀ + mean_Φ * vol_elem.vol_patch.h₀ᴵ
        
        for Φ₊ in Φ₊_arr
            if abs(Φ₊) > 9.21034036724
                ratio = 1 - abs(Φ₊)/9.21034036724
                Δsing += 0.25 * k_min * ratio * ratio
            end
        end
        
        @static if DOMAINS
            # gaussian quadrature for the integral of Φ(Φ-1)
            temp = 0.0
            Φ_interp = lagrange_interpolation(Δs, Φ_arr)
            ∇Φ_interp = derivative_of_lagrange_interpolation(Δs, Φ_arr)
            function en_Φ_interp(x)
                Φ = Φ_interp(x)
                ∇Φ = ∇Φ_interp(x)
                return Φ * (1.0 - Φ) + 0.5 * A * ∇Φ * ∇Φ
            end
            temp, err = quadgk(en_Φ_interp, 0.0, Δs, order=3)
            temp /= Δs
            Φing += DOMAINS_ENF_STR * temp
            @static if !IS_PLANAR
                Φing /= 2 * π * r_
            end
            @static if ADD_LINE_TENSION
                Φing += LINE_TENSION_FIX * temp
            end 
        else
            Φ₀ = vol_elem.vol_patch.Φ₀
            ∇Φ = (Φ_arr[end] - Φ_arr[1]) / Δs
            Φing = (mean_Φ*log(mean_Φ/Φ₀) + (1 - mean_Φ)*log((1 - mean_Φ)/(1-Φ₀)) + CHI*mean_Φ*(Φ₀-mean_Φ) + Kₗᵢₙₑ*∇Φ*∇Φ)*4.114/0.639 # APL for POPC from nanodisc paper.
            # normalize Φing to element's area deformation
            Φing /= 1.0 + 0.5*((h₀ - l0) / sqrt(l0*h₀) + (h₀ - l1) / sqrt(l1*h₀))
        end
    end

    @static if DOMAINS
        if only_Φ
            return Φing * Δs + Δsing, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5
        end
    end
    
    bending = 0.5 * κ * ((divd + Jₛ) * (divd + Jₛ) - Jₛ * Jₛ)
    # hkk tilt definition
    tilting = 0.5 * kₜ * norm_2(N_ - d_ / dot(N_, d_))
    # following Rytham stretch definition
    stretching = 0.25 * Kₐ * ((l0 - h₀) * (l0 - h₀) / l0 + (l1 - h₀) * (l1 - h₀) / l1) / h₀
    # gaussian energy
    g_curving = κᵍ*det∇d

    @static if IS_Z_SELF_INTERACTION
        if vol_elem.vol_patch.self_interactions && (N_[1] > 0.5)
            #selfing = - Φ_ * (0.5*tanh(100.0*(r_-1.4))+0.5) * ϵ_Z_SELF_INTERACTION * exp(- 2.0 * r_ / σ_Z_SELF_INTERACTION) * 4.0 * (N_[1] - 0.5) * (N_[1] - 0.5) # the tanh turns off the interaction near r=1.5 nm
            selfing = - mean_Φ * ϵ_Z_SELF_INTERACTION * exp(- 2.0 * r_ / σ_Z_SELF_INTERACTION) * 4.0 * (N_[1] - 0.5) * (N_[1] - 0.5)
        else
            selfing = 0.0
        end
    else
        selfing = 0.0
    end    

    if Δs < vol_elem.min_Δs
        ratio = 1 - Δs / vol_elem.min_Δs
        Δsing += 0.5 * k_min * ratio * ratio
    else
        if Δs > vol_elem.max_Δs
            ratio = 1 - Δs / vol_elem.max_Δs
            Δsing += 0.5 * k_min * ratio * ratio
        end
    end
    
    Δs_mid = sqrt((r3 - r4)*(r3 - r4) + (z3 - z4)*(z3 - z4))
    ratio = 1 - Δs_mid / vol_elem.Δs_mid
    Δsing += 0.25 * k_min2 * ratio * ratio

    if Δs_mid > 1.1 * vol_elem.Δs_mid
        ratio = 1 - Δs_mid / (1.1 * vol_elem.Δs_mid)
        Δsing += 0.25 * k_min * ratio * ratio
    end

    if Δs_mid < 0.9 * vol_elem.Δs_mid
        ratio = 1 - Δs_mid / (0.9 * vol_elem.Δs_mid)
        Δsing += 0.25 * k_min * ratio * ratio
    end

    
    @static if IS_UNSTABLE
        tau_mid = SVector{2,Float64}(vol_elem.coords[4].r - vol_elem.coords[3].r, vol_elem.coords[4].z - vol_elem.coords[3].z)
        tau_mid /= norm(tau_mid)
    
        cosθₛₘ = abs(dot(tau_, tau_mid))
    
        if cosθₛₘ < 0.877
            ratio = 1.0 - cosθₛₘ / 0.877
            Δsing += 0.5 * k_min * ratio * ratio * 10.0 
        end    
        
        # compare convex hull area with the area of the volume element
        # if they are different, then the volume element is not convex and thus unrealistic.
        ch_area, area = jarvis_march(vol_elem.coords[1:4])
        if ch_area != area
            ratio = 1 - area / ch_area
            Δsing += 0.5 * k_min * ratio * ratio
        end

        @static if !IS_PLANAR
            r_min = R_min * vol_elem.Δs_mid

            for i in [1,2]
                ri = vol_elem.coords[i].r
                if !vol_elem.coords[i].fixed_r && ri < r_min
                    ratio = 1 - ri / r_min
                    Δsing += 0.5 * k_min * ratio * ratio * exp(1.0/(1.0-ratio))
                end
            end
        end
    end

    @static if !IS_PLANAR
        ΔS = 2 * π * r_ * Δs
    else
        ΔS = Δs
    end

    return (bending + g_curving + tilting + stretching + Φing + selfing) * ΔS + Δsing, bending, g_curving, tilting, stretching, Φing, selfing, Δsing, Δs, ΔS, mean_Φ
end

function calc_dE(coord::Coordinate{VolumePatch,VolumeElement{VolumePatch}}, ; only_Φ::Bool= false, separate_Δsing = false)
    dE = 0.0
    dE_Δsing = 0.0
    for (vol_patch, vol_elems) in zip(coord.vol_patches, coord.vol_elems)
        dv = 0.0
        dΦ = 0.0
        dA = 0.0
        for vol_elem in vol_elems
            elem_en = calc_energy(vol_elem, only_Φ=only_Φ)
            dE += separate_Δsing ? elem_en[ENERGY_MAP["total"]] - elem_en[ENERGY_MAP["dsing"]] - vol_elem.energy[ENERGY_MAP["total"]] + vol_elem.energy[ENERGY_MAP["dsing"]] : elem_en[ENERGY_MAP["total"]] - vol_elem.energy[ENERGY_MAP["total"]]
            dE_Δsing += elem_en[ENERGY_MAP["dsing"]] - vol_elem.energy[ENERGY_MAP["dsing"]]
            dv += calc_volume(vol_elem) - vol_elem.volume
            dΦ += elem_en[ENERGY_MAP["mean_Φ"]]*elem_en[ENERGY_MAP["dS"]] - vol_elem.energy[ENERGY_MAP["mean_Φ"]]*vol_elem.energy[ENERGY_MAP["dS"]]
            dA += elem_en[ENERGY_MAP["dS"]] - vol_elem.energy[ENERGY_MAP["dS"]]
        end
        if vol_patch.constrain_Φ
            dΦ_ =  (vol_patch.total_Φ + dΦ)/(vol_patch.total_area + dA) - vol_patch.Φ₀
            λᵩdΦ_ = vol_patch.λᵩ * dΦ_
            if λᵩdΦ_ < 0.0
                λᵩdΦ_ *= exp(-50.0*abs(dΦ_))
            end
            dΦ_old =  (vol_patch.total_Φ)/(vol_patch.total_area) - vol_patch.Φ₀
            λᵩdΦ_old = vol_patch.λᵩ * dΦ_old
            if λᵩdΦ_old < 0.0
                λᵩdΦ_old *= exp(-50.0*abs(dΦ_old))
            end
            dE += λᵩdΦ_ - λᵩdΦ_old + 0.5 * vol_patch.μᵩ * (dΦ_ * dΦ_ - dΦ_old * dΦ_old)
        end
        if only_Φ
            continue
        end
        if vol_patch.constrain_v
            dv_ = (vol_patch.total_volume + dv - vol_patch.v₀) / vol_patch.v₀
            λᵥdv_ = vol_patch.λᵥ * dv_
            if λᵥdv_ < 0.0
                λᵥdv_ *= exp(-50.0*abs(dv_))
            end
            dv_old = (vol_patch.total_volume - vol_patch.v₀) / vol_patch.v₀
            λᵥdv_old = vol_patch.λᵥ * dv_old
            if λᵥdv_old < 0.0
                λᵥdv_old *= exp(-50.0*abs(dv_old))
            end
            dE += λᵥdv_ - λᵥdv_old + 0.5 * vol_patch.μᵥ * (dv_ * dv_ - dv_old * dv_old)
        end
    end
    return separate_Δsing ? (dE, dE_Δsing) : dE
end

function calc_energy(vol_patches::Vector{VolumePatch},which::String="total")
    @static if DOMAINS
        which_dict = Dict("total" => [ENERGY_MAP["bending"], ENERGY_MAP["g_curving"], ENERGY_MAP["tilting"], ENERGY_MAP["stretching"], ENERGY_MAP["selfing"]], "bending" => [ENERGY_MAP["bending"]], "g_curving" => [ENERGY_MAP["g_curving"]], "tilting" => [ENERGY_MAP["tilting"]], "stretching" => [ENERGY_MAP["stretching"]], "Φing" => [ENERGY_MAP["Φing"]], "selfing" => [ENERGY_MAP["selfing"]], "dsing" => [ENERGY_MAP["dsing"]])
    else
        which_dict = Dict("total" => [ENERGY_MAP["bending"], ENERGY_MAP["g_curving"], ENERGY_MAP["tilting"], ENERGY_MAP["stretching"], ENERGY_MAP["Φing"], ENERGY_MAP["selfing"]], "bending" => [ENERGY_MAP["bending"]], "g_curving" => [ENERGY_MAP["g_curving"]], "tilting" => [ENERGY_MAP["tilting"]], "stretching" => [ENERGY_MAP["stretching"]], "Φing" => [ENERGY_MAP["Φing"]], "selfing" => [ENERGY_MAP["selfing"]], "dsing" => [ENERGY_MAP["dsing"]])
    end
    en = 0.0
    for vol_patch in vol_patches
        for vol_elem in vol_patch.vol_elems
            if which ≠ "dsing"
                en += sum(vol_elem.energy[which_dict[which]]) * vol_elem.energy[10]
            else
                en += sum(vol_elem.energy[which_dict[which]])
            end
        end
    end
    return en
end

calc_energy(vol_patch::VolumePatch,which::String="total") = calc_energy([vol_patch], which)
