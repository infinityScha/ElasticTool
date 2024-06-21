# finite-differentiation spacing
const Δx = 1.0e-5

# construct functions for gradient
function gradient(coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}})::Vector{Float64}
    grad = Float64[]
    for coord in coords
        if !coord.fixed_z
            prev_z = coord.z
            coord.z += Δx
            dE = calc_dE(coord)
            coord.z = prev_z
            push!(grad, dE)
        end
        if !coord.fixed_r
            prev_r = coord.r
            coord.r += Δx
            dE = calc_dE(coord)
            coord.r = prev_r
            push!(grad, dE)
        end
        if !coord.fixed_Φ₊
            prev_Φ₊ = coord.Φ₊
            coord.Φ₊ += Δx
            dE = calc_dE(coord)
            coord.Φ₊ = prev_Φ₊
            push!(grad, dE)
        end
    end
    return grad / Δx
end

function dg(coord::Coordinate{VolumePatch,VolumeElement{VolumePatch}})::Vector{Float64}
    dg_arr = Float64[]
    if !coord.fixed_z
        prev_z = coord.z
        coord.z += Δx
        dE0 = calc_dE(coord)
        coord.z = prev_z - Δx
        dE1 = calc_dE(coord)
        coord.z = prev_z
        push!(dg_arr, dE0 - dE1)
    end
    if !coord.fixed_r
        prev_r = coord.r
        coord.r += Δx
        dE0 = calc_dE(coord)
        coord.r = prev_r - Δx
        dE1 = calc_dE(coord)
        coord.r = prev_r
        push!(dg_arr, dE0 - dE1)
    end
    if !coord.fixed_Φ₊
        prev_Φ₊ = coord.Φ₊
        coord.Φ₊ += Δx
        dE0 = calc_dE(coord)
        coord.Φ₊ = prev_Φ₊ - Δx
        dE1 = calc_dE(coord)
        coord.Φ₊ = prev_Φ₊
        push!(dg_arr, dE0 - dE1)
    end
    return dg_arr
end

function gradient_parallel(indep_coords::Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}})::Vector{Float64}
    grad_ = Vector{Float64}[]
    for coords in indep_coords
        temp = Vector{Vector{Float64}}(undef, length(coords))
        @threads for i in 1:length(coords)
            temp[i] = dg(coords[i])
        end
        push!(grad_, vcat(temp...))
    end
    return vcat(grad_...) / Δx * 0.5
end

function dg_Φ(coord::Coordinate{VolumePatch,VolumeElement{VolumePatch}})::Vector{Float64}
    dg_arr = Float64[]
    if !coord.fixed_z
        push!(dg_arr, 0.0)
    end
    if !coord.fixed_r
        push!(dg_arr, 0.0)
    end
    if !coord.fixed_Φ₊
        prev_Φ₊ = coord.Φ₊
        coord.Φ₊ += Δx
        dE0 = calc_dE(coord)
        coord.Φ₊ = prev_Φ₊ - Δx
        dE1 = calc_dE(coord)
        coord.Φ₊ = prev_Φ₊
        push!(dg_arr, dE0 - dE1)
    end
    return dg_arr
end

function gradient_Φ_parallel(indep_coords::Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}})::Vector{Float64}
    grad_ = Vector{Float64}[]
    for coords in indep_coords
        temp = Vector{Vector{Float64}}(undef, length(coords))
        @threads for i in 1:length(coords)
            temp[i] = dg_Φ(coords[i])
        end
        push!(grad_, vcat(temp...))
    end
    return vcat(grad_...) / Δx * 0.5
end

function dg_neb(coord::Coordinate{VolumePatch,VolumeElement{VolumePatch}})::Tuple{Vector{Float64}, Vector{Float64}}
    dg_arr = Float64[]
    dg_Δsing_arr = Float64[]
    if !coord.fixed_z
        prev_z = coord.z
        coord.z += Δx
        dE0, dE_Δsing0 = calc_dE(coord, separate_Δsing = true)
        coord.z = prev_z - Δx
        dE1, dE_Δsing1 = calc_dE(coord, separate_Δsing = true)
        coord.z = prev_z
        push!(dg_arr, dE0 - dE1)
        push!(dg_Δsing_arr, dE_Δsing0 - dE_Δsing1)
    end
    if !coord.fixed_r
        prev_r = coord.r
        coord.r += Δx
        dE0, dE_Δsing0 = calc_dE(coord, separate_Δsing = true)
        coord.r = prev_r - Δx
        dE1, dE_Δsing1 = calc_dE(coord, separate_Δsing = true)
        coord.r = prev_r
        push!(dg_arr, dE0 - dE1)
        push!(dg_Δsing_arr, dE_Δsing0 - dE_Δsing1)
    end
    if !coord.fixed_Φ₊
        prev_Φ₊ = coord.Φ₊
        coord.Φ₊ += Δx
        dE0, dE_Δsing0 = calc_dE(coord, separate_Δsing = true)
        coord.Φ₊ = prev_Φ₊ - Δx
        dE1, dE_Δsing1 = calc_dE(coord, separate_Δsing = true)
        coord.Φ₊ = prev_Φ₊
        push!(dg_arr, dE0 - dE1)
        push!(dg_Δsing_arr, dE_Δsing0 - dE_Δsing1)
    end
    return (dg_arr, dg_Δsing_arr)
end

function gradient_neb_parallel(indep_coords::Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}})::Tuple{Vector{Float64}, Vector{Float64}}
    grad_ = Vector{Float64}[]
    grad_Δsing = Vector{Float64}[]
    for coords in indep_coords
        temp = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(coords))
        @threads for i in 1:length(coords)
            temp[i] = dg_neb(coords[i])
        end
        push!(grad_, vcat([t[1] for t in temp]...))
        push!(grad_Δsing, vcat([t[2] for t in temp]...))
    end
    return (vcat(grad_...) / Δx * 0.5, vcat(grad_Δsing...) / Δx * 0.5)
end