const ENERGY_MAP = Dict("total"=>1, "bending"=>2, "g_curving"=>3, "tilting"=>4, "stretching"=>5, "Φing"=>6, "selfing"=>7, "dsing"=>8, "ds"=>9, "dS"=>10, "mean_Φ"=>11)

counter_coordinate = 1
counter_volume_element = 1
counter_volume_patch = 1

function Φ(Φ₊::Float64)
    return 0.5*(tanh(Φ₊)+1.0)
end

function Φ₊(Φ::Float64)
    return atanh(2.0*(Φ-0.5))
end

abstract type AbsCoordinate end
abstract type AbsVolumeElement end
abstract type AbsVolumePatch end

mutable struct Coordinate{T<:AbsVolumePatch,P<:AbsVolumeElement} <: AbsCoordinate
    const idx::Int64
    r::Float64
    z::Float64
    Φ₊::Float64
    const vol_patches::Array{T,1}
    const vol_elems::Array{Array{P,1}}
    const fixed_r::Bool
    const fixed_z::Bool
    const fixed_Φ₊::Bool
end

mutable struct VolumeElement{T<:AbsVolumePatch} <: AbsVolumeElement
    const idx::Int64
    const coords::Array{Coordinate{T,VolumeElement{T}},1} # surface0, surface1, midplane1, midplane0
    const vol_patch::T
    const Δs_mid::Float64
    const max_Δs::Float64
    const min_Δs::Float64
    volume::Float64
    energy::NTuple{11,Float64}
end

mutable struct VolumePatch <: AbsVolumePatch
    const idx::Int64
    const vol_elems::Array{VolumeElement{VolumePatch},1}
    const constrain_v::Bool
    const v₀::Float64
    const constrain_s::Bool
    const use_Φ::Bool
    const constrain_Φ::Bool
    const Φ₀::Float64
    total_volume::Float64
    total_area::Float64
    total_Φ::Float64
    μᵥ::Float64
    λᵥ::Float64
    μᵩ::Float64
    λᵩ::Float64
    const κ::Float64
    const κᵍ::Float64
    const Jₛ::Float64
    const kₜ::Float64
    const Kₐ::Float64
    const h₀::Float64
    const κᴵ::Float64 
    const κᵍᴵ::Float64
    const Jₛᴵ::Float64
    const kₜᴵ::Float64
    const Kₐᴵ::Float64
    const h₀ᴵ::Float64
    const self_interactions::Bool
    total_energy::Float64
end

"""
Creates a volume patch for a leaflet without a Φ field over it.
The volume constraint has to be explicitly defined.
"""
function LeafletPatch(κ::Float64, κᵍ::Float64, Jₛ::Float64, kₜ::Float64, Kₐ::Float64, h₀::Float64, ; constrain_v::Bool=false, v₀::Float64=0.0, self_interactions::Bool=false)
    global counter_volume_patch
    obj = VolumePatch(
        counter_volume_patch,
        VolumeElement[],
        constrain_v,
        v₀,
        false,false,false,0.0,0.0,0.0,0.0,
        1.0e2, 0.0, # μᵥ, λᵥ
        1.0e6, 0.0, # μᵩ, λᵩ
        κ,κᵍ,Jₛ,kₜ,Kₐ,h₀,
        0.0,0.0,0.0,0.0,0.0,0.0,
        self_interactions,
        0.0
    )
    counter_volume_patch += 1
    return obj
end

"""
Creates a volume patch for a leaflet !with! a Φ field over it.
The volume and Φ constraints has to be explicitly defined.
"""
function LeafletPatch(κ::Float64, κᵍ::Float64, Jₛ::Float64, kₜ::Float64, Kₐ::Float64, h₀::Float64, κᴵ::Float64, Jₛᴵ::Float64, κᵍᴵ::Float64, kₜᴵ::Float64, Kₐᴵ::Float64, h₀ᴵ::Float64, ; constrain_v::Bool=false, v₀::Float64=0.0, constrain_Φ::Bool=false, Φ₀::Float64=0.0, self_interactions::Bool=false)
    global counter_volume_patch
    obj = VolumePatch(
        counter_volume_patch,
        VolumeElement[],
        constrain_v,
        v₀,
        false,
        true,
        constrain_Φ, 
        Φ₀,
        0.0,0.0,0.0,
        1.0e2, 0.0, # μᵥ, λᵥ
        1.0e6, 0.0, # μᵩ, λᵩ
        κ,Jₛ,κᵍ,kₜ,Kₐ,h₀,
        κᴵ,Jₛᴵ,κᵍᴵ,kₜᴵ,Kₐᴵ,h₀ᴵ,
        self_interactions,
        0.0
    )
    counter_volume_patch += 1
    return obj
end

"""
Creates a volumetrically constrained water volume patch.
"""
function WaterPatch(v₀::Float64)
    global counter_volume_patch
    obj = VolumePatch(
        counter_volume_patch,
        VolumeElement[],
        true,
        v₀,
        false,false,false,0.0,0.0,0.0,0.0,
        1.0e2, # μᵥ
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,false,0.0
    )
    counter_volume_patch += 1
    return obj
end

"""
Creates a length (s) constrained volume patch.
"""
function AreaPatch(s₀::Float64)
    global counter_volume_patch
    obj = VolumePatch(
        counter_volume_patch,
        VolumeElement[],
        true, s₀,
        true,
        false,false,0.0,0.0,0.0,0.0,
        1.0e2, # μᵥ
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,false,0.0
    )
    counter_volume_patch += 1
    return obj
end

function VolumePatchLike(vol_patch::VolumePatch)
    global counter_volume_patch
    obj = VolumePatch(
        vol_patch.idx,
        VolumeElement[],
        vol_patch.constrain_v,
        vol_patch.v₀,
        vol_patch.constrain_s,
        vol_patch.use_Φ,
        vol_patch.constrain_Φ, 
        vol_patch.Φ₀,
        0.0,0.0,0.0,
        vol_patch.μᵥ,vol_patch.λᵥ,vol_patch.μᵩ,vol_patch.λᵩ,
        vol_patch.κ,vol_patch.κᵍ,vol_patch.Jₛ,vol_patch.kₜ,vol_patch.Kₐ,vol_patch.h₀,
        vol_patch.κᴵ,vol_patch.κᵍᴵ,vol_patch.Jₛᴵ,vol_patch.kₜᴵ,vol_patch.Kₐᴵ,vol_patch.h₀ᴵ,
        vol_patch.self_interactions,
        0.0
    )
    counter_volume_patch += 1
    return obj
end

function VolumeElement(coords::Array{Coordinate{VolumePatch,VolumeElement{VolumePatch}},1}, vol_patch::VolumePatch, Δs_mid::Float64)::VolumeElement{VolumePatch}
    global counter_volume_element
    obj = VolumeElement{VolumePatch}(
        counter_volume_element,
        coords,
        vol_patch,
        Δs_mid,
        Δs_mid * 2.5,
        Δs_mid / 2.5,
        0.0,
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    for coord in coords
        i = find_same_idx(coord.vol_patches, vol_patch)
        if i != -1
            push!(coord.vol_elems[i], obj)
            continue
        end
        push!(coord.vol_patches, vol_patch)
        push!(coord.vol_elems, [obj])
    end
    push!(vol_patch.vol_elems, obj)
    counter_volume_element += 1
    return obj
end

function VolumeElement(coords::Array{Coordinate{VolumePatch,VolumeElement{VolumePatch}},1}, vol_patch::VolumePatch)::VolumeElement{VolumePatch}
    global counter_volume_element
    if length(coords) > 2
        dr = coords[3].r - coords[4].r
        dz = coords[3].z - coords[4].z
        Δs_mid = sqrt(dr * dr + dz * dz)
    else
        Δs_mid = 0.0
    end
    return VolumeElement(coords, vol_patch, Δs_mid)
end

"""
    Creates a surface coordinate with (r,z,Φ) for cylindrical symmetry and (x,z,Φ) for planar symmetry.
    Fixed dimensions has to be set explicitly.
    Note: In most cases, fixed dimension are kept only for the coordinates in the edge of the remeshing patches. Refer to RemeshingPatch documentation to see exceptions to this rule.
"""
function SurfaceCoordinate(r::Float64, z::Float64, Φ::Float64, ; fixed_r::Bool=false, fixed_z::Bool=false, fixed_Φ::Bool=false)::Coordinate{VolumePatch,VolumeElement{VolumePatch}}
    global counter_coordinate
    obj = Coordinate{VolumePatch,VolumeElement{VolumePatch}}(
        counter_coordinate,
        r,
        z,
        Φ₊(Φ),
        VolumePatch[],
        Array{VolumeElement{VolumePatch}}[],
        fixed_r,
        fixed_z,
        fixed_Φ,
    )
    counter_coordinate += 1
    return obj
end

function SurfaceCoordinate(r::Float64, z::Float64, ; fixed_r::Bool=false, fixed_z::Bool=false)::Coordinate{VolumePatch,VolumeElement{VolumePatch}}
    global counter_coordinate
    obj = Coordinate{VolumePatch,VolumeElement{VolumePatch}}(
        counter_coordinate,
        r,
        z,
        0.0,
        VolumePatch[],
        Array{VolumeElement{VolumePatch}}[],
        fixed_r,
        fixed_z,
        true,
    )
    counter_coordinate += 1
    return obj
end

"""
    Creates a surface coordinate which only has a Φ field, used for composition interpolation within volume elements.
"""
function Φ_Coordinate(Φ::Float64)::Coordinate{VolumePatch,VolumeElement{VolumePatch}}
    global counter_coordinate
    obj = Coordinate{VolumePatch,VolumeElement{VolumePatch}}(
        counter_coordinate,
        0.0,
        0.0,
        Φ₊(Φ),
        VolumePatch[],
        Array{VolumeElement{VolumePatch}}[],
        true,
        true,
        false,
    )
    counter_coordinate += 1
    return obj
end

"""
    Creates a midplane coordinate with (r,z,Φ) for cylindrical symmetry and (x,z,Φ) for planar symmetry.
    Fixed dimensions has to be set explicitly.
    Note: In most cases, fixed dimension are kept only for the coordinates in the edge of the remeshing patches. Refer to RemeshingPatch documentation to see exceptions to this rule.
"""
function MidplaneCoordinate(r::Float64, z::Float64, ; fixed_r::Bool=false, fixed_z::Bool=false)::Coordinate{VolumePatch,VolumeElement{VolumePatch}}
    global counter_coordinate
    obj = Coordinate{VolumePatch,VolumeElement{VolumePatch}}(
        counter_coordinate,
        r,
        z,
        0.0,
        VolumePatch[],
        Array{VolumeElement{VolumePatch}}[],
        fixed_r,
        fixed_z,
        true,
    )
    counter_coordinate += 1
    return obj
end

function CoordinateLike(coord::Coordinate{VolumePatch,VolumeElement{VolumePatch}})::Coordinate{VolumePatch,VolumeElement{VolumePatch}}
    global counter_coordinate
    obj = Coordinate{VolumePatch,VolumeElement{VolumePatch}}(
        counter_coordinate,
        coord.r,
        coord.z,
        coord.Φ₊,
        VolumePatch[],
        Array{VolumeElement{VolumePatch}}[],
        coord.fixed_r,
        coord.fixed_z,
        coord.fixed_Φ₊,
    )
    counter_coordinate += 1
    return obj
end

function find_same_idx(arr::Array{T,1}, indexed_element::T) where {T}
    for i in 1:length(arr)
        if arr[i].idx == indexed_element.idx
            return i
        end
    end
    return -1
end

function omit_duplicates!(arr::Array{Union{T,Nothing},1}) where {T}
    to_omit = Int64[]
    for i in 1:length(arr)
        if arr[i] === nothing
            push!(to_omit, i)
        end
        if i in to_omit
            continue
        end
        for j in i+1:length(arr)
            if arr[j] === nothing
                continue
            end
            if arr[i].idx == arr[j].idx
                push!(to_omit, j)
            end
        end
    end
    deleteat!(arr, sort!(unique!(to_omit)))
end

function omit_duplicates!(arr::Array{T,1}) where {T}
    to_omit = Int64[]
    for i in 1:length(arr)
        if i in to_omit
            continue
        end
        for j in i+1:length(arr)
            if arr[i].idx == arr[j].idx
                push!(to_omit, j)
            end
        end
    end
    deleteat!(arr, sort!(unique!(to_omit)))
end

function shared_only(arr0::Array{T,1}, arr1::Array{T,1}) where {T}
    arr = T[]
    for elem in arr0
        idx = find_same_idx(arr1, elem)
        if idx == -1
            continue
        end
        push!(arr, elem)
    end
    return arr
end

function x_to_coords!(x::Array{Float64,1}, coords::Array{Coordinate{VolumePatch,VolumeElement{VolumePatch}},1}, vol_patches::Array{VolumePatch,1}, ; update_en::Bool=true)
    counter_ = 1
    for coord in coords
        if !coord.fixed_z
            coord.z = x[counter_]
            counter_ += 1
        end
        if !coord.fixed_r
            coord.r = x[counter_]
            counter_ += 1
        end
        if !coord.fixed_Φ₊
            coord.Φ₊ = x[counter_]
            counter_ += 1
        end
    end

    if !update_en
        return
    end

    for vol_patch in vol_patches
        update!(vol_patch)
    end
end

function pos(coords::Array{Coordinate{VolumePatch,VolumeElement{VolumePatch}},1})
    pos_arr = Float64[]
    for coord in coords
        if !coord.fixed_z
            push!(pos_arr, coord.z)
        end
        if !coord.fixed_r
            push!(pos_arr, coord.r)
        end
        if !coord.fixed_Φ₊
            push!(pos_arr, coord.Φ₊)
        end
    end
    return pos_arr
end