using Dierckx

"""
    RemeshingPatch(upper_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},
                   mid_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},
                   lower_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},
                   resolution::Float64, upper_vol_patch::VolumePatch, lower_vol_patch::VolumePatch, 
                   upper_water_patch::VolumePatch=Nothing, lower_water_patch::VolumePatch=Nothing,
                   fixed_midplane_r::Bool=false, fixed_midplane_z::Bool=false, midplane_area_patch::Union{VolumePatch,Nothing}=Nothing)

Remeshing patch defines the set of coordinates designated for remeshing in between. It also creates the volume elements and everything so it makes life way easier when given to the 'remesh' fuction

# Arguments
- `upper_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}`: An optional vector of `Coordinate` objects representing the upper coordinates of the patch (must be provided if the lower coordinates are not provided). These coordinates should be ordered such that the volume elements created in their direction have a positive orientation.
- `mid_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}`: A vector of `Coordinate` objects representing the mid coordinates of the patch. These coordinates should be in the same direction as the upper coordinates.
- `lower_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}`: An optional vector of `Coordinate` objects representing the lower coordinates of the patch (must be provided if the upper coordinates are not provided). These coordinates should be in the same direction as the upper coordinates.
- `resolution::Float64`: A floating-point number specifying the desired resolution of the midplane segments.
- `upper_vol_patch::VolumePatch`: The `VolumePatch` object representing the upper volume patch, must be provided if the upper coordinates are provided.
- `lower_vol_patch::VolumePatch`: The `VolumePatch` object representing the lower volume patch, must be provided if the lower coordinates are provided.
- `upper_water_patch::VolumePatch=Nothing`: An optional `VolumePatch` object representing the water volume patch above the upper coordinates.
- `lower_water_patch::VolumePatch=Nothing`: An optional `VolumePatch` object representing the water volume patch below the lower coordinates.
- `fixed_midplane_r::Bool`: A boolean specifying if the midplane r-coordinates are fixed.
- `fixed_midplane_z::Bool`: A boolean specifying if the midplane z-coordinates are fixed.
- `midplane_area_patch::Union{VolumePatch,Nothing}` : An optional `VolumePatch` object representing the midplane area patch such that the midplane area could be controlled.
- `upper_coord_Φ_order::Int64=0`: An integer specifying the order of the interpolation of the upper coordinates in Φ. Default value is 0.
- `lower_coord_Φ_order::Int64=0`: An integer specifying the order of the interpolation of the lower coordinates in Φ. Default value is 0.
"""
struct RemeshingPatch
    upper_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}
    mid_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}
    lower_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}
    resolution::Float64
    upper_vol_patch::Union{VolumePatch,Nothing}
    lower_vol_patch::Union{VolumePatch,Nothing}
    upper_water_patch::Union{VolumePatch,Nothing}
    lower_water_patch::Union{VolumePatch,Nothing}
    fixed_midplane_r::Bool
    fixed_midplane_z::Bool
    midplane_area_patch::Union{VolumePatch,Nothing}
    upper_coord_Φ_order::Int64
    lower_coord_Φ_order::Int64
    function RemeshingPatch(upper_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}=nothing,
        mid_coords::Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}=Coordinate{VolumePatch,VolumeElement}[],
        lower_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}=nothing,
        resolution::Float64=0.1,
        upper_vol_patch::Union{VolumePatch,Nothing}=nothing,
        lower_vol_patch::Union{VolumePatch,Nothing}=nothing,
        upper_water_patch::Union{VolumePatch,Nothing}=nothing,
        lower_water_patch::Union{VolumePatch,Nothing}=nothing,
        ; fixed_midplane_r::Bool=false,
        fixed_midplane_z::Bool=false,
        midplane_area_patch::Union{VolumePatch,Nothing}=nothing,
        upper_coord_Φ_order::Int64=0,
        lower_coord_Φ_order::Int64=0)

        new(upper_coords, mid_coords, lower_coords, resolution, upper_vol_patch, lower_vol_patch, upper_water_patch, lower_water_patch, fixed_midplane_r, fixed_midplane_z, midplane_area_patch, upper_coord_Φ_order, lower_coord_Φ_order)
    end
end

function change_rem_patch_resolution(rem_patch::RemeshingPatch, resolution::Float64)
    return RemeshingPatch(rem_patch.upper_coords, rem_patch.mid_coords, rem_patch.lower_coords, resolution, rem_patch.upper_vol_patch, rem_patch.lower_vol_patch, rem_patch.upper_water_patch, rem_patch.lower_water_patch, fixed_midplane_r = rem_patch.fixed_midplane_r, fixed_midplane_z = rem_patch.fixed_midplane_z, midplane_area_patch = rem_patch.midplane_area_patch, upper_coord_Φ_order = rem_patch.upper_coord_Φ_order, lower_coord_Φ_order = rem_patch.lower_coord_Φ_order)
end

function _fix_Φ(Φ::Float64)
    Φ = max(min(Φ, 0.99999999), 0.00000001)
    if Φ > 0.999
        Φ = 0.99999999
    elseif Φ < 0.001
        Φ = 0.00000001
    end
    return Φ
end

function _Φ(Φ₊::Float64)
    Φ = 0.5*(tanh(Φ₊)+1.0)
    return _fix_Φ(Φ)
end

function _Φ₊(Φ::Float64)
    Φ₊ = atanh(2.0*(_fix_Φ(Φ)-0.5))
    return Φ₊
end

"""
Remeshes a set of RemeshingPatches using a specified smoothing parameter.

## Arguments
- `rem_patches::Vector{RemeshingPatch}`: A vector of `RemeshingPatch` objects representing the patches to be remeshed.
- `s::Float64 = 0.0`: A floating-point value representing the smoothing parameter. Default value is 0.0.
- 'N::Vector{Int64} = Int64[]': The number of points for each of the remeshed patches. Used instead of the mesh resolution if set to any value other than -1.

## Returns
- A vector of `Coordinate` objects representing the remeshed coordinates.
- A list of vectors of `Coordinate` objects of the remeshed coordinates separated to subsets of independent coordinates for parallelized gradient calculation.
- A vector of `VolumePatch` objects representing the remeshed volume patches.
- A vector of `RemeshingPatch` objects representing the remeshed patches.
"""
function remesh(rem_patches::Vector{RemeshingPatch}, s::Float64=0.0, ; N::Vector{Int64}=Int64[])
    global counter_coordinate, counter_volume_element, counter_volume_patch

    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # create new coords, elements, patches
    new_coords_arr = Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]
    
    # create the endpoints of the remeshing parts first!
    edge_coords = Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]
    for rem_patch in rem_patches
        push!(edge_coords, rem_patch.mid_coords[1])
        push!(edge_coords, rem_patch.mid_coords[end])
        if rem_patch.upper_coords ≠ nothing
            push!(edge_coords, rem_patch.upper_coords[1])
            push!(edge_coords, rem_patch.upper_coords[end])
        end
        if rem_patch.lower_coords ≠ nothing
            push!(edge_coords, rem_patch.lower_coords[1])
            push!(edge_coords, rem_patch.lower_coords[end])
        end
    end

    omit_duplicates!(edge_coords)
    # create the coords as dict 
    # Note: it allegedly creates a surface coordinate here, THIS IS OK!!!!
    new_edge_coords = Dict(coord.idx => SurfaceCoordinate(coord.r, coord.z, coord.fixed_Φ₊ ? Φ(coord.Φ₊) : _Φ(coord.Φ₊), fixed_r=coord.fixed_r, fixed_z=coord.fixed_z, fixed_Φ=coord.fixed_Φ₊) for coord in edge_coords)
    new_coords_arr = vcat(new_coords_arr, collect(values(new_edge_coords)))
    
    # define new volume patches by dicts
    vol_patches = vcat([[rem_patch.upper_vol_patch, rem_patch.lower_vol_patch] for rem_patch in rem_patches]...)
    omit_duplicates!(vol_patches)
    new_vol_patches = Dict(vol_patch.idx => VolumePatchLike(vol_patch) for vol_patch in vol_patches)

    # create the water patches
    water_patches = VolumePatch[]
    for rem_patch in rem_patches
        if rem_patch.upper_water_patch ≠ nothing
            push!(water_patches, rem_patch.upper_water_patch)
        end
        if rem_patch.lower_water_patch ≠ nothing
            push!(water_patches, rem_patch.lower_water_patch)
        end
    end
    omit_duplicates!(water_patches)
    new_water_patches = Dict(water_patch.idx => VolumePatchLike(water_patch) for water_patch in water_patches)

    # create the midplane area patches
    midplane_area_patches = VolumePatch[]
    for rem_patch in rem_patches
        if rem_patch.midplane_area_patch ≠ nothing
            push!(midplane_area_patches, rem_patch.midplane_area_patch)
        end
    end
    omit_duplicates!(midplane_area_patches)
    new_midplane_area_patches = Dict(midplane_area_patch.idx => VolumePatchLike(midplane_area_patch) for midplane_area_patch in midplane_area_patches)

    indep_sur_coords_arr = [Coordinate{VolumePatch,VolumeElement{VolumePatch}}[], Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]]
    indep_mid_coords_arr = [Coordinate{VolumePatch,VolumeElement{VolumePatch}}[], Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]]
    indep_upper_Φ_coords_arr = [Coordinate{VolumePatch,VolumeElement{VolumePatch}}[] for _ in 1:maximum([rem_patch.upper_coord_Φ_order for rem_patch in rem_patches])]
    indep_lower_Φ_coords_arr = [Coordinate{VolumePatch,VolumeElement{VolumePatch}}[] for _ in 1:maximum([rem_patch.lower_coord_Φ_order for rem_patch in rem_patches])]

     use_N  = length(N) > 0

    # remesh each patch
    new_rem_patches = RemeshingPatch[]
    for (j, rem_patch) in enumerate(rem_patches)
        # extract the original 
        ds = [sqrt((coord0.r - coord1.r) * (coord0.r - coord1.r) + (coord0.z - coord1.z) * (coord0.z - coord1.z)) for (coord0, coord1) in zip(rem_patch.mid_coords[1:end-1], rem_patch.mid_coords[2:end])]
        pushfirst!(ds, 0.0)
        ds = cumsum(ds)
        mid_rs = [coord.r for coord in rem_patch.mid_coords]
        mid_zs = [coord.z for coord in rem_patch.mid_coords]
        # create the cubic spline interpolation of the upper/mid/lower coords.
        spl_mid_rs = Spline1D(ds, mid_rs, s=s)
        spl_mid_zs = Spline1D(ds, mid_zs, s=s)
        # choose the number of coordinates to be odd, this is needed for the parallelization
        N_ = (use_N && N[j] ≠ -1) ? N[j] : round(Int, ds[end] / rem_patch.resolution)
        N_ = iseven(N_) ? N_+1 : N_
        # create the coords
        new_ds = LinRange(0.0, ds[end], N_)[2:end-1]
        new_mid_coords = [MidplaneCoordinate(spl_mid_rs(ds_), spl_mid_zs(ds_), fixed_r=rem_patch.fixed_midplane_r, fixed_z=rem_patch.fixed_midplane_z) for ds_ in new_ds]
        new_coords_arr = vcat(new_coords_arr, new_mid_coords)
        pushfirst!(new_mid_coords, new_edge_coords[rem_patch.mid_coords[1].idx])
        push!(new_mid_coords, new_edge_coords[rem_patch.mid_coords[end].idx])

        # divide the mid_coords to 2 independent sets and add them to the currect subsets
        indep_mid_coords_arr[1] = vcat(indep_mid_coords_arr[1], new_mid_coords[1:2:end])
        indep_mid_coords_arr[2] = vcat(indep_mid_coords_arr[2], new_mid_coords[2:2:end])    

        new_upper_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing} = nothing
        new_upper_vol_patch::Union{VolumePatch,Nothing} = nothing
        if rem_patch.upper_coords ≠ nothing
            upper_rs = [coord.r for coord in rem_patch.upper_coords]
            upper_zs = [coord.z for coord in rem_patch.upper_coords]
            upper_Φ₊s = [coord.Φ₊ for coord in rem_patch.upper_coords]
            spl_upper_rs = Spline1D(ds, upper_rs, s=s)
            spl_upper_zs = Spline1D(ds, upper_zs, s=s)
            spl_upper_Φs = Spline1D(ds, _Φ.(upper_Φ₊s), s=s)
            new_upper_coords = [SurfaceCoordinate(spl_upper_rs(ds_), spl_upper_zs(ds_), _fix_Φ(spl_upper_Φs(ds_)), fixed_Φ=!rem_patch.upper_vol_patch.use_Φ) for ds_ in new_ds]
            new_coords_arr = vcat(new_coords_arr, new_upper_coords)
            pushfirst!(new_upper_coords, new_edge_coords[rem_patch.upper_coords[1].idx])
            push!(new_upper_coords, new_edge_coords[rem_patch.upper_coords[end].idx])

            # divide the upper_coords to 2 independent sets and add them to the currect subsets
            indep_sur_coords_arr[1] = vcat(indep_sur_coords_arr[1], new_upper_coords[1:2:end])
            indep_sur_coords_arr[2] = vcat(indep_sur_coords_arr[2], new_upper_coords[2:2:end])

            # create the volume patch
            new_upper_vol_patch = new_vol_patches[rem_patch.upper_vol_patch.idx]

            if rem_patch.upper_coord_Φ_order > 0
                new_upper_coords_Φ = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
                relevant_s = LinRange(0.0, ds[end], N_)
                for (s0, s1) in zip(relevant_s[1:end-1], relevant_s[2:end])
                    new_upper_coords_Φ_ = Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]
                    for i in 1:rem_patch.upper_coord_Φ_order
                        s_ = s0 + (s1 - s0) * i/(rem_patch.upper_coord_Φ_order+1)
                        new_coord = Φ_Coordinate(_fix_Φ(spl_upper_Φs(s_)))
                        push!(new_upper_coords_Φ_, new_coord)
                        push!(indep_upper_Φ_coords_arr[i], new_coord)
                    end
                    push!(new_upper_coords_Φ, new_upper_coords_Φ_)
                end
                [VolumeElement(vcat([c0, c1, c2, c3], c4_arr), new_upper_vol_patch) for (c0, c1, c2, c3, c4_arr) in zip(new_upper_coords[1:end-1], new_upper_coords[2:end], new_mid_coords[2:end], new_mid_coords[1:end-1], new_upper_coords_Φ)]
            else
                # create the volume elements and patch them together
                # note that the orientation is positive!
                 [VolumeElement([c0, c1, c2, c3], new_upper_vol_patch) for (c0, c1, c2, c3) in zip(new_upper_coords[1:end-1], new_upper_coords[2:end], new_mid_coords[2:end], new_mid_coords[1:end-1])]
            end
        end

        new_lower_coords::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing} = nothing
        new_lower_vol_patch::Union{VolumePatch,Nothing} = nothing
        if rem_patch.lower_coords ≠ nothing
            lower_rs = [coord.r for coord in rem_patch.lower_coords]
            lower_zs = [coord.z for coord in rem_patch.lower_coords]
            lower_Φ₊s = [coord.Φ₊ for coord in rem_patch.lower_coords]
            spl_lower_rs = Spline1D(ds, lower_rs, s=s)
            spl_lower_zs = Spline1D(ds, lower_zs, s=s)
            spl_lower_Φs = Spline1D(ds, _Φ.(lower_Φ₊s), s=s)
            new_lower_coords = [SurfaceCoordinate(spl_lower_rs(ds_), spl_lower_zs(ds_), _fix_Φ(spl_lower_Φs(ds_)), fixed_Φ=!rem_patch.lower_vol_patch.use_Φ) for ds_ in new_ds]
            new_coords_arr = vcat(new_coords_arr, new_lower_coords)
            pushfirst!(new_lower_coords, new_edge_coords[rem_patch.lower_coords[1].idx])
            push!(new_lower_coords, new_edge_coords[rem_patch.lower_coords[end].idx])
            
            # divide the upper_coords to 2 independent sets and add them to the currect subsets
            indep_sur_coords_arr[1] = vcat(indep_sur_coords_arr[1], new_lower_coords[1:2:end])
            indep_sur_coords_arr[2] = vcat(indep_sur_coords_arr[2], new_lower_coords[2:2:end])
                        
            new_lower_vol_patch = new_vol_patches[rem_patch.lower_vol_patch.idx]

            if rem_patch.lower_coord_Φ_order > 0
                new_lower_coords_Φ = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
                relevant_s = LinRange(0.0, ds[end], N_)
                for (s0, s1) in zip(relevant_s[1:end-1], relevant_s[2:end])
                    new_lower_coords_Φ_ = Coordinate{VolumePatch,VolumeElement{VolumePatch}}[]
                    for i in 1:rem_patch.lower_coord_Φ_order
                        s_ = s0 + (s1 - s0) * i/(rem_patch.lower_coord_Φ_order+1)
                        new_coord = Φ_Coordinate(_fix_Φ(spl_lower_Φs(s_)))
                        push!(new_lower_coords_Φ_, new_coord)
                        push!(indep_lower_Φ_coords_arr[i], new_coord)
                    end
                    push!(new_lower_coords_Φ, new_lower_coords_Φ_)
                end
                # print the type of new_lower_coords_Φ
                [VolumeElement(vcat([c0, c1, c2, c3], c4_arr), new_lower_vol_patch) for (c0, c1, c2, c3, c4_arr) in zip(new_lower_coords[2:end], new_lower_coords[1:end-1], new_mid_coords[1:end-1], new_mid_coords[2:end], new_lower_coords_Φ)]
            else
                # create the volume elements and patch them together
                # note that the orientation is positive!
                [VolumeElement([c0, c1, c2, c3], new_lower_vol_patch) for (c0, c1, c2, c3) in zip(new_lower_coords[2:end], new_lower_coords[1:end-1], new_mid_coords[1:end-1], new_mid_coords[2:end])]
            end
        end
        
        # remove duplicates for the independent coordinates
        for i in [1, 2]
            omit_duplicates!(indep_mid_coords_arr[i])
            omit_duplicates!(indep_sur_coords_arr[i])
        end

        # make 2-element volume elements for water patches if needed
        new_upper_water_patch::Union{VolumePatch,Nothing} = nothing
        if rem_patch.upper_water_patch ≠ nothing
            new_upper_water_patch = new_water_patches[rem_patch.upper_water_patch.idx]
            [VolumeElement([c0, c1], new_upper_water_patch) for (c0, c1) in zip(new_upper_coords[2:end], new_upper_coords[1:end-1])]
        end
        new_lower_water_patch::Union{VolumePatch,Nothing} = nothing
        if rem_patch.lower_water_patch ≠ nothing
            new_lower_water_patch = new_water_patches[rem_patch.lower_water_patch.idx]
            [VolumeElement([c0, c1], new_lower_water_patch) for (c0, c1) in zip(new_lower_coords[1:end-1], new_lower_coords[2:end])]
        end

        # add 2-element volume elements for water patches if needed for midplane if needed
        new_midplane_area_patch::Union{VolumePatch,Nothing} = nothing
        if rem_patch.midplane_area_patch ≠ nothing
            new_midplane_area_patch = new_midplane_area_patches[rem_patch.midplane_area_patch.idx]
            [VolumeElement([c0, c1], new_midplane_area_patch) for (c0, c1) in zip(new_mid_coords[1:end-1], new_mid_coords[2:end])]
        end

        # finally, create the new remeshing patch and add it to the remeshing patches
        push!(new_rem_patches, RemeshingPatch(new_upper_coords, new_mid_coords, new_lower_coords, rem_patch.resolution, new_upper_vol_patch, new_lower_vol_patch, new_upper_water_patch, new_lower_water_patch, fixed_midplane_r = rem_patch.fixed_midplane_r, fixed_midplane_z = rem_patch.fixed_midplane_z, midplane_area_patch = new_midplane_area_patch, upper_coord_Φ_order = rem_patch.upper_coord_Φ_order, lower_coord_Φ_order = rem_patch.lower_coord_Φ_order))
    end

    # add all independent subsets of coords to an array
    indep_coords_arr = vcat(indep_mid_coords_arr, indep_sur_coords_arr, indep_upper_Φ_coords_arr, indep_lower_Φ_coords_arr)
    # set new_coords_arr in the same order as indep_coords_arr
    new_coords_arr = vcat(indep_coords_arr...)

    # update the volume patches now that they are ready!
    new_vol_patches = vcat(collect(values(new_vol_patches)), collect(values(new_water_patches)), collect(values(new_midplane_area_patches)))
    for vol_patch in new_vol_patches
        update!(vol_patch)
    end

    return new_coords_arr, indep_coords_arr, new_vol_patches, new_rem_patches
end

function remesh_N(rem_patches::Vector{RemeshingPatch}, s::Float64=0.0, ; N::Int64=-1)
    return remesh(rem_patches, s, N=[N for _ in 1:length(rem_patches)])
end

function split_rem_patch(rem_patches::Vector{RemeshingPatch}, which::Int64, at::Int64)
    rem_patch = rem_patches[which]
    upper_A::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}, upper_B::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing} = nothing, nothing
    if rem_patch.upper_coords ≠ nothing
        upper_A, upper_B = rem_patch.upper_coords[1:at], rem_patch.upper_coords[at:end]
    end
    lower_A::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing}, lower_B::Union{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}},Nothing} = nothing, nothing
    if rem_patch.lower_coords ≠ nothing
        lower_A, lower_B = rem_patch.lower_coords[1:at], rem_patch.lower_coords[at:end]
    end
    mid_A, mid_B = rem_patch.mid_coords[1:at], rem_patch.mid_coords[at:end]
    rem_patch_A = RemeshingPatch(upper_A, mid_A, lower_A, rem_patch.resolution, rem_patch.upper_vol_patch, rem_patch.lower_vol_patch, rem_patch.upper_water_patch, rem_patch.lower_water_patch, fixed_midplane_r = rem_patch.fixed_midplane_r, fixed_midplane_z = rem_patch.fixed_midplane_z, midplane_area_patch = rem_patch.midplane_area_patch, upper_coord_Φ_order = rem_patch.upper_coord_Φ_order, lower_coord_Φ_order = rem_patch.lower_coord_Φ_order)
    rem_patch_B = RemeshingPatch(upper_B, mid_B, lower_B, rem_patch.resolution, rem_patch.upper_vol_patch, rem_patch.lower_vol_patch, rem_patch.upper_water_patch, rem_patch.lower_water_patch, fixed_midplane_r = rem_patch.fixed_midplane_r, fixed_midplane_z = rem_patch.fixed_midplane_z, midplane_area_patch = rem_patch.midplane_area_patch, upper_coord_Φ_order = rem_patch.upper_coord_Φ_order, lower_coord_Φ_order = rem_patch.lower_coord_Φ_order)
    new_rem_patches = RemeshingPatch[]
    for (i, rem_patch) in enumerate(rem_patches)
        if i ≠ which
            push!(new_rem_patches, rem_patch)
            continue
        end
        push!(new_rem_patches, rem_patch_A)
        push!(new_rem_patches, rem_patch_B)
    end
    return new_rem_patches
end

function change_rem_patch_area_patch(rem_patch::RemeshingPatch, midplane_area_patch::Union{VolumePatch, Nothing})
    return RemeshingPatch(rem_patch.upper_coords, rem_patch.mid_coords, rem_patch.lower_coords, rem_patch.resolution, rem_patch.upper_vol_patch, rem_patch.lower_vol_patch, rem_patch.upper_water_patch, rem_patch.lower_water_patch, fixed_midplane_r = rem_patch.fixed_midplane_r, fixed_midplane_z = rem_patch.fixed_midplane_z, midplane_area_patch = midplane_area_patch === nothing ? nothing : VolumePatchLike(midplane_area_patch), upper_coord_Φ_order = rem_patch.upper_coord_Φ_order, lower_coord_Φ_order = rem_patch.lower_coord_Φ_order)
end
