function _save_coords(coords::Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}, filename::AbstractString)
    file = open(filename, "w")
    println(file, "idx     r           z           Φ           fixed_r fixed_z fixed_Φ")
    for coord in coords
        @printf(file, "%-7d %-11.5f %-11.5f %-8.9f %1d       %1d        %1d\n",
        coord.idx, coord.r, coord.z, Φ(coord.Φ₊),
        coord.fixed_r, coord.fixed_z, coord.fixed_Φ₊)
    end
    close(file)
end

function _save_vol_patches(patches::Vector{VolumePatch}, filename::AbstractString)
    file = open(filename, "w")
    for patch in patches
        println(file, "VOLPATCH $(patch.idx)")
        println(file, "CONSTR_V $(patch.constrain_v)")
        println(file, "V0 $(patch.v₀)")
        println(file, "CONSTR_S $(patch.constrain_s)")
        println(file, "USE_PHI $(patch.use_Φ)")
        println(file, "CONSTR_PHI $(patch.constrain_Φ)")
        println(file, "PHI0 $(patch.Φ₀)")
        println(file, "MU_V $(patch.μᵥ)")
        println(file, "LAMBDA_V $(patch.λᵥ)")
        println(file, "MU_PHI $(patch.μᵩ)")
        println(file, "LAMBDA_PHI $(patch.λᵩ)")
        println(file, "KAPPA $(patch.κ)")
        println(file, "JS $(patch.Jₛ)")
        println(file, "KG $(patch.κᵍ)")
        println(file, "KT $(patch.kₜ)")
        println(file, "KA $(patch.Kₐ)")
        println(file, "H0 $(patch.h₀)")
        println(file, "KAPPA_I $(patch.κᴵ)")
        println(file, "JS_I $(patch.Jₛᴵ)")
        println(file, "KG_I $(patch.κᵍᴵ)")
        println(file, "KT_I $(patch.kₜᴵ)")
        println(file, "KA_I $(patch.Kₐᴵ)")
        println(file, "H0_I $(patch.h₀ᴵ)")
        println(file, "SELF_INTERACTIONS $(patch.self_interactions)")
        println(file, "END VOLPATCH $(patch.idx)")
    end
    close(file)
end

function _save_rem_patches(patches::Vector{RemeshingPatch}, filename::AbstractString)
    file = open(filename, "w")
    for (patch_number, patch) in enumerate(patches)
        println(file, "REMPATCH ", patch_number)
        if patch.upper_coords === nothing
            println(file, "UPPER_COORDS NOTHING")
        else
            coord_indices = [coord.idx for coord in patch.upper_coords]
            coord_indices_str = join(coord_indices, " ")
            println(file, "UPPER_COORDS ", coord_indices_str)
        end
        mid_indices = [coord.idx for coord in patch.mid_coords]
        mid_indices_str = join(mid_indices, " ")
        println(file, "MID_COORDS ", mid_indices_str)
        if patch.lower_coords === nothing
            println(file, "LOWER_COORDS NOTHING")
        else
            coord_indices = [coord.idx for coord in patch.lower_coords]
            coord_indices_str = join(coord_indices, " ")
            println(file, "LOWER_COORDS ", coord_indices_str)
        end
        println(file, "RESOLUTION ", patch.resolution)
        if patch.upper_vol_patch === nothing
            println(file, "UPPER_VOL_PATCH NOTHING")
        else
            println(file, "UPPER_VOL_PATCH ", patch.upper_vol_patch.idx)
        end
        if patch.lower_vol_patch === nothing
            println(file, "LOWER_VOL_PATCH NOTHING")
        else
            println(file, "LOWER_VOL_PATCH ", patch.lower_vol_patch.idx)
        end
        if patch.upper_water_patch === nothing
            println(file, "UPPER_WATER_PATCH NOTHING")
        else
            println(file, "UPPER_WATER_PATCH ", patch.upper_water_patch.idx)
        end
        if patch.lower_water_patch === nothing
            println(file, "LOWER_WATER_PATCH NOTHING")
        else
            println(file, "LOWER_WATER_PATCH ", patch.lower_water_patch.idx)
        end
        println(file, "FIXED_MIDPLANE_R ", patch.fixed_midplane_r)
        println(file, "FIXED_MIDPLANE_Z ", patch.fixed_midplane_z)
        if patch.midplane_area_patch === nothing
            println(file, "MIDPLANE_AREA_PATCH NOTHING")
        else
            println(file, "MIDPLANE_AREA_PATCH ", patch.midplane_area_patch.idx)
        end
        println(file, "END REMPATCH ", patch_number)
    end
    close(file)
end

function _load_coords(filename::AbstractString)
    coords = Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}()

    file = open(filename, "r")
    # Skip the header line
    readline(file)

    for line in eachline(file)
        parts = split(line)
        idx = parse(Int, parts[1])
        r = parse(Float64, parts[2])
        z = parse(Float64, parts[3])
        Φ = parse(Float64, parts[4])
        fixed_r = parse(Int, parts[5])
        fixed_z = parse(Int, parts[6])
        fixed_Φ₊ = parse(Int, parts[7])

        coord = Coordinate{VolumePatch, VolumeElement{VolumePatch}}(idx, r, z, Φ₊(Φ), VolumePatch[], VolumeElement{VolumePatch}[], fixed_r, fixed_z, fixed_Φ₊)
        push!(coords, coord)
    end

    close(file)
    return [coord.idx for coord in coords], coords
end

function _load_vol_patches(filename::AbstractString)
    patches = Vector{VolumePatch}()

    local idx
    local constrain_v
    local v₀
    local constrain_s
    local use_Φ
    local constrain_Φ
    local Φ₀
    local μᵥ
    local λᵥ
    local μᵩ
    local λᵩ
    local κ
    local Jₛ
    local κᵍ
    local kₜ
    local Kₐ
    local h₀
    local κᴵ
    local Jₛᴵ
    local κᵍᴵ
    local kₜᴵ
    local Kₐᴵ
    local h₀ᴵ
    local self_interactions

    open(filename, "r") do file
        for line in eachline(file)
            parts = split(line)
            if length(parts) == 0
                continue
            end
            keyword = parts[1]
            if keyword == "VOLPATCH"
                idx = parse(Int, parts[2])
            elseif keyword == "CONSTR_V"
                constrain_v = parse(Bool, parts[2])
            elseif keyword == "V0"
                v₀ = parse(Float64, parts[2])
            elseif keyword == "CONSTR_S"
                constrain_s = parse(Bool, parts[2])
            elseif keyword == "USE_PHI"
                use_Φ = parse(Bool, parts[2])
            elseif keyword == "CONSTR_PHI"
                constrain_Φ = parse(Bool, parts[2])
            elseif keyword == "PHI0"
                Φ₀ = parse(Float64, parts[2])
            elseif keyword == "MU_V"
                μᵥ = parse(Float64, parts[2])
            elseif keyword == "LAMBDA_V"
                λᵥ = parse(Float64, parts[2])
            elseif keyword == "MU_PHI"
                μᵩ = parse(Float64, parts[2])
            elseif keyword == "LAMBDA_PHI"
                λᵩ = parse(Float64, parts[2])
            elseif keyword == "KAPPA"
                κ = parse(Float64, parts[2])
            elseif keyword == "JS"
                Jₛ = parse(Float64, parts[2])
            elseif keyword == "KG"
                κᵍ = parse(Float64, parts[2])
            elseif keyword == "KT"
                kₜ = parse(Float64, parts[2])
            elseif keyword == "KA"
                Kₐ = parse(Float64, parts[2])
            elseif keyword == "H0"
                h₀ = parse(Float64, parts[2])
            elseif keyword == "KAPPA_I"
                κᴵ = parse(Float64, parts[2])
            elseif keyword == "JS_I"
                Jₛᴵ = parse(Float64, parts[2])
            elseif keyword == "KG_I"
                κᵍᴵ = parse(Float64, parts[2])
            elseif keyword == "KT_I"
                kₜᴵ = parse(Float64, parts[2])
            elseif keyword == "KA_I"
                Kₐᴵ = parse(Float64, parts[2])
            elseif keyword == "H0_I"
                h₀ᴵ = parse(Float64, parts[2])
            elseif keyword == "SELF_INTERACTIONS"
                self_interactions = parse(Bool, parts[2])
            elseif keyword == "END"
                patch = VolumePatch(idx, VolumeElement{VolumePatch}[], constrain_v, v₀, constrain_s, use_Φ, constrain_Φ, Φ₀, 0.0, 0.0, 0.0, μᵥ, λᵥ, μᵩ, λᵩ, κ, κᵍ, Jₛ, kₜ, Kₐ, h₀, κᴵ, κᵍᴵ, Jₛᴵ, kₜᴵ, Kₐᴵ, h₀ᴵ, self_interactions, 0.0)
                push!(patches, patch)
            end
        end
    end
    return [patch.idx for patch in patches], patches
end

function save_state(coords::Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}, vol_patches::Vector{VolumePatch}, rem_patches::Vector{RemeshingPatch},coords_filename::AbstractString, vol_patches_filename::AbstractString, rem_patch_filename::AbstractString)
    _save_coords(coords, coords_filename)
    _save_vol_patches(vol_patches, vol_patches_filename)
    _save_rem_patches(rem_patches, rem_patch_filename)
end

function save_string(coords_arr::Vector{Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}}, vol_patches_arr::Vector{Vector{VolumePatch}}, rem_patches_arr::Vector{Vector{RemeshingPatch}}, dir::AbstractString)
    if !ispath(dir)
        mkpath(dir)
    end
    dir = dir*"/"
    M = length(coords_arr)
    for i in 1:M
        save_state(coords_arr[i], vol_patches_arr[i], rem_patches_arr[i], dir*string(i)*".coords", dir*string(i)*".vol_patches", dir*string(i)*".rem_patches")
    end
end

function load_state(coords_filename::AbstractString, vol_patches_filename::AbstractString, rem_patch_filename::AbstractString)
    patches = Vector{RemeshingPatch}()
    
    coord_indices, coords = _load_coords(coords_filename)
    vol_patch_indices, vol_patches = _load_vol_patches(vol_patches_filename)

   local upper_coords
   local mid_coords
   local lower_coords 
   local resolution 
   local upper_vol_patch 
   local lower_vol_patch 
   local upper_water_patch 
   local lower_water_patch 
   local fixed_midplane_r 
   local fixed_midplane_z
   local midplane_area_patch

   midplane_area_patch = nothing

    open(rem_patch_filename, "r") do file
        for line in eachline(file)
            parts = split(line)
            
            keyword = parts[1]
            
            if keyword == "REMPATCH"
                continue
            elseif keyword == "UPPER_COORDS"
                if parts[2] == "NOTHING"
                    upper_coords = nothing
                else
                    indices = [findfirst(==(idx), coord_indices) for idx in parse.(Int, parts[2:end])] 
                    upper_coords = [coords[i] for i in indices]
                end
            elseif keyword == "MID_COORDS"
                indices = [findfirst(==(idx), coord_indices) for idx in parse.(Int, parts[2:end])] 
                mid_coords = [coords[i] for i in indices]
            elseif keyword == "LOWER_COORDS"
                if parts[2] == "NOTHING"
                    lower_coords = nothing
                else
                    indices = [findfirst(==(idx), coord_indices) for idx in parse.(Int, parts[2:end])] 
                    lower_coords = [coords[i] for i in indices]
                end
            elseif keyword == "RESOLUTION"
                resolution = parse(Float64, parts[2])
            elseif keyword == "UPPER_VOL_PATCH"
                if parts[2] == "NOTHING"
                    upper_vol_patch = nothing
                else
                    upper_vol_patch = vol_patches[findfirst(==(parse(Int, parts[2])), vol_patch_indices)]
                end
            elseif keyword == "LOWER_VOL_PATCH"
                if parts[2] == "NOTHING"
                    lower_vol_patch = nothing
                else
                    lower_vol_patch = vol_patches[findfirst(==(parse(Int, parts[2])), vol_patch_indices)]
                end
            elseif keyword == "UPPER_WATER_PATCH"
                if parts[2] == "NOTHING"
                    upper_water_patch = nothing
                else
                    upper_water_patch = vol_patches[findfirst(==(parse(Int, parts[2])), vol_patch_indices)]
                end
            elseif keyword == "LOWER_WATER_PATCH"
                if parts[2] == "NOTHING"
                    lower_water_patch = nothing
                else
                    lower_water_patch = vol_patches[findfirst(==(parse(Int, parts[2])), vol_patch_indices)]
                end
            elseif keyword == "FIXED_MIDPLANE_R"
                fixed_midplane_r = parse(Bool, parts[2])
            elseif keyword == "FIXED_MIDPLANE_Z"
                fixed_midplane_z = parse(Bool, parts[2])
            elseif keyword == "MIDPLANE_AREA_PATCH"
                if parts[2] == "NOTHING"
                    midplane_area_patch = nothing
                else
                    midplane_area_patch = vol_patches[findfirst(==(parse(Int, parts[2])), vol_patch_indices)]
                end
            elseif keyword == "END"
                patch = RemeshingPatch(upper_coords, mid_coords, lower_coords, resolution, upper_vol_patch, lower_vol_patch, upper_water_patch, lower_water_patch, fixed_midplane_r=fixed_midplane_r, fixed_midplane_z=fixed_midplane_z, midplane_area_patch=midplane_area_patch)
                push!(patches, patch)
            end
        end
    end
    # in the end remesh the patches to clean
    return remesh(patches)
end

function load_string(dir::AbstractString, M::Int64)
    dir = dir*"/"
    coords_arr = Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}[]
    indep_coords_arr = Vector{Vector{Coordinate{VolumePatch, VolumeElement{VolumePatch}}}}[]
    vol_patches_arr = Vector{VolumePatch}[]
    rem_patches_arr = Vector{RemeshingPatch}[]
    for i in 1:M
        coords, indep_coords, vol_patches, rem_patches = load_state(dir*string(i)*".coords", dir*string(i)*".vol_patches", dir*string(i)*".rem_patches")
        push!(coords_arr, coords)
        push!(indep_coords_arr, indep_coords)
        push!(vol_patches_arr, vol_patches)
        push!(rem_patches_arr, rem_patches)
    end
    return coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr
end

function get_patch_property(filename::AbstractString, patch_idx::Int64, property::AbstractString)
    file_contents = read(filename, String)

    patch_start = findfirst("VOLPATCH $patch_idx", file_contents)[1]
    patch_end = findfirst("END VOLPATCH $patch_idx", file_contents)[1]

    if patch_start === nothing || patch_end === nothing
        patch_start = findfirst("REMPATCH $patch_idx", file_contents)[1]
        patch_end = findfirst("END REMPATCH $patch_idx", file_contents)[1]    
    end

    if patch_start === nothing || patch_end === nothing
        error("Patch with index $patch_idx not found in the file.")
    end
    
    patch_string = file_contents[patch_start:patch_end]
    property_start = findfirst(property, patch_string)[1]
    
    if property_start === nothing
        error("Property '$property' not found in the patch $patch_idx.")
    end
    
    property_end = findnext("\n", patch_string, property_start)[1]
    property_line = patch_string[property_start:property_end]

    property_value = split(property_line)[2]
    
    return property_value
end

function change_coord_property(filename::AbstractString, coord_idx::Int64, property::AbstractString, new_value)        
    file_contents = read(filename, String)

    idxs, rs, zs, Φs, fixed_rs, fixed_zs, fixed_Φs = Int[], Float64[], Float64[], Float64[], Int[], Int[], Int[]
    for line in eachline(IOBuffer(file_contents))
        parts = split(line)
        if length(parts) == 0 || parts[1] == "idx"
            continue
        end
        push!(idxs, parse(Int, parts[1]))
        push!(rs, parse(Float64, parts[2]))
        push!(zs, parse(Float64, parts[3]))
        push!(Φs, parse(Float64, parts[4]))
        push!(fixed_rs, parse(Int, parts[5]))
        push!(fixed_zs, parse(Int, parts[6]))
        push!(fixed_Φs, parse(Int, parts[7]))
    end

    coord_idx = findfirst(==(coord_idx), idxs)
    if coord_idx === nothing
        error("Coordinate with index $coord_idx not found in the file.")
    end

    coord_idx = coord_idx[1]

    if property == "r"
        rs[coord_idx] = new_value
    elseif property == "z"
        zs[coord_idx] = new_value
    elseif property == "Φ"
        Φs[coord_idx] = new_value
    elseif property == "fixed_r"
        fixed_rs[coord_idx] = new_value
    elseif property == "fixed_z"
        fixed_zs[coord_idx] = new_value
    elseif property == "fixed_Φ"
        fixed_Φs[coord_idx] = new_value
    else
        error("Property '$property' not found in the coordinate $coord_idx.")
    end

    # save the new coordinates
    open(filename, "w") do file
        println(file, "idx     r           z           Φ           fixed_r fixed_z fixed_Φ")
        for (idx, r, z, Φ, fixed_r, fixed_z, fixed_Φ) in zip(idxs, rs, zs, Φs, fixed_rs, fixed_zs, fixed_Φs)
            @printf(file, "%-7d %-11.5f %-11.5f %-8.9f %1d       %1d        %1d\n"
            , idx, r, z, Φ, fixed_r, fixed_z, fixed_Φ)
        end
    end
end

function change_patch_property(filename::AbstractString, patch_idx::Int64, property::AbstractString, new_value::AbstractString)
    property = property*" "
    file_contents = read(filename, String)

    patch_start = findfirst("VOLPATCH $patch_idx", file_contents)
    patch_end = findfirst("END VOLPATCH $patch_idx", file_contents)

    if patch_start === nothing || patch_end === nothing
        patch_start = findfirst("REMPATCH $patch_idx", file_contents)
        patch_end = findfirst("END REMPATCH $patch_idx", file_contents)    
    end

    if patch_start === nothing || patch_end === nothing
        error("Patch with index $patch_idx not found in the file.")
    end

    patch_start = patch_start[1]
    patch_end = patch_end[1]
    
    patch_string = file_contents[patch_start:patch_end]
    property_start = findfirst(property, patch_string)
    
    if property_start === nothing
        error("Property '$property' not found in the patch $patch_idx.")
    end

    property_start = property_start[1]
    
    property_end = findnext("\n", patch_string, property_start)[1]
    property_line = patch_string[property_start:property_end]
    
    # Modify the property value here
    
    new_patch_string = replace(patch_string, property_line => property*" "*new_value*"\n")
    new_file_contents = replace(file_contents, patch_string => new_patch_string)
    
    open(filename, "w") do file
        write(file, new_file_contents)
    end
end

change_patch_property(filename::AbstractString, patch_idx::Int64, property::AbstractString, new_value) = change_patch_property(filename, patch_idx, property, string(new_value))

function _change_patch_property_in_string(dir::AbstractString, M::Int64, is_volpatch::Bool, patch_idx::Int64, property::AbstractString, new_value)
    extension = is_volpatch ? ".vol_patches" : ".rem_patches"
    for i in 1:M
        change_patch_property(dir*"/"*string(i)*extension, patch_idx, property, new_value)
    end
end

change_vol_patch_property_in_string(dir::AbstractString, M::Int64, patch_idx::Int64, property::AbstractString, new_value) = _change_patch_property_in_string(dir, M, true, patch_idx, property, new_value)
change_rem_patch_property_in_string(dir::AbstractString, M::Int64, patch_idx::Int64, property::AbstractString, new_value) = _change_patch_property_in_string(dir, M, false, patch_idx, property, new_value)
