
WHICH_SYSTEM = ARGS[1]


const MAIN_PATH = "results/$WHICH_SYSTEM/path/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 150.0
const PLOTTING_MIN_Z = -15.0
const Φ_MAX = WHICH_SYSTEM == "50" ? 0.525 : 0.125
const Φ_MIN = WHICH_SYSTEM == "50" ? 0.475 : 0.075

# is the system unstable? then setting it to yes will improve stability
const IS_UNSTABLE = true # false or true
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

M = 10

include("../../ElasticTool.jl")

if WHICH_SYSTEM == "mixed"
    # load the initial and final configurations
    coords0, indep_coords0, vol_patches0, rem_patches0 = load_state("results/mixed/initial-new/states/final.coords", "results/mixed/initial-new/states/final.vol_patches", "results/mixed/initial-new/states/final.rem_patches")
    coords1, indep_coords1, vol_patches1, rem_patches1 = load_state("results/mixed/final-new/states/final.coords", "results/mixed/final-new/states/final.vol_patches", "results/mixed/final-new/states/final.rem_patches")
    # interpolate between the initial and final configurations
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = interpolate(rem_patches1, rem_patches0, M, 2, 1.0/M, keep_edges=true)
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = coords_arr[end:-1:1], indep_coords_arr[end:-1:1], vol_patches_arr[end:-1:1], rem_patches_arr[end:-1:1]
    # save the initial path
    save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "initial_path/")
    # load the initial path
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string(MAIN_PATH * "initial_path/", M)

    function split_domain_at_neck()
        final_domain = rem_patches_arr[end][1]
        ds_arr = ds_arr_over_coords(final_domain)
        r_arr = [coord.r for coord in final_domain.mid_coords]
        i = 2
        # find the 2nd minimum in r_arr
        # first get to the maximum
        while r_arr[i] > r_arr[i-1]
            i += 1
        end
        # then get to the minimum
        while r_arr[i] < r_arr[i-1]
            i += 1
        end
        area_half = 0.5*ds_arr[i]
        new_rem_patches_arr = Vector{RemeshingPatch}[]
        new_area_patch = AreaPatch(area_half)
        new_area_patch_2 = AreaPatch(area_half)
        for rem_patches in rem_patches_arr
            ds_arr = ds_arr_over_coords(rem_patches[1])
            split_idx = argmin(abs.(ds_arr .- area_half))
            new_rem_patches = split_rem_patch(rem_patches, 1, split_idx)
            new_rem_patches[1] = change_rem_patch_area_patch(new_rem_patches[1], new_area_patch)

            ds_arr = ds_arr_over_coords(new_rem_patches[2])
            split_idx = argmin(abs.(ds_arr .- area_half))
            new_rem_patches = split_rem_patch(new_rem_patches, 2, split_idx)
            new_rem_patches[2] = change_rem_patch_area_patch(new_rem_patches[2], new_area_patch_2)

            push!(new_rem_patches_arr, new_rem_patches)
        end

        # plot each rem_patch of the endpoint
        new_rem_patches_end = new_rem_patches_arr[end]
        for new_rem_patch in new_rem_patches_end
            mid_coords = new_rem_patch.mid_coords
            r_arr = [coord.r for coord in mid_coords]
            z_arr = [coord.z for coord in mid_coords]
            plt.plot(r_arr, z_arr, label="$i")
        end
        # save fig
        plt.savefig(MAIN_PATH * "split_domain_at_neck.png")
        plt.close()

        temp_coords_arr = Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}[]
        temp_indep_coords_arr = Vector{Vector{Coordinate{VolumePatch,VolumeElement{VolumePatch}}}}[]
        temp_vol_patches_arr = Vector{VolumePatch}[]
        temp_rem_patches_arr = Vector{RemeshingPatch}[]
        
        for rem_patches in new_rem_patches_arr
            coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)
            push!(temp_coords_arr, coords)
            push!(temp_indep_coords_arr, indep_coords)
            push!(temp_vol_patches_arr, vol_patches)
            push!(temp_rem_patches_arr, rem_patches)
        end

        return temp_coords_arr, temp_indep_coords_arr, temp_vol_patches_arr, temp_rem_patches_arr
    end

    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = split_domain_at_neck()
    save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "initial_path_new")
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string(MAIN_PATH * "initial_path_new", M)
else
    # copy the path found for mixed
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string("results/mixed/path/final_path/", M)
    save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "initial_path/")
    # change the string according to the system
    if WHICH_SYSTEM == "upper"
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS", 0.04)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "CONSTR_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "USE_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.0)
    elseif WHICH_SYSTEM == "both"
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS", 0.02)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "CONSTR_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "USE_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 2, "JS", -0.02)
    elseif WHICH_SYSTEM == "50"
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.08)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.5)
    elseif WHICH_SYSTEM == "none"
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS", 0.0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "CONSTR_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "USE_PHI", 0)
        change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.0)
    end
    # iterate over all the coordinate files of the string and change the column of phi to 0 and set phi to being fixed or 0.5 depending on the system
    for i in 1:M
        # load the file as a DataFrame
        df = CSV.read(MAIN_PATH * "initial_path/$i.coords", DataFrame)
        # change the column of phi to 0.5
        df.Φ .= 0.5
        if WHICH_SYSTEM != "50"
            # set phi to being fixed
            df.fixed_Φ .= 1
        end
        # save the DataFrame in the same format as the one loaded
        CSV.write(MAIN_PATH * "initial_path/$i.coords", df, delim=' ', header=true)
    end
    # load the initial path
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string(MAIN_PATH * "initial_path/", M)
end
neb_coords_place_arr = [NEBCoordinatePlace(1, false), NEBCoordinatePlace(2, false)]
# minimize the path
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = neb_method(coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr, neb_coords_place_arr, initial_step_num = 0, initial_ΔEₘₐₓ=1.0, final_ΔEₘₐₓ=1.0e-3, remesh_spacing=1, μᵥ=5.0e4, μᵩ=5.0e6, mult_val=0.999)

save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "final_path/")
println("Final path saved.")