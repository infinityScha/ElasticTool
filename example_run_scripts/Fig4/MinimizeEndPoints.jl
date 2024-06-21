WHICH_SYSTEM = ARGS[1]
WHICH_ENDPOINT = ARGS[2]


const MAIN_PATH = "results/$WHICH_SYSTEM/$WHICH_ENDPOINT/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 130.0
const PLOTTING_MIN_Z = -15.0
const Φ_MAX = 0.65
const Φ_MIN = 0.35

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

include("../../ElasticTool.jl")

# make the initial configuration
# the main part of the run
# STEP 1: MINIMIZE EDGES
if WHICH_SYSTEM == "mixed"
    rem_patches_open = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = true, is_binary = true, Φ₀=0.1, Jₛᵘ=0.4, r_tot=50.0)[end]
    rem_patches_bud = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = false, is_binary = true, Φ₀=0.1, Jₛᵘ=0.4, r_tot=50.0)[end]
elseif WHICH_SYSTEM == "upper"
    rem_patches_open = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = true, is_binary = false, Jₛ_upper=0.04, r_tot=50.0)[end]
    rem_patches_bud = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = false, is_binary = false, Jₛ_upper=0.04, r_tot=50.0)[end]
elseif WHICH_SYSTEM == "both"
    rem_patches_open = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = true, is_binary = false, Jₛ_upper=0.02, Jₛ_lower=-0.02, r_tot=50.0)[end]
    rem_patches_bud = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.98, is_open = false, is_binary = false, Jₛ_upper=0.02, Jₛ_lower=-0.02, r_tot=50.0)[end]
end

if WHICH_ENDPOINT == "final"
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = interpolate(rem_patches_open, rem_patches_bud, 10, 2, 0.15, keep_edges=false)
    coords, indep_coords, vol_patches, rem_patches = coords_arr[end], indep_coords_arr[end], vol_patches_arr[end], rem_patches_arr[end]
else
    coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches_open)
end

function anchor_the_membrane(coords, vol_patches, rem_patches)
    rem_patch = rem_patches[1]
    ds_arr = ds_arr_over_coords(rem_patch)
    ds_arr = -ds_arr .+ ds_arr[end]
    # find point closest to 20
    split_idx = argmin(abs.(ds_arr .- 20.0))
    rem_patches = split_rem_patch(rem_patches, 1, split_idx)
    new_area_patch = AreaPatch(20.0)
    rem_patches[2] = change_rem_patch_area_patch(rem_patches[2], new_area_patch)        

    coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)

    coord_to_anchor = rem_patches[2].upper_coords[1]
    # find by how much to shift all coordinates so that the anchored point is at 0
    shift_z = - coord_to_anchor.z
    # shift all coordinates z by shift_z
    for coord in coords
        coord.z += shift_z
    end

    # save the initial configuration
    save_state(coords, vol_patches, rem_patches, MAIN_PATH * "states/final.coords", MAIN_PATH * "states/final.vol_patches", MAIN_PATH * "states/final.rem_patches")

    # change the anchored point to have a fixed z-coordinate
    change_coord_property(MAIN_PATH * "states/final.coords", coord_to_anchor.idx, "fixed_z", 1)
    
    # reload the initial configuration
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH * "states/final.coords", MAIN_PATH * "states/final.vol_patches", MAIN_PATH * "states/final.rem_patches")

    return coords, indep_coords, vol_patches, rem_patches
end

function unanchor_the_membrane(coords, vol_patches, rem_patches)
    rem_patch = rem_patches[1]

    coord_to_unanchor = rem_patches[1].lower_coords[end]

    # shift this point to 0
    shift_z = - coord_to_unanchor.z

    for coord in coords
        coord.z += shift_z
    end

    # save the initial configuration
    save_state(coords, vol_patches, rem_patches, MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches") 

    # set the previouly anchored point at the bottom to not have a fixed z-coordinate
    change_coord_property(MAIN_PATH * "initial.coords", coord_to_unanchor.idx, "fixed_z", 0)

    # reload the initial configuration
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")

    return coords, indep_coords, vol_patches, rem_patches
end

coords, indep_coords, vol_patches, rem_patches = unanchor_the_membrane(coords, vol_patches, rem_patches)

"""
# for initial configurations slowly reduce the volume
if WHICH_ENDPOINT == "initial"
    initial_vol_w = parse(Float64, get_patch_property(MAIN_PATH * "initial.vol_patches", 3, "V0")) / 0.98
    initial_vol_u = parse(Float64, get_patch_property(MAIN_PATH * "initial.vol_patches", 1, "V0")) / 1.1
    change_patch_property(MAIN_PATH * "initial.vol_patches", 3, "V0", initial_vol_w)
    change_patch_property(MAIN_PATH * "initial.vol_patches", 1, "V0", initial_vol_u)
    multiplies1 = range(1.0, 0.98, 10)
    multiplies2 = range(1.0, 1.1, 10)
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")
    for (mult1, mult2) in zip(multiplies1, multiplies2)
        global coords, indep_coords, vol_patches, rem_patches
        save_state(coords, vol_patches, rem_patches, MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")
        change_patch_property(MAIN_PATH * "initial.vol_patches", 3, "V0", initial_vol_w * mult1)
        change_patch_property(MAIN_PATH * "initial.vol_patches", 1, "V0", initial_vol_u * mult2)
        coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")
        coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, max_steps=300)
    end
end
"""

which_to_adjust = [1,2,1,2,1,2,1,2]
when_to_adjust = [500,500,1000,1000,2000,2000,4000,4000]
adjust_resolution_to = [0.85,0.85,0.75,0.75,0.625,0.625,0.5,0.5]
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, which_to_adjust, when_to_adjust, adjust_resolution_to, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4)

# anchor_the_membrane
coords, indep_coords, vol_patches, rem_patches = anchor_the_membrane(coords, vol_patches, rem_patches)
# minimize with the anchored membrane
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4)

# remove the fig folder if it exists
if isdir(MAIN_PATH*"figs")
    rm(MAIN_PATH*"figs", recursive=true)
end
# create the fig folder
mkdir(MAIN_PATH*"figs")
# save the final configuration plot 
plot_coords(vol_patches)