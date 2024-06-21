const MAIN_PATH = "results/mixed/res_0.5/path/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 150.0
const PLOTTING_MIN_Z = -15.0
const Φ_MAX = 0.2
const Φ_MIN = 0.0

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

# copy the path found for mixed
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string("results/mixed/path_reset2/final_path/", M)

# load the initial path
save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "initial_path_res0.5")
# change all rem patches to have resolution of 0.5 nm
for i in 1:length(rem_patches_arr[1])
    change_rem_patch_property_in_string(MAIN_PATH * "initial_path_res0.5/", M, i, "RESOLUTION", 0.5)
end
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string(MAIN_PATH * "initial_path_res0.5", M)

neb_coords_place_arr = [NEBCoordinatePlace(1, false), NEBCoordinatePlace(2, false)]
# minimize the path
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = neb_method(coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr, neb_coords_place_arr, initial_step_num = 0, initial_ΔEₘₐₓ=1.0, final_ΔEₘₐₓ=2.5e-4, remesh_spacing=1, μᵥ=5.0e4, μᵩ=5.0e6, mult_val=0.999)

save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "final_path/")
println("Final path saved.")