const MAIN_PATH = "results/basic_structure/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 100.0
const PLOTTING_MIN_Z = -70.0
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
# for energy report
const ENERGY_FACTOR = 1.0

include("../../ElasticTool.jl")

# the main part of the run
# STEP 1: MINIMIZE EDGES
rem_patches_open = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.95, is_open = true, is_binary = true, Φ₀=0.1/1.1, Jₛᵘ=0.0, r_tot=50.0)[end]
rem_patches_bud = vesicle_manuscript(Vupper_extra=0.1, Vwater=0.95, is_open = false, is_binary = true, Φ₀=0.1/1.1, Jₛᵘ=0.0, r_tot=50.0)[end]
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = interpolate(rem_patches_open, rem_patches_bud, 15, 2, 0.1, keep_edges=false)
coords, indep_coords, vol_patches, rem_patches = coords_arr[end], indep_coords_arr[end], vol_patches_arr[end], rem_patches_arr[end]
coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)
which_to_adjust = [1,1,1,1]
when_to_adjust = [500,1000,2000,3000]
adjust_resolution_to = [0.8,0.6,0.4,0.3]
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, which_to_adjust, when_to_adjust, adjust_resolution_to, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999)
