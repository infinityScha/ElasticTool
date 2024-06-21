WHICH_SYSTEM = ARGS[1]

ν = WHICH_SYSTEM == "open" ? 1.0 : 0.91
m₀_4π = WHICH_SYSTEM == "open" ? 2.0 : 3.0

const MAIN_PATH = "results/$WHICH_SYSTEM/$ν/$m₀_4π/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 70.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 120.0
const PLOTTING_MIN_Z = -70.0
const Φ_MAX = 0.1
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

include("../../ElasticTool.jl")

# the main part of the run
# STEP 1: MINIMIZE EDGES
rem_patches_open = vesicle_manuscript_ADE(ν=ν, m₀_4π=m₀_4π, is_open = true)[end]
rem_patches_bud = vesicle_manuscript_ADE(ν=ν, m₀_4π=m₀_4π, is_open = false)[end]
if WHICH_SYSTEM != "open"
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = interpolate(rem_patches_open, rem_patches_bud, 7, 2, 0.2, keep_edges=false)
    coords, indep_coords, vol_patches, rem_patches = coords_arr[end], indep_coords_arr[end], vol_patches_arr[end], rem_patches_arr[end]
else
    coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches_open)
end
coords, indep_coords, vol_patches, rem_patches = remesh(rem_patches)

which_to_adjust = [1]
when_to_adjust = [5]
adjust_resolution_to = [0.5]
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, which_to_adjust, when_to_adjust, adjust_resolution_to, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4)

# remove the fig folder if it exists
if isdir(MAIN_PATH*"figs")
    rm(MAIN_PATH*"figs", recursive=true)
end
# create the fig folder
mkdir(MAIN_PATH*"figs")
# save the final configuration plot 
plot_coords(vol_patches)