WHICH_SYSTEM = ARGS[1]
WHICH_ENDPOINT = ARGS[2]


const MAIN_PATH = "results/$WHICH_SYSTEM/$WHICH_ENDPOINT-new/"

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

# copy the last final state
coords, indep_coords, vol_patches, rem_patches = load_state("results/$WHICH_SYSTEM/$WHICH_ENDPOINT/states/final.coords", "results/$WHICH_SYSTEM/$WHICH_ENDPOINT/states/final.vol_patches", "results/$WHICH_SYSTEM/$WHICH_ENDPOINT/states/final.rem_patches")
# save as the new initial state
save_state(coords, vol_patches, rem_patches, MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")
# change the volume of the water patch (idx 3) by 0.96/0.98
V0_old = get_patch_property(MAIN_PATH*"initial.vol_patches", 3, "V0")
V0_new = parse(Float64, V0_old) * 0.96/0.98
change_patch_property(MAIN_PATH*"initial.vol_patches", 3, "V0", V0_new)
# load the new initial state
coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH * "initial.coords", MAIN_PATH * "initial.vol_patches", MAIN_PATH * "initial.rem_patches")
# minimize
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4)

# remove the fig folder if it exists
if isdir(MAIN_PATH*"figs")
    rm(MAIN_PATH*"figs", recursive=true)
end
# create the fig folder
mkdir(MAIN_PATH*"figs")
# save the final configuration plot 
plot_coords(vol_patches)