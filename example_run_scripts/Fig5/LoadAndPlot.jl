# get from args 

const MAIN_PATH = "results/plots/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 60.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 175.0
const PLOTTING_MIN_Z = -10.0
const Φ_MAX = 0.1
const Φ_MIN = 0.0

# is the system unstable? then setting it to yes will improve stability
const IS_UNSTABLE = false # false or true
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

for which in ["open/0.76/1.6", "close/0.82/2.0", "open/0.97/1.6"]
    coords, indep_coords, vol_patches, rem_patches = load_state("results/$which/states/final.coords", "results/$which/states/final.vol_patches", "results/$which/states/final.rem_patches")
    # find minimal z of coords
    min_z = minimum([coord.z for coord in coords])
    # shift all coords to have minimal z at 0
    for coord in coords
        coord.z -= min_z
    end
    coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = [coords], [indep_coords], [vol_patches], [rem_patches]
    plot_reaction(rem_patches_arr)
end