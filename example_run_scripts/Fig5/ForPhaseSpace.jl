# get from args 
WHICH_SYSYTEM, ν_PREV, ν_CURR, m₀_4π_PREV, m₀_4π_CURR = ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5]

const MAIN_PATH = "results/$WHICH_SYSYTEM/$ν_CURR/$m₀_4π_CURR/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 120.0
const PLOTTING_MIN_Z = -70.0
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

# check if the run was already done (final.coords in states folder), if so, don't rerun
if isfile(MAIN_PATH*"states/final.coords")
    println("The run was already done, running a bit of the final energy to get the final energy value")
    # check if energy.csv has only one line, if so, then the run was not finished
    if length(readlines(MAIN_PATH*"energy.csv")) == 1
        coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
        # run the minimization
        coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=1.0e4, μᵩ=1.0e4, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, initial_ΔEₘₐₓ=1e-3, max_steps=200)
    end
else
    # the main part of the run
    # copy the previous state to the current state folder
    PREV_PATH = "results/$WHICH_SYSYTEM/$ν_PREV/$m₀_4π_PREV/"
    cp(PREV_PATH*"states/final.coords", MAIN_PATH*"initial.coords", force=true)
    cp(PREV_PATH*"states/final.vol_patches", MAIN_PATH*"initial.vol_patches", force=true)
    cp(PREV_PATH*"states/final.rem_patches", MAIN_PATH*"initial.rem_patches", force=true)
    # change previous state resolution to 0.5 nm
    change_patch_property(MAIN_PATH*"initial.rem_patches", 1, "RESOLUTION", 0.5)
    Vs = [parse(Float64, get_patch_property(MAIN_PATH * "initial.vol_patches", i, "V0")) for i in 1:3]
    # find the 2nd largest index
    order1 = sortperm(Vs, rev=true)
    # make an initial configuration of the ADE model to get all the required volumes from
    vol_patches_open = vesicle_manuscript_ADE(ν=parse(Float64, ν_CURR), m₀_4π=parse(Float64, m₀_4π_CURR), is_open = true)[3]
    order2 = sortperm([vol_patches_open[i].v₀ for i in 1:3], rev=true)
    for (oi1, oi2) in zip(order1, order2)
        change_patch_property(MAIN_PATH*"initial.vol_patches", oi1, "V0", vol_patches_open[oi2].v₀)
    end
    # load the previous state
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"initial.coords", MAIN_PATH*"initial.vol_patches", MAIN_PATH*"initial.rem_patches")
    # run the minimization
    coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=5.0e4, μᵩ=1.0e4, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, initial_ΔEₘₐₓ=1.0)
end
# remove the fig folder if it exists
if isdir(MAIN_PATH*"figs")
    rm(MAIN_PATH*"figs", recursive=true)
end
# create the fig folder
mkdir(MAIN_PATH*"figs")
# main analysis of results
# load the final state
coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
# plot it and save it
plot_coords(vol_patches)