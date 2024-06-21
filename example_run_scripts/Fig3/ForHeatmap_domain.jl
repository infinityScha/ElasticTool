# get from args 
JS_PREV, JS_CURR, CHI_PREV, CHI_CURR = ARGS[1], ARGS[2], ARGS[3], ARGS[4]

const MAIN_PATH = "results_domain/$JS_CURR/$CHI_CURR/"

# set ploting range here for easier control
const PLOTTING_MAX_R = 100.0
const PLOTTING_MIN_R = 0.0 
const PLOTTING_MAX_Z = 120.0
const PLOTTING_MIN_Z = -70.0
const Φ_MAX = 0.2
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
const CHI = parse(Float64, CHI_CURR)
# smooth line parameter
const Kₗᵢₙₑ = 0.2 * 0.639 # in kT / nm, this gives a penalty of 0.2 kT / nm for a naive linear interface with width of 0.8 nm
# for energy report
const ENERGY_FACTOR = 1.0

include("../../ElasticTool.jl")


PREV_PATH = "results_domain/$JS_PREV/$CHI_PREV/"

# check if the run was already done (final.coords in states folder), if so, don't rerun
if false #isfile(MAIN_PATH*"states/final.coords")
    println("The run was already done, running a bit of the final energy to get the final energy value")
    # check if energy.csv has only one line, if so, then the run was not finished
    if length(readlines(MAIN_PATH*"energy.csv")) == 1
        coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
        # run the minimization
        coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=1.0e4, μᵩ=1.0e4, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, initial_ΔEₘₐₓ=1e-3, max_steps=200)
    end
else
    # copy the previous state to the current state folder
    cp(PREV_PATH*"states/final.coords", MAIN_PATH*"initial.coords", force=true)
    cp(PREV_PATH*"states/final.vol_patches", MAIN_PATH*"initial.vol_patches", force=true)
    cp(PREV_PATH*"states/final.rem_patches", MAIN_PATH*"initial.rem_patches", force=true)
    # set Js value
    change_patch_property(MAIN_PATH*"initial.vol_patches", 1, "JS_I", parse(Float64, JS_CURR))
    # load the previous state
    #coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"initial.coords", MAIN_PATH*"initial.vol_patches", MAIN_PATH*"initial.rem_patches")
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
    # find the index of the last coordinates in the first rem patch
    idx0 = rem_patches[1].upper_coords[end].idx
    idx1 = rem_patches[1].lower_coords[end].idx
    idx2 = rem_patches[1].mid_coords[end].idx
    # resave final state 
    save_state(coords, vol_patches, rem_patches, MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
    # free the z coordinate of the last coordinates in the first rem patch
    change_coord_property(MAIN_PATH*"states/final.coords", idx0, "fixed_z", 0)
    change_coord_property(MAIN_PATH*"states/final.coords", idx1, "fixed_z", 0)
    change_coord_property(MAIN_PATH*"states/final.coords", idx2, "fixed_z", 0)
    change_patch_property(MAIN_PATH*"states/final.rem_patches", 1, "RESOLUTION", 0.3)
    coords, indep_coords, vol_patches, rem_patches = load_state(MAIN_PATH*"states/final.coords", MAIN_PATH*"states/final.vol_patches", MAIN_PATH*"states/final.rem_patches")
    # change previous state resolution to 0.5 nm
    coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, initial_ΔEₘₐₓ=1e-2)
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
"""# get all the r and phi values of the upper leaflet of the first rem_patch
rem_patch = rem_patches[1]
upper_rs = [coord.r for coord in rem_patch.upper_coords]
upper_zs = [coord.z for coord in rem_patch.upper_coords]
upper_Φ₊s = [coord.Φ₊ for coord in rem_patch.upper_coords]
upper_Φs = [Φ(Φ₊) for Φ₊ in upper_Φ₊s]
# get curve length representation
ds = [sqrt((r0 - r1)^2 + (z0 - z1)^2) for (r0, r1, z0, z1) in zip(upper_rs[1:end-1], upper_rs[2:end], upper_zs[1:end-1], upper_zs[2:end])]
dA = [0.5*(r0+r1)*ds_ for (r0, r1, ds_) in zip(upper_rs[1:end-1], upper_rs[2:end], ds)]
dn = [0.5*(Φ0+Φ1)*dA_ for (Φ0, Φ1, dA_) in zip(upper_Φs[1:end-1], upper_Φs[2:end], dA)]
pushfirst!(ds, 0.0)
pushfirst!(dA, 0.0)
pushfirst!(dn, 0.0)
s = cumsum(ds)
A_ = cumsum(dA)
n = cumsum(dn)
# find the first minimum of r, start by finding the first maximum
k = 1
while upper_rs[k] < upper_rs[k+1]
    global k
    k += 1
end
# then find the first minimum
while upper_rs[k] > upper_rs[k+1]
    global k
    k += 1
end
# get the total area up to the first minimum
A₀ = A_[k]
# get the total n up to the first minimum
n₀ = n[k]
Φ₀ = n₀/A₀
# get the total area from the first minimum
A₁ = A_[end] - A₀
# get the total n from the first minimum
n₁ = n[end] - n₀
Φ₁ = n₁/A₁
# delete any lines that have the same chi and js values, use pure julia
lines = readlines("enrichment_domain.csv")
push!(lines, "$CHI_CURR, $JS_CURR, $Φ₀, $Φ₁")
chis = [split(line, ",")[1] for line in lines[2:end]]
jss = [split(line, ",")[2] for line in lines[2:end]]
# omit duplicate lines, keep the latter line
inds = []
for i in 1:length(chis)
    if i in inds
        continue
    end
    for j in i+1:length(chis)
        if chis[i] == chis[j] && jss[i] == jss[j]
            push!(inds, i)
        end
    end
end
for i in inds
    deleteat!(lines, i)
end
open("enrichment_domain.csv", "w") do io
    # save the results to enrichment.csv
    for line in lines
        println(io, line)
    end
end"""