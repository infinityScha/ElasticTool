const MAIN_PATH = "results_domain2/-0.05/2.6/"

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
const CHI = 2.6
# smooth line parameter
const Kₗᵢₙₑ = 0.2 * 0.639
# for energy report
const ENERGY_FACTOR = 1.0

include("../../ElasticTool.jl")

# the main part of the run
INIT_PATH = "results/-0.05/0/"
coords, indep_coords, vol_patches, rem_patches = load_state(INIT_PATH*"initial.coords", INIT_PATH*"initial.vol_patches", INIT_PATH*"initial.rem_patches")
# get all the r and phi values of the upper leaflet of the first rem_patch
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
k *= 6
k = floor(Int, k)
new_Φ₀ = 0.02
new_Φ₁ = 0.98
# change the initial state to the new Φ₀ and Φ₁
for coord in rem_patch.upper_coords[1:k]
    coord.Φ₊ = Φ₊(new_Φ₀)
end
for coord in rem_patch.upper_coords[k+1:end]
    coord.Φ₊ = Φ₊(new_Φ₁)
end
which_to_adjust = [1,1,1]
when_to_adjust = [100,200,300]
adjust_resolution_to = [0.8,0.65,0.5]
coords, indep_coords, vol_patches, rem_patches = minimize(coords, indep_coords, vol_patches, rem_patches, which_to_adjust, when_to_adjust, adjust_resolution_to, μᵥ=5.0e4, μᵩ=5.0e6, remesh_spacing=1, mult_val=0.999, final_ΔEₘₐₓ=2.5e-4, initial_ΔEₘₐₓ=1.0)
