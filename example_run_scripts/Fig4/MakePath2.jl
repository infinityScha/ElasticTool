
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

using DelimitedFiles

# copy the path found for mixed
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string("results/mixed/path_reset2/final_path/", M)
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
elseif WHICH_SYSTEM == "20"
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.2)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.2)
elseif WHICH_SYSTEM == "40"
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.1)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.4)
elseif WHICH_SYSTEM == "none"
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS_I", 0.0)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "JS", 0.0)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "CONSTR_PHI", 0)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "USE_PHI", 0)
    change_vol_patch_property_in_string(MAIN_PATH * "initial_path/", M, 1, "PHI0", 0.0)
end

if WHICH_SYSTEM != "mixed"
    # iterate over all the coordinate files of the string and change the column of phi to 0 and set phi to being fixed or 0.5 depending on the system
    for i in 1:M
        # remove padding
        lines = readlines(MAIN_PATH * "initial_path/$i.coords")
        lines = [replace(line, r"\s+" => " ") for line in lines]
        # save the unpadded file
        open(MAIN_PATH * "initial_path/$i.coords", "w") do io
            for line in lines
                println(io, line)
            end
        end
        # load the unpadded file as a DataFrame
        df = CSV.read(MAIN_PATH * "initial_path/$i.coords", DataFrame, delim=' ', header=true)
        # change the column of phi to 0.5
        df.Φ .= 0.5
        if WHICH_SYSTEM != "50"
            # set phi to being fixed
            df.fixed_Φ .= 1
        end
        # save the DataFrame in the same format as the one loaded
        CSV.write(MAIN_PATH * "initial_path/$i.coords", df, delim=' ', header=true)
    end
end
# load the initial path
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = load_string(MAIN_PATH * "initial_path/", M)

neb_coords_place_arr = [NEBCoordinatePlace(1, false), NEBCoordinatePlace(2, false)]
# minimize the path
coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr = neb_method(coords_arr, indep_coords_arr, vol_patches_arr, rem_patches_arr, neb_coords_place_arr, initial_step_num = 0, initial_ΔEₘₐₓ=1.0e-2, final_ΔEₘₐₓ=5.0e-4, remesh_spacing=1, μᵥ=5.0e4, μᵩ=5.0e6, mult_val=0.999)

# remove the fig folder if it exists
if isdir(MAIN_PATH*"figs")
    rm(MAIN_PATH*"figs", recursive=true)
end
# create the fig folder
mkdir(MAIN_PATH*"figs")

# save the final configuration plot
for i in 1:M
    plot_coords(vol_patches_arr[i])
end
# save reaction coordinate plot
plot_reaction(rem_patches_arr)

save_string(coords_arr, vol_patches_arr, rem_patches_arr, MAIN_PATH * "final_path/")
println("Final path saved.")