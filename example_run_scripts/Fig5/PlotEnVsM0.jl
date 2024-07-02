const MAIN_PATH = "results/filler/"

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

using DataFrames

#νs = [0.76, 0.79, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97]
#m₀s = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

νs = [0.76]
m₀s = [1.65, 1.7, 1.8, 1.9, 2.0]

function get_energy(which, ν, m₀)
    println("results/$which/$ν/$m₀/energy.csv")
    # read energy.csv
    res = CSV.read("results/$which/$ν/$m₀/energy.csv", DataFrame)
    # get energy, last line
    return res[end, :energy]
end

for ν in νs
    poly_arr = []
    m₀_dict, en_dict = Dict(), Dict()
    for which in ["open", "close"]
        m₀_dict[which], en_dict[which] = [], []
        for m₀ in m₀s
            # load coords of final config
            coords, indep_coords, vol_patches, rem_patches = load_state("results/$which/$ν/$m₀/states/final.coords", "results/$which/$ν/$m₀/states/final.vol_patches", "results/$which/$ν/$m₀/states/final.rem_patches")
            # check if the state has a neck of less than 7 nm
            rem_patch = rem_patches[1]
            upper_rs = [coord.r for coord in rem_patch.upper_coords]
            found_min = false
            # find the first maximum of r
            k = 1
            while upper_rs[k] < upper_rs[k+1]
                k += 1
            end
            # now find the first minimum of r
            for K in k:length(upper_rs)-10
                if upper_rs[K] < upper_rs[K+1]
                    found_min = true
                    k = K
                    break
                end
            end
            has_neck = found_min && (upper_rs[k] < 20.0) #15.0)
            
            #println(has_neck, " ", which, " ", ν, " ", m₀, " ", k, " ", upper_rs[k])
            
           if has_neck && (which == "open")
                continue
           end

           if !has_neck && (which == "close")
                continue
           end

            m₀_dict[which] = push!(m₀_dict[which], m₀)
            en_dict[which] = push!(en_dict[which], get_energy(which, ν, m₀))
        end
    end
    m₀_for_diff, en_diff = [], []
    for m₀ in m₀s
        if m₀ in m₀_dict["open"] && m₀ in m₀_dict["close"]
            loc_open = findfirst(m₀_dict["open"] .== m₀)
            loc_close = findfirst(m₀_dict["close"] .== m₀)
            en_diff = push!(en_diff, en_dict["open"][loc_open] - en_dict["close"][loc_close])
            m₀_for_diff = push!(m₀_for_diff, m₀)
        end
    end
    plt.scatter(m₀_for_diff, en_diff, label="ΔE")
    # find 3 closest points to zero in en_diff
    m₀_for_diff_copy, en_diff_copy = copy(m₀_for_diff), copy(en_diff)
    m₀_for_for_fit, en_for_fit = [], []
    for i in 1:3
        idx = argmin(abs.(en_diff_copy))
        m₀ = m₀_for_diff_copy[idx]
        en = en_diff_copy[idx]
        m₀_for_for_fit = push!(m₀_for_for_fit, m₀)
        en_for_fit = push!(en_for_fit, en)
        deleteat!(m₀_for_diff_copy, idx)
        deleteat!(en_diff_copy, idx)
    end
    p = np.polyfit(m₀_for_for_fit, en_for_fit, 2)
    m₀_fit = range(1.5, 3.0, length=100)
    en_fit = p[1] .* m₀_fit.^2 .+ p[2] .* m₀_fit .+ p[3]
    plt.plot(m₀_fit, en_fit, label="quadratic fit near the root", linestyle="--")   

    # find roots of the quadratic equation
    m₀_int = (-p[2] + sqrt(p[2]^2 - 4 * p[1] * p[3])) / (2 * p[1])
    # keep only 3 digits
    m₀_int = round(m₀_int, digits=3)
    en_int = 0.0

    plt.scatter([m₀_int], [en_int], label="intersection", color="black")
    plt.xlabel("m₀/4π")
    plt.ylabel("Energy [kT]")
    plt.xlim(1.5, 3.0)
    plt.title("ν = $ν, l₀ = 1.5 nm, intersection at m₀/4π = $m₀_int")
    plt.legend()
    plt.savefig("results/energies/$ν.pdf")
    plt.clf()
end
