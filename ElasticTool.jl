using LinearAlgebra
using SparseArrays
using StaticArrays
using Optim, LineSearches
using Statistics 
using Printf
using Base.Threads
using ProgressMeter
using CSV

# load the modules
include("Geometry.jl")
include("Energetics.jl")
include("Diff.jl")
include("Remesh.jl")
include("SaveStates.jl")
include("Plotting.jl")
include("Builder.jl")
include("Minimizer.jl")
include("Minimizer_NEB.jl")
