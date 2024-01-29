include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

#=
HZAM.run_HZAM_sets_complete("FIRST POST")
=#

@load "output.JLD2" output

println(output.population_overlap)

println(length(output.population_data.genotypes_F))