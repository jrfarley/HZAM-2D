include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

#=
HZAM.run_HZAM_sets_complete("FIRST POST")
=#

outcome_array = HZAM.load_overlap_data_from_folder("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/FIRST POST")
println(outcome_array["full_pleiotropy"])