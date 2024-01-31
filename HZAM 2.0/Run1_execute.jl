include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots

println("Number of threads is $(Threads.nthreads()) ")

HZAM.run_HZAM_sets_complete("Run1_31Jan2024")
