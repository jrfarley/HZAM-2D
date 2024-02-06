include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

@load "HZAM_simulation_run_gen1000_SC0.0_Whyb0.9_SAM1000.0.jld2" outcome
outcome1 = outcome
@load "HZAM_simulation_run_gen1000_SC0.0_Whyb1.0_SAM1.0.jld2" outcome
outcome2 = outcome

println(outcome1.sim_params)
println(outcome1.bimodality)


println(outcome2.sim_params)
println(outcome2.bimodality)