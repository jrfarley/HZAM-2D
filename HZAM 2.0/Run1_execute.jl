
using Distributed 
addprocs(4)

@everywhere begin
    using Pkg
    Pkg.activate("HZAM")
    Pkg.instantiate()
    Pkg.precompile()
  end

@everywhere begin
    include("HZAM/src/HZAM.jl")

    import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
end


println("Number of threads is $(Threads.nthreads()) ")

@time HZAM.run_HZAM_sets_complete(
    "Run1_31Jan2024",
    w_hyb_set_of_run=[1],
    S_AM_set_of_run=[10],
    max_generations=100,
    K_total = 10000
)
