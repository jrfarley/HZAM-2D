
using Distributed
rmprocs(4)
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



HZAM.run_HZAM_sets_complete(
    "Run1_02Feb2024",
    set_numbers=[1,2,3,8,9]
)
