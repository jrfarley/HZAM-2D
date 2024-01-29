include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

#takes about  1000 s
output = HZAM.run_one_HZAM_sim(0.95, 500, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=40000, max_generations=50,
    do_plot=true, plot_int=10,
    total_loci=12,
    female_mating_trait_loci=1:3,
    male_mating_trait_loci=4:6,
    hybrid_survival_loci=7:9,
    per_reject_cost=0,
    gene_plot=false,
    save_plot=false,
    track_population_data=false
)

@save "output.JLD2" output