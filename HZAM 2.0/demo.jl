include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

#takes about  1000 s
outcome, phenotype_counts = HZAM.run_one_HZAM_sim(0.5, 300, 1.3; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=30000, max_generations=1000,
    do_plot=true, plot_int=10,
    total_loci=6,
    female_mating_trait_loci=1:6,
    male_mating_trait_loci=1:6,
    hybrid_survival_loci=1:6,
    per_reject_cost=0,
    sigma_disp=0.03
)

