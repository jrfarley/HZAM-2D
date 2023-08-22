include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

K = 40000

println(K)
#=
fitnesses, pd, loci, tracking_data = HZAM.run_one_HZAM_sim(0.92, 300, 1, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=400,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=10,
    total_loci=8,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=1:4,
    competition_trait_loci=1:4,
    hybrid_survival_loci=1:4,
    per_reject_cost=0.0,
    gene_plot=false,
    save_plot=false)

@save "mates_per_phenotype.JLD2" fitnesses
HZAM.plot_fitnesses(fitnesses)

#@save "tracking6.JLD2" tracking_data

=#

HZAM.run_one_HZAM_sim(0.8, 400, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=2000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=1,
    total_loci=16,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:8,
    competition_trait_loci=1:4,
    hybrid_survival_loci=1:4,
    per_reject_cost=0.0,
    gene_plot=true,
    save_plot=true
)
