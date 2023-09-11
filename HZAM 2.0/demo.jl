include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format


outcome = HZAM.run_one_HZAM_sim(0.4, 300, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=1000, max_generations=1000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=10,
    total_loci=16,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:8,
    competition_trait_loci=9:12,
    hybrid_survival_loci=9:12,
    per_reject_cost=0,
    gene_plot=false,
    save_plot=true,
    track_phenotypes=true
)
