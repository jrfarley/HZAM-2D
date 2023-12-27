include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format


outcome = HZAM.run_one_HZAM_sim(0.8, 50, 1.0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=25000, max_generations=100,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=false, plot_int=10,
    total_loci=8,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=1:4,
    competition_trait_loci=1:4,
    hybrid_survival_loci=1:4,
    per_reject_cost=0,
    gene_plot=false,
    save_plot=true,
    track_population_data=true
)
println(outcome.population_overlap)
println(outcome.hybrid_zone_width)

@save "file.JLD2" outcome

HZAM.plot_population_tracking_data("file.JLD2")
