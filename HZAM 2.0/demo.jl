include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format


outcome = HZAM.run_one_HZAM_sim(0.8, 50, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=5000, max_generations=300,
    sigma_disp=0.03, sigma_comp=0.01, do_plot=true, plot_int=10,
    total_loci=7,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:7,
    competition_trait_loci=1:4,
    hybrid_survival_loci=1:4,
    per_reject_cost=0.01,
    gene_plot=false,
    save_plot=false,
    track_population_data=false
)
println("Overlap: ", outcome.population_overlap)
println("Width : ", outcome.hybrid_zone_width)

@save "file2.JLD2" outcome

readline()

