include("HZAM/src/HZAM.jl")

import .HZAM

K = 10000

println(K) 


outcome, pd, loci = HZAM.run_one_HZAM_sim(0.7, 100, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=2000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=10,
    total_loci=20,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:8,
    competition_trait_loci=9:12,
    hybrid_survival_loci=5:8,
    per_reject_cost=0.0)

HZAM.plot_fitnesses(outcome.fitness_per_phenotype)

#ClineWidths(0.256065f0, 0.25054f0, 0.279335f0, 0.27077502f0, 0.269615f0, 0.19386499f0, 0.28293502f0)