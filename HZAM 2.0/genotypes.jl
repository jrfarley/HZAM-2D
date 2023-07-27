include("HZAM/src/HZAM.jl")

import .HZAM

K = 20000

loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15, 16], [17, 18, 19, 20], 1:4, 5:8, 9:12, 13:16))
#=
outcome, pd = HZAM.run_one_HZAM_sim(0.8, 100, 1, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=1000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=10,
    total_loci=20,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:8,
    competition_trait_loci=9:12,
    hybrid_survival_loci=13:16)


=#
filepath = "HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/genotypes_ecolDiff0_w_hyb0.8_S_AM100.jld2"
#=
HZAM.save_genotypes(pd, filepath)


HZAM.check_male_mating_trait(filepath)
=#
#=HZAM.check_competition_trait(filepath)
HZAM.check_neutral_trait(filepath)
HZAM.check_male_mating_trait(filepath)
HZAM.check_female_mating_trait(filepath)
HZAM.check_hybrid_survival(filepath)=#

HZAM.plot_trait_correlations_all(filepath, loci, "Trait correlation ecolDiff0_w_hyb0.8_S_AM100")
println(HZAM.find_extinct_alleles(filepath))