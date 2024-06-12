using JLD2

include("HZAM-J_beta.jl")
#include("make_density_plots.jl")
#=
genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(0, 1, 0, 1.1, do_plot=false, K_total =1000)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb0_S_AM1_ecolDiff1_full_pleiotropy" hybrid_indices locations


genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(0.8, 1, 0, 1.1, do_plot=false, K_total =4000, max_generations = 1000)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb0.8_S_AM1_ecolDiff1_full_pleiotropy" hybrid_indices locations


genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(0.95, 1, 0, 1.1, do_plot=false, K_total =4000, max_generations = 1000)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb0.95_S_AM1_ecolDiff1_full_pleiotropy" hybrid_indices locations


genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(1, 1, 0, 1.1, do_plot=false, K_total =10000, max_generations = 1000)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb1_S_AM1_ecolDiff1_full_pleiotropy" hybrid_indices locations

genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(1, Inf, 0, 1.1, do_plot=false, K_total =10000, max_generations = 1000)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb1_S_AM_inf_ecolDiff1_full_pleiotropy" hybrid_indices locations


=#


genotypes_F, locations_F, genotypes_M, locations_M = run_one_HZAM_sim(0.9, 100, 0, 1.1, do_plot=false, K_total =10000, max_generations = 100, sigma_disp=0.07)

female_traits = calc_traits_additive(genotypes_F[:, 1:3, :])
male_traits = calc_traits_additive(genotypes_M[:, 1:3, :])

hybrid_indices = [female_traits; male_traits]
locations = [locations_F; locations_M]

@save "w_hyb0.8_S_AM100_full_pleiotropy" hybrid_indices locations

