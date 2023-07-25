include("HZAM/src/HZAM.jl")

import .HZAM

for i in 0:0.1:1
    set_name = string("t8_ecolDiff", string(i))
    outcomes, sim_params = HZAM.run_HZAM_set(
        set_name,
        1.1,
        i,
        K_total=20000,
        max_generations=1000,
        total_loci=30,
        female_mating_trait_loci=1:6,
        male_mating_trait_loci=7:12,
        competition_trait_loci=13:18,
        hybrid_survival_loci=19:24
    )

    HZAM.make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", set_name, outcomes, sim_params)
end

#HZAM.make_and_save_figs_gene_flow("HZAM_Sym_Julia_results_GitIgnore", "t6_gf", outcomes, sim_params)

#outcomes, sim_params = HZAM.run_HZAM_set("test_fylesystem", 1.1, 0.5, K_total=20000, max_generations=1)

#HZAM.make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", "trial_3_ecolDiff_0.5", outcomes, w_hybs, S_AMs)

#=outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 0.5, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_0.5", outcomes, w_hybs, S_AMs)

outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 1.0, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_1", outcomes, w_hybs, S_AMs)=#

