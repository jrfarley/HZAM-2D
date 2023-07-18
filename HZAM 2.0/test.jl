include("HZAM/src/HZAM.jl")

import .HZAM

#outcomes, sim_params = HZAM.run_HZAM_set("testing_filesystem", 1.1, 0.5, K_total=20000, max_generations=1)

sim_params, outcomes = HZAM.load_from_file("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/testing_filesystem")

HZAM.make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", "testing_filesystem", outcomes, sim_params)

#outcomes, sim_params = HZAM.run_HZAM_set("test_fylesystem", 1.1, 0.5, K_total=20000, max_generations=1)

#HZAM.make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", "trial_3_ecolDiff_0.5", outcomes, w_hybs, S_AMs)

#=outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 0.5, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_0.5", outcomes, w_hybs, S_AMs)

outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 1.0, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_1", outcomes, w_hybs, S_AMs)=#

