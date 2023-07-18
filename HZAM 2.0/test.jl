include("HZAM/src/HZAM.jl")

import .HZAM

outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("trial_3_ecolDiff_0", 1.1, 0, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_3_ecolDiff_0", outcomes, w_hybs, S_AMs)

#=outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 0.5, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_0.5", outcomes, w_hybs, S_AMs)

outcomes, w_hybs, S_AMs = HZAM.run_HZAM_set_ecolDiff("AAAAAAAAAAAAAAAAA", 1.1, 1.0, K_total=20000, max_generations=100)

HZAM.make_and_save_figs2("HZAM_Sym_Julia_results_GitIgnore", "trial_2_ecolDiff_1", outcomes, w_hybs, S_AMs)=#

