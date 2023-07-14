include("HZAM/src/HZAM.jl")

import .HZAM

outcomes = HZAM.run_HZAM_set("testing", 0, 1.1, [1]; max_generations=1)

HZAM.make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", "RunName", outcomes)

