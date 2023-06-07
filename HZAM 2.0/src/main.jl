include("HZAM-J_beta_jack.jl")

make_and_save_figs("HZAM_Sym_Julia_results_GitIgnore", "test", convert_to_cat_array(run_HZAM_set("test", 0, 1.1, collect(1:5))))