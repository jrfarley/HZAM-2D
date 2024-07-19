include("HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using Plots
using Statistics

#outcome_array=HZAM.load_from_folder("HZAM_Sym_Julia_results_GitIgnore/HZAM_simulation_outcomes_Feb2024_GitIgnore/no_pleiotropy")

#HZAM.plot_overlap(outcome_array)
#HZAM.plot_bimodality(outcome_array)

#println(HZAM.get_new_cline_widths(outcome_array))
#HZAM.methods_plot()
HZAM.compare_phenotype_plots()

#=
HZAM.run_HZAM_set_num_loci(
    "num_loci_comparison_sc_0.05",
    0.05
)
    =#
#=    
outcome_array=HZAM.load_from_folder_num_loci("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes2/num_loci_comparison_sc_0.05")
HZAM.plot_num_loci(outcome_array)
readline()
=#