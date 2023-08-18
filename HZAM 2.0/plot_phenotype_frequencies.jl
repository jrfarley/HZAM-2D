include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2

@load "HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/phenotype_frequencies/phenotypes2.JLD2" phenotypes

HZAM.plot_phenotypes(phenotypes)