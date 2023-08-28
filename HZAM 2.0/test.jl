include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

K = 40000
@load "blending.JLD2" outcome

phenotypes = outcome.phenotypes

HZAM.plot_phenotypes(phenotypes)
