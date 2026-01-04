include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

outcome = HZAM.run_one_HZAM_sim(1, 10, 1.1; # these values are 
	# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
	K_total = 20000, max_generations = 1500,
	do_plot = true, plot_int = 10,
	total_loci = 6,
	female_mating_trait_loci = 1:3,
	male_mating_trait_loci = 1:3,
	hybrid_survival_loci = 1:3,
	per_reject_cost = 0.01,
	sigma_disp = 0.03,
	sigma_comp = 0.01,
	gene_plot = false,
	cline_width_loci = "mmt",
)


@save "HZAM-J_2D_results/example_outcomes/demo.jld2" outcome
