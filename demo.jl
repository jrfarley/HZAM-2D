include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format


t1=time()

outcome =
	HZAM.run_one_HZAM_sim(0.9, 1000, 1.1; # these values are 
		# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
		K_total = 20000, max_generations = 1500,
		do_plot = true, plot_int = 10,
		total_loci = 9,
		female_mating_trait_loci = 1:3,
		male_mating_trait_loci = 4:6,
		hybrid_survival_loci = 1:3,
		per_reject_cost = 0.0,
		sigma_disp = 0.03,
		sigma_comp = 0.01,
		gene_plot=false,
		cline_width_loci="mmt"
	)

println(time()-t1)

@save "HZAM_simulation_run_gen1500_SC0.0_Whyb0.9_SAM1000.0.jld2" outcome
