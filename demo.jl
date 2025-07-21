include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

t1=time()
outcome =
	HZAM.run_one_HZAM_sim(0.8, 10, 1.1; # these values are 
		# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
		K_total = 20000, max_generations = 4000,
		do_plot = true, plot_int = 10,
		total_loci = 3,
		female_mating_trait_loci = [2],
		male_mating_trait_loci = [1],
		hybrid_survival_loci = [1],
		per_reject_cost = 0,
		sigma_disp = 0.03,
		sigma_comp = 0.01,
		gene_plot=false
	)

#println(time()-t1)

println(outcome.hybrid_zone_width)

@save "magic_cue.jld2" outcome
