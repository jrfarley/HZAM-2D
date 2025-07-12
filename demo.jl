include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

t1=time()
outcome =
	HZAM.run_one_HZAM_sim(1, 100, 1.1; # these values are 
		# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
		K_total = 20000, max_generations = 100,
		do_plot = false, plot_int = 1,
		total_loci = 6,
		female_mating_trait_loci = 1:6,
		male_mating_trait_loci = 1:6,
		hybrid_survival_loci = 1:6,
		per_reject_cost = 0.05,
		sigma_disp = 0.03,
		sigma_comp = 0.01,
		track_population_data = true,
		exit_early = false,
	)

println(time()-t1)



@save "sample_outcome_full_pleiotropy.jld2" outcome
