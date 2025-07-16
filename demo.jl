include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

t1=time()
outcome =
	HZAM.run_one_HZAM_sim(1, 300, 1.1; # these values are 
		# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
		K_total = 20000, max_generations = 3000,
		do_plot = true, plot_int = 10,
		total_loci = 12,
		female_mating_trait_loci = 1:4,
		male_mating_trait_loci = 5:8,
		hybrid_survival_loci = 9:12,
		per_reject_cost = 0.1,
		sigma_disp = 0.03,
		sigma_comp = 0.01,
		gene_plot=false
	)

println(time()-t1)



@save "sample_outcome_full_pleiotropy2.jld2" outcome
