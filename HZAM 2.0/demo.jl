include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

outcome =
    HZAM.run_one_HZAM_sim(1, 2, 1.1; # these values are 
        # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
        K_total=30000, max_generations=2000,
        do_plot=true, plot_int=10,
        total_loci=1,
        female_mating_trait_loci=[1],
        male_mating_trait_loci=[1],
        hybrid_survival_loci=[1],
        per_reject_cost=0.01,
        sigma_disp=0.03,
        sigma_comp=0.01,
        track_population_data=false,
        exit_early = false
    )
#=
HZAM.run_HZAM_set(
	"testing_july19",
	3,
	[1],
	[2],
	[3],
	0,
	[0, 0.5, 1],
	[1, 10, 100],
    K_total=10000
)=#