include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

outcome =
    HZAM.run_one_HZAM_sim(0.7, 300, 1.1; # these values are 
        # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
        K_total=1000, max_generations=1000,
        do_plot=true, plot_int=10,
        total_loci=9,
        female_mating_trait_loci=1:9,
        male_mating_trait_loci=1:9,
        hybrid_survival_loci=1:9,
        per_reject_cost=0,
        sigma_disp=0.5,
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