using Distributed
interrupt()
rmprocs(20)
addprocs(4)

@everywhere begin
	using Pkg
	Pkg.activate("HZAM")
	Pkg.instantiate()
	Pkg.precompile()
end

@everywhere begin
	include("HZAM/src/HZAM.jl")
	include("make_figures/make_main_plot.jl")
	import .HZAM

	# creates a folder to store the simulation results in
	HZAM.set_results_folder("HZAM-J_2D_results/NEW_HZAM_RUN/")
end

# setups 6 and 7 (low/high search cost with full pleiotropy) are omitted 
set_numbers = [1, 2, 3, 4, 5, 8, 9, 10, 11]

# run 3 replicates of each simulation setup
for replicate in 1:3
	HZAM.run_HZAM_sets_complete_one_locus("replicate_$replicate"; set_numbers)
end


completed = false
global replicate = 4

# continue running simulations until each set of simulation parameters has an outcome 
# type that has occurred 3 times
while !completed
	completed = HZAM.first_to_three(HZAM.results_folder, "replicate_$replicate"; set_numbers, categorization_loci = "mmt")
	global replicate += 1
end

# save the data in a condensed format that's easier to work with
HZAM.extract_all("HZAM-J_2D_results/NEW_HZAM_RUN", "HZAM-J_2D_results_categorized/NEW_HZAM_RUN"; set_numbers = [1, 2, 3, 4, 5, 8, 9, 10, 11])
