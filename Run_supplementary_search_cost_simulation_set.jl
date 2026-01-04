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

S_AM_set = [3, 10, 30, 100, 300, 1000]
sc_set = [0.0, 0.001, 0.003, 0.01, 0.02, 0.05, 0.1]

# run 4 replicates of each simulation setup

# One locus per trait
for i in 1:4
	HZAM.run_HZAM_set_search_cost_supplement(
		"one_mmt_one_fmt_w_hyb_095_$i",
		6,
		[1],
		[2],
		4:6,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 0.95,
	)

	HZAM.run_HZAM_set_search_cost_supplement(
		"one_mmt_one_fmt_w_hyb_1_$i",
		6,
		[1],
		[2],
		4:6,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 1.0,
	)
end

# Three fmt loci
for i in 1:4
	HZAM.run_HZAM_set_search_cost_supplement(
		"one_mmt_three_fmt_w_hyb_095_$i",
		9,
		1:3,
		[4],
		7:9,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 0.95,
	)

	HZAM.run_HZAM_set_search_cost_supplement(
		"one_mmt_three_fmt_w_hyb_1_$i",
		9,
		1:3,
		[4],
		7:9,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 1.0,
	)
end

# Multivariate preference
for i in 1:4
	HZAM.run_HZAM_set_search_cost_supplement(
		"multivariate_preference_w_hyb_095_$i",
		9,
		1:3,
		4:6,
		7:9,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 0.95,
		mating_preference_type = "multivariate",
	)

	HZAM.run_HZAM_set_search_cost_supplement(
		"multivariate_preference_w_hyb_1_$i",
		9,
		1:3,
		4:6,
		7:9,
		sc_set,
		S_AM_set,
		"HZAM-J_2D_results/supplement";
		w_hyb = 1.0,
		mating_preference_type = "multivariate",
	)
end

# save the data in a condensed format that's easier to work with
folders = [
	"one_mmt_one_fmt_w_hyb_095",
	"one_mmt_one_fmt_w_hyb_1",
	"one_mmt_three_fmt_w_hyb_095",
	"one_mmt_three_fmt_w_hyb_1",
	"multivariate_preference_w_hyb_095",
	"multivariate_preference_w_hyb_1",
]

for i in 1:4
	for j in 1:6
		extract_data("HZAM-J_2D_results/supplement", "$(folders[j])_$i", "HZAM-J_2D_results_categorized/supplement", "mmt")
	end
end



function extract_data(source_folder, source_setup, output_folder, categorization_loci)
	source_dir = string(source_folder, "/", source_setup)
	println(source_dir)
	files=[]
	try
		files = readdir(source_dir)
	catch
		return
	end
	outcome_array = fill(
		7, length(sc_set), length(S_AM_set),
	)
	parameters_array = Array{Union{HZAM.DataAnalysis.SimParams, Missing}, 2}(
		undef, length(sc_set), length(S_AM_set),
	)
	for file in files
		path = string(source_dir, "/", file)
		if occursin(".jld2", path)
			@load path outcome

			outcome_array[
				indexin(outcome.sim_params.per_reject_cost, sc_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = HZAM.categorize(outcome; categorization_loci)
			parameters_array[
				indexin(outcome.sim_params.per_reject_cost, sc_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = outcome.sim_params
		end
	end
	results_folder = mkpath(output_folder)
	println(source_setup)
	println(outcome_array)
	@save string(results_folder, "/", source_setup, ".JLD2") outcome_array parameters_array
end