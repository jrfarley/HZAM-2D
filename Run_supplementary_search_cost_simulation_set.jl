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
all_jobs = []
# setups 6 and 7 (low/high search cost with full pleiotropy) are omitted 
set_numbers = [1, 2, 3, 4, 5, 8, 9, 10, 11]

S_AM_set = [3, 10, 30, 100, 300, 1000]
sc_set = [0.0, 0.001, 0.003, 0.01, 0.02, 0.05, 0.1]



function setup_search_cost_supplement_set(
	set_name::String,
	total_loci::Int,
	female_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	male_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	hybrid_survival_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	per_reject_cost_set::Union{UnitRange{<:Real}, Vector{<:Real}},
	S_AM_set::Union{UnitRange{<:Real}, Vector{<:Real}},
	dir::String;
	intrinsic_R::Real = 1.1, K_total::Int = 20000, max_generations::Int = 1500,
	survival_fitness_method::String = "epistasis", sigma_disp = 0.03,
	cline_width_loci = "mmt",
	w_hyb = 1,
	mating_preference_type = "additive",
)
	dir = mkpath(string(dir, "/", set_name))

	index_array = Array{Missing, 2}(
		undef, length(per_reject_cost_set), length(S_AM_set),
	)

	# Loop through the different simulation sets
	for i in CartesianIndices(index_array)
		sc = per_reject_cost_set[i[1]]
		S_AM = S_AM_set[i[2]]

		run_name =
			string(
				"HZAM_simulation_run_gen", max_generations, "_SC",
				sc, "_Whyb", w_hyb, "_SAM", S_AM,
			)

		filename = string(dir, "/", run_name, ".jld2")

		if isfile(filename)
			println("Already completed: $run_name")
			continue
		end

		p = HZAM.JobParams(
			w_hyb,
			S_AM,
			intrinsic_R,
			K_total,
			max_generations,
			total_loci,
			female_mating_trait_loci,
			male_mating_trait_loci,
			hybrid_survival_loci,
			survival_fitness_method,
			sc,
			sigma_disp,
			run_name,
			cline_width_loci,
			mating_preference_type,
			filename,
		)
		push!(all_jobs, p)
	end # of S_AM loop
end

function queue_jobs()
	# One FMT locus
	for i in 1:4
		setup_search_cost_supplement_set(
			"three_mmt_one_fmt_w_hyb_09_$i",
			9,
			[1],
			4:6,
			7:9,
			sc_set,
			S_AM_set,
			"HZAM-J_2D_results/supplement";
			w_hyb = 0.9,
		)

		setup_search_cost_supplement_set(
			"three_mmt_one_fmt_w_hyb_1_$i",
			9,
			[1],
			4:6,
			7:9,
			sc_set,
			S_AM_set,
			"HZAM-J_2D_results/supplement";
			w_hyb = 1.0,
		)
	end

	# Three fmt loci
	for i in 1:4
		setup_search_cost_supplement_set(
			"three_mmt_three_fmt_w_hyb_09_$i",
			9,
			1:3,
			4:6,
			7:9,
			sc_set,
			S_AM_set,
			"HZAM-J_2D_results/supplement";
			w_hyb = 0.9,
		)

		setup_search_cost_supplement_set(
			"three_mmt_three_fmt_w_hyb_1_$i",
			9,
			1:3,
			4:6,
			7:9,
			sc_set,
			S_AM_set,
			"HZAM-J_2D_results/supplement";
			w_hyb = 1.0,
		)
	end

	# Multivariate preference
	for i in 1:4
		setup_search_cost_supplement_set(
			"multivariate_preference_w_hyb_09_$i",
			9,
			1:3,
			4:6,
			7:9,
			sc_set,
			S_AM_set,
			"HZAM-J_2D_results/supplement";
			w_hyb = 0.9,
			mating_preference_type = "multivariate",
		)

		#=setup_search_cost_supplement_set(
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
		)=#
	end
end

@everywhere function worker(params)
	for p in params
		HZAM.job(p)
	end
end

# launch the RemoteChannel
params = RemoteChannel(()->Channel(Inf))

queue_jobs()
println(all_jobs)
for p in all_jobs
	put!(params, p)
end
close(params)

@sync begin
	for w in workers()
		@async remotecall_wait(x -> invokelatest(worker, x), w, params)
	end
end



#=
# save the data in a condensed format that's easier to work with
folders = [
	"three_mmt_one_fmt_w_hyb_09",
	"three_mmt_one_fmt_w_hyb_1",
	"three_mmt_three_fmt_w_hyb_09",
	"three_mmt_three_fmt_w_hyb_1",
	"multivariate_preference_w_hyb_09",
	"multivariate_preference_w_hyb_1",
]

for i in 1:4
	for j in 1:5
		extract_data("HZAM-J_2D_results/supplement", "$(folders[j])_$i", "HZAM-J_2D_results_categorized/supplement", "mmt")
	end
end
=#


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