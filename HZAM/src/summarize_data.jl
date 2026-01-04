using ColorSchemes
using Colors
using Dates
using JLD2
using LsqFit: curve_fit
using Statistics
using MannKendall

"Compute the mean value of a vector."
mean(itr) = sum(itr) / length(itr)

"The default directory where all files are saved."
global results_folder = mkpath("HZAM-J_2D_results")


"The set of hybrid fitnesses (w_hyb) values that will be run"
global w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]


"The set of assortative mating strengths (S_AM) values that will be run"
global S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]

"Evenly spaced x coordinates over the range."
global spaced_locations = collect(Float32, 0:0.01:1)

global set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
	"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
	"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
	"low_reject_separate_fmt", "high_reject_separate_fmt"]

"""
	run_HZAM_sets_complete_three_loci(;
		set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
		w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
		S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
		max_generations::Integer = 1500,
		K_total::Integer = 20000,
	)

Run the main set of simulations with three loci per trait and save the outcomes to `results_folder`.

# Arguments

- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: which sets to run (the sets 
are full pleiotropy, no pleiotropy, magic cue, magic preference, etc.).
- `w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which w_hyb values to use.
- `S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which S_AM values to use.
- `max_generations::Integer`: how many generations to run the simulations for.
- `K_total::Integer`: the carrying capacity of the simulated range.
"""
function run_HZAM_sets_complete_three_loci(run_name::String;
	set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
	w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
	S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
	max_generations::Integer = 1500,
	K_total::Integer = 20000,
	cline_width_loci = "mmt",
)
	dir = "HZAM-J_2D_results/$(run_name)_three_loci"

	set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
		"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
		"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
		"low_reject_separate_fmt", "high_reject_separate_fmt"]
	total_loci = [6, 12, 9, 9, 9, 6, 6, 9, 9, 9, 9]
	female_mating_trait_loci = [1:3, 4:6, 1:3, 4:6, 4:6, 1:3, 1:3, 4:6, 4:6, 4:6, 4:6]
	male_mating_trait_loci = [1:3, 7:9, 4:6, 1:3, 4:6, 1:3, 1:3, 7:9, 7:9, 1:3, 1:3]
	hybrid_survival_loci = [1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3]
	per_reject_cost = [0, 0, 0, 0, 0, 0.01, 0.05, 0.01, 0.05, 0.01, 0.05]


	println("$(nprocs()) threads running")
	@async @distributed for i in set_numbers
		run_HZAM_set(
			set_names[i],
			total_loci[i],
			female_mating_trait_loci[i],
			male_mating_trait_loci[i],
			hybrid_survival_loci[i],
			per_reject_cost[i],
			w_hyb_set_of_run,
			S_AM_set_of_run,
			dir;
			max_generations,
			K_total,
			cline_width_loci,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end


"""
	run_HZAM_sets_complete_one_locus(;
		set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
		w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
		S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
		max_generations::Integer = 1500,
		K_total::Integer = 20000,
	)

Run the main set of simulations with one locus per trait and save the outcomes to `results_folder`.

# Arguments

- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: which sets to run (the sets 
are full pleiotropy, no pleiotropy, magic cue, magic preference, etc.).
- `w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which w_hyb values to use.
- `S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which S_AM values to use.
- `max_generations::Integer`: how many generations to run the simulations for.
- `K_total::Integer`: the carrying capacity of the simulated range.
"""
function run_HZAM_sets_complete_one_locus(run_name;
	set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
	w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
	S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
	max_generations::Integer = 1500,
	K_total::Integer = 20000,
)
	dir = "HZAM-J_2D_results/$(run_name)_one_locus"


	set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
		"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
		"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
		"low_reject_separate_fmt", "high_reject_separate_fmt"]
	total_loci = [2, 4, 3, 3, 3, 2, 2, 3, 3, 3, 3]
	female_mating_trait_loci = [[1], [2], [1], [2], [2], [1], [1], [2], [2], [2], [2]]
	male_mating_trait_loci = [[1], [3], [2], [1], [2], [1], [1], [3], [3], [1], [1]]
	hybrid_survival_loci = [[1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1]]
	per_reject_cost = [0, 0, 0, 0, 0, 0.01, 0.05, 0.01, 0.05, 0.01, 0.05]

	println("$(nprocs()) threads running")
	@async @distributed for i in set_numbers
		run_HZAM_set(
			set_names[i],
			total_loci[i],
			female_mating_trait_loci[i],
			male_mating_trait_loci[i],
			hybrid_survival_loci[i],
			per_reject_cost[i],
			w_hyb_set_of_run,
			S_AM_set_of_run,
			dir;
			max_generations,
			K_total,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end

"""
	run_HZAM_sets_complete_nine_loci(;
		set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
		w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
		S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
		max_generations::Integer = 1500,
		K_total::Integer = 20000,
	)

Run the main set of simulations with nine loci per trait and save the outcomes to `results_folder`.

# Arguments

- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: which sets to run (the sets 
are full pleiotropy, no pleiotropy, magic cue, magic preference, etc.).
- `w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which w_hyb values to use.
- `S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which S_AM values to use.
- `max_generations::Integer`: how many generations to run the simulations for.
- `K_total::Integer`: the carrying capacity of the simulated range.
"""
function run_HZAM_sets_complete_nine_loci(run_name;
	set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
	w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
	S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
	max_generations::Integer = 1500,
	K_total::Integer = 20000,
	cline_width_loci = "mmt",
)
	dir = "HZAM-J_2D_results/$(run_name)_nine_loci"
	dir = "HZAM-J_2D_results/$(run_name)_nine_loci"

	set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
		"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
		"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
		"low_reject_separate_fmt", "high_reject_separate_fmt"]
	total_loci = [18, 36, 27, 27, 27, 18, 18, 27, 27, 27, 27]
	female_mating_trait_loci = [1:9, 10:18, 1:9, 10:18, 10:18, 1:9, 1:9, 10:18, 10:18, 10:18, 10:18]
	male_mating_trait_loci = [1:9, 19:27, 10:18, 1:9, 10:18, 1:9, 1:9, 19:27, 19:27, 1:9, 1:9]
	hybrid_survival_loci = [1:9, 1:9, 1:9, 1:9, 1:9, 1:9, 1:9, 1:9, 1:9, 1:9, 1:9]
	per_reject_cost = [0, 0, 0, 0, 0, 0.01, 0.05, 0.01, 0.05, 0.01, 0.05]

	println("$(nprocs()) threads running")
	@async @distributed for i in set_numbers
		run_HZAM_set(
			set_names[i],
			total_loci[i],
			female_mating_trait_loci[i],
			male_mating_trait_loci[i],
			hybrid_survival_loci[i],
			per_reject_cost[i],
			w_hyb_set_of_run,
			S_AM_set_of_run,
			dir;
			max_generations,
			K_total,
			cline_width_loci,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end

"""
	categorize(outcome)

Classify the simulation outcome into one of six outcome types based on the cline width, bimodality, and population overlap.
"""
function categorize(outcome; categorization_loci = "mmt")::Integer
	# if there is no simulation outcome return 7 (not one of the 6 outcome types)
	if ismissing(outcome)
		return 7
	end

	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}} = [1]

	if categorization_loci == "fmt"
		loci = outcome.sim_params.female_mating_trait_loci
	else
		loci = outcome.sim_params.male_mating_trait_loci
	end

	# return 6 when one of the initial populations is extinct, 5 when the two populations are blended
	if is_one_species_extinct(outcome; categorization_loci)
		return 6
	elseif is_blended(outcome)
		return 5
	end

	# calculate the percentage of the simulated range occupied by both populations
	overlap = Population.calc_species_overlap(
		outcome.population_data.population,
		0.03,
		0.01,
		loci,
	)[1]

	# calculate the hybrid zone bimodality	
	bimodality = DataAnalysis.calc_bimodality(outcome.population_data, loci)

	# bimodality determines if the simulation outcome is overlap, bimodal, or Unimodal
	# overlap determines if the simulation outcome is narrow overlap or broad overlap
	if bimodality > 0.95
		return overlap < 0.3 ? 3 : 4
	elseif bimodality < 0.5
		return 2
	end

	return 1
end

"""
	is_one_species_extinct(outcome::HZAM.DataAnalysis.OutputData; categorization_loci = "mmt")::Bool

Determine if the simulation outcome shows the extinction of one population.
"""
function is_one_species_extinct(outcome; categorization_loci = "mmt")::Bool
	# phenotype_counts stores the number of individuals with each categorization_loci phenotype for each generation
	# num_loci is the number of loci in the trait of interest (as determined by categorization_loci)
	if categorization_loci=="fmt"
		phenotype_counts = outcome.population_tracking_data[1]
		num_loci = length(outcome.sim_params.female_mating_trait_loci)
	else
		phenotype_counts = outcome.population_tracking_data[2]
		num_loci = length(outcome.sim_params.male_mating_trait_loci)
	end

	# phenotype counts at the end of the simulation
	final_phenotype_counts = phenotype_counts[end]

	species_A_proportion = sum(final_phenotype_counts[1:(Integer(floor(0.2*num_loci))+1)])

	species_B_proportion = sum(final_phenotype_counts[(Integer(ceil(1.8*num_loci))+1):length(final_phenotype_counts)])

	return (species_A_proportion == 0 && species_B_proportion > 0.5) || (species_B_proportion == 0 && species_A_proportion>0.5)
end

"""
	is_blended(outcome::HZAM.DataAnalysis.OutputData)::Bool

Determine if a simulation outcome is blended.
"""
function is_blended(outcome)::Bool
	# Simulation outcomes with bimodality>0.95 are not blended. Otherwise they are if the 
	# cline width ever excedes the width of the simulated range.
	if outcome.bimodality > 0.95
		return false
	elseif maximum(outcome.hybrid_zone_width) > 1
		return true
	end

	# Perform a Mann-Kendall test to determine if the cline width is increasing
	cline_widths = outcome.hybrid_zone_width[11:end]
	test = mk_original_test(cline_widths)

	return test.trend=="increasing"
end


"""
	run_HZAM_set(
		set_name::String,
		total_loci::Int,
		female_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		male_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		hybrid_survival_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		per_reject_cost::Real,
		w_hyb_set::Union{UnitRange{<:Real}, Vector{<:Real}},
		S_AM_set::Union{UnitRange{<:Real}, Vector{<:Real}},
		dir::String;
		intrinsic_R::Real = 1.1, K_total::Int = 20000, max_generations::Int = 1500,
		survival_fitness_method::String = "epistasis", sigma_disp = 0.03,
	)

Run the simulation for the given combinations of hybrid fitness and assortative mating 
strength and store the outcome of each simulation in a JLD2 file.

# Arguments
- `set_name::String`: the name assigned to the set of simulations.
- `intrinsic_R::Real`: the intrinsic growth rate.
- `K_total::Integer=20000`: the carrying capacity of the environment.
- `max_generations::Integer=1500`: the number of generations that the simulation will run.
- `total_loci::Integer=6`: the total number of loci in the genome.
- `female_mating_trait_loci=1:3`: the loci specifying the female's mate preference.
- `male_mating_trait_loci=1:3`: the loci specifying the male's mating trait.
- `hybrid_survival_loci=1:3`: the loci specifying the probability of survival to adulthood.
- `survival_fitness_method:String="epistasis"`: the method used to calculate the probability 
of survival to adulthood.
- `per_reject_cost=0`: the fitness loss of female per male rejected (due to search time, 
etc.). Can take values of 0 to 1.
- `sigma_disp=0.05`: the standard deviation of the normal distribution determining how far 
offspring will disperse from their mothers.
"""
function run_HZAM_set(
	set_name::String,
	total_loci::Int,
	female_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	male_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	hybrid_survival_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	per_reject_cost::Real,
	w_hyb_set::Union{UnitRange{<:Real}, Vector{<:Real}},
	S_AM_set::Union{UnitRange{<:Real}, Vector{<:Real}},
	dir::String;
	intrinsic_R::Real = 1.1, K_total::Int = 20000, max_generations::Int = 1500,
	survival_fitness_method::String = "epistasis", sigma_disp = 0.03,
	cline_width_loci = "mmt",
)
	dir = mkpath(string(dir, "/", set_name))

	# Loop through the different simulation sets
	for i in eachindex(w_hyb_set)
		for j in eachindex(S_AM_set)
			w_hyb = w_hyb_set[i]
			S_AM = S_AM_set[j]

			run_name =
				string(
					"HZAM_simulation_run_gen", max_generations, "_SC",
					per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM,
				)

			println("Beginning: $run_name")
			# run one simulation by calling the function defined above:
			outcome = run_one_HZAM_sim(
				w_hyb,
				S_AM,
				intrinsic_R;
				K_total,
				max_generations,
				total_loci,
				female_mating_trait_loci,
				male_mating_trait_loci,
				hybrid_survival_loci,
				survival_fitness_method,
				per_reject_cost,
				sigma_disp,
				do_plot = false,
				run_name = run_name,
				cline_width_loci,
			)
			println("Ending: $run_name")

			if !isnothing(outcome)
				run_name_full_gen =
					string(
						"HZAM_simulation_run_gen", max_generations, "_SC",
						per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM,
					)
				filename = string(dir, "/", run_name_full_gen, ".jld2")
				@save filename outcome
			end
		end # of S_AM loop
	end # of w_hyb loop   
end

function run_HZAM_set_search_cost_supplement(
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
	println(per_reject_cost_set)
	println(S_AM_set)

	index_array = Array{Missing, 2}(
		undef, length(per_reject_cost_set), length(S_AM_set),
	)

	# Loop through the different simulation sets
	@async @distributed for i in CartesianIndices(index_array)
		sc = per_reject_cost_set[i[1]]
		S_AM = S_AM_set[i[2]]

		run_name =
			string(
				"HZAM_simulation_run_gen", max_generations, "_SC",
				sc, "_Whyb", w_hyb, "_SAM", S_AM,
			)

		println("Beginning: $run_name")
		# run one simulation by calling the function defined above:
		outcome = run_one_HZAM_sim(
			w_hyb,
			S_AM,
			intrinsic_R;
			K_total,
			max_generations,
			total_loci,
			female_mating_trait_loci,
			male_mating_trait_loci,
			hybrid_survival_loci,
			survival_fitness_method,
			per_reject_cost = sc,
			sigma_disp,
			do_plot = false,
			run_name = run_name,
			cline_width_loci,
		)
		println("Ending: $run_name")

		if !isnothing(outcome)
			run_name_full_gen =
				string(
					"HZAM_simulation_run_gen", max_generations, "_SC",
					sc, "_Whyb", w_hyb, "_SAM", S_AM,
				)
			filename = string(dir, "/", run_name_full_gen, ".jld2")
			@save filename outcome
		end
	end # of S_AM loop
end


"""
	set_results_folder(dir::String)

Set the output parent directory for all methods in summarize_data.jl.
"""
function set_results_folder(dir::String)
	global results_folder = mkpath(dir)
end


"""
	read_directory(dir::String, set::String)

Load the data from each simulation output style in each run folder in the given directory 
for the given simulation set (i.e. "full_pleiotropy").
"""
function read_directory(dir::String, set::String)
	input_folders = readdir(dir, join = true)
	outcome_arrays = []
	for folder in input_folders
		try
			push!(outcome_arrays, load_from_folder("$(folder)/$(set)"))
		catch
		end
	end
	return outcome_arrays
end

"""
	load_from_folder(dir::String)

Load the data from each simulation output file stored in the given directory into an array 
organized by hybrid fitness and assortative mating strength.
"""
function load_from_folder(dir::String)
	outcome_array = Array{Union{DataAnalysis.OutputData, Missing}, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	for i in eachindex(outcome_array)
		outcome_array[i] = missing
	end
	try
		files = readdir(dir)
		for file in files
			path = string(dir, "/", file)
			if occursin(".jld2", path)
				@load path outcome

				println(outcome.sim_params)

				if outcome.sim_params.S_AM in S_AM_set
					outcome_array[
						indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
						indexin(outcome.sim_params.S_AM, S_AM_set)[1],
					] = outcome
				end
			end
		end
	catch
	end
	return outcome_array
end

"""
	plot_fitnesses(fitnesses::Vector{<:Dict})

Produce a plot of fitnesses per phenotype over time.

# Arguments
- `fitnesses::Vector{<:Dict}`: number of offspring per phenotype each generation.
"""
function plot_fitnesses(fitnesses::Vector{<:Dict})
	dir = mkpath(string(results_folder, "/plots/"))
	@save string(dir, "fitness.jld2") fitnesses
	generations = collect(eachindex(fitnesses))
	num_generations = length(generations)
	phenotypes = collect(keys(fitnesses[1]))
	num_phenotypes = length(fitnesses[1])
	xs = vcat([fill(gen, num_phenotypes) for gen in eachindex(fitnesses)]...)
	ys = vcat([collect(keys(f)) for f in fitnesses]...)
	zs = vcat([collect(values(f)) for f in fitnesses]...)

	"""
	Determines which phenotype has the highest fitness each generation.
	"""
	function maximum_fitness()
		return map(
			f -> reduce((x, y) -> f[x] ≥ f[y] ? x : y, keys(f)),
			fitnesses,
		)
	end

	"""
	Compute the running average of the fitnesses for a given phenotype.
	"""
	function average_fitness(phenotype)
		n = 10
		average_fitnesses = []
		for i in eachindex(fitnesses)
			fitness_values =
				[fitnesses[j][phenotype] for j in i:min(length(fitnesses), i+n-1)]
			push!(average_fitness, mean(fitness_values))
		end
		return average_fitnesses
	end

	fig = Figure(resolution = (1800, 1200), figure_padding = 60)
	ax = Axis(fig[1, 1], xlabel = "generation", ylabel = "phenotype")
	scatter!(ax, xs, ys)
	display(fig)
	display(heatmap(xs, ys, zs, colormap = :grayC))
	readline()
end

"""
	plot_population_tracking_data(filepath::String)

Create plots of the population size, hybrid zone width, hybrid index, and population overlap 
vs time.

# Arguments
- `filepath::String`:: the filepath for the file containing the outcome of the simulation.
"""
function plot_population_tracking_data(filepath::String)
	@load filepath outcome
	population_tracking_data = outcome.population_tracking_data
	sim_params = outcome.sim_params

	overlaps = [t.overlap for t in population_tracking_data]
	hybridity = [t.hybridity for t in population_tracking_data]
	widths = [t.width for t in population_tracking_data]
	populations = [t.population_size for t in population_tracking_data]

	xs = collect(eachindex(population_tracking_data))
	ax = Vector(undef, 4)
	points = Vector(undef, 4)

	fontsize_theme = Theme(fontsize = 40)
	set_theme!(fontsize_theme)  # set the standard font size
	fig = Figure(resolution = (1800, 1200), figure_padding = 60)
	# create the axis and labels
	ax[1] = Axis(
		fig[1, 1],
		ylabel = "# individuals",
		yticklabelsize = 30.0f0,
	)
	ax[2] = Axis(
		fig[2, 1],
		ylabel = "cline width",
		limits = ((nothing, nothing), (-0.05, 1.05)),
	)
	ax[3] = Axis(
		fig[3, 1],
		ylabel = "hybridity",
		yticklabelsize = 30.0f0,
	)
	ax[4] = Axis(
		fig[4, 1],
		ylabel = "overlap",
		xlabel = "generation",
		limits = ((nothing, nothing), (-0.05, 1.05)),
	)
	points[1] = scatter!(ax[1], xs, populations)
	points[2] = scatter!(ax[2], xs, widths)
	points[3] = scatter!(ax[3], xs, hybridity)
	points[4] = scatter!(ax[4], xs, overlaps)

	yspace = 1.5 * maximum(tight_yticklabel_spacing!, [ax[1], ax[2]])

	ax[1].yticklabelspace = yspace
	ax[2].yticklabelspace = yspace
	ax[3].yticklabelspace = yspace
	ax[4].yticklabelspace = yspace

	display(fig)

	readline()

	dir = mkpath(string(results_folder, "/plots/"))

	save(string(
		dir,
		"_S_AM",
		sim_params.S_AM,
		"_w_hyb",
		sim_params.w_hyb,
		".png",
		fig,
	))

	readline()
end

"""
	summarize_gene_correlations(source_dir::String; filename="gene_correlations")

Create plots of the gene correlations between different traits vs S_AM and w_hyb for 
different combinations of the same loci controlling different traits.

source_dir should point to a folder containing magic_preference.jld2, magic_cue.jld2, 
search_cost.jld2, and no_magic.jld2. Each file should contain an outcome array from a 
simulation set.

# Arguments
- `source_dir::String`: the folder where the outcomes of all the simulations are stored.
- `file_name="gene_correlations`: the name of the image file for the plot. 
"""
function summarize_gene_correlations(source_dir::String; filename = "gene_correlations")
	names = ["magic_preference", "magic_cue", "search_cost", "no_magic"]
	xlabelnames = ["Female mating trait", "Male mating trait", "Neutral trait"]
	ylabelnames = ["Magic preference", "Magic cue", "Search cost", "No pleiotropy"]
	xlabels = []
	ylabels = []
	correlations = Matrix(undef, 4, 3)

	magic_loci = [[1, 5], [3, 5], [5, 6], [5, 6]]
	fmt_loci = [[2], 1:2, 1:2, 1:2]
	mmt_loci = [3:4, [4], 3:4, 3:4]
	neutral_loci = [6:9, 6:9, 7:10, 7:10]

	for i in 1:4
		filename = string(source_dir, "/", names[i], ".jld2")
		@load filename sim_params outcome_array
		correlations[i, 1] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					fmt_loci[i],
				),
				outcome_array,
			)...)

		correlations[i, 2] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					mmt_loci[i],
				),
				outcome_array,
			)...)

		correlations[i, 3] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					neutral_loci[i],
				),
				outcome_array,
			)...)

	end

	w_hyb_set = string.([0.95, 0.9, 0.7, 0.5])
	S_AM_set = string.([1, 10, 100, 1000])
	xticks = [4, 3, 2, 1]
	yticks = [1, 2, 3, 4]

	xs = vcat(fill(xticks, 4)...)
	ys = vcat([fill(s, 4) for s in yticks]...)

	ax = Matrix(undef, 4, 4)
	hm = Matrix(undef, 4, 4)


	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)  # this sets the standard font size

	fig = Figure(
		resolution = (1800, 1200),
		figure_padding = 60,
		colorrange = (0, 1),
		colormap = :curl,
	)

	println(correlations)

	for j in 1:4
		for i in 1:3
			ax[j, i] = Axis(
				fig[j, i],
				xlabel = "w_hyb",
				ylabel = "S_AM",
				xticks = (xticks, w_hyb_set),
				yticks = (yticks, S_AM_set),
				xticklabelsize = 20,
				yticklabelsize = 20,
				xlabelsize = 15,
				ylabelsize = 15,
			)
			hm[j, i] = heatmap!(ax[j, i], xs, ys, correlations[j, i], colorrange = (-1, 1))
			println(correlations[j, i][16])
		end
		println("")
	end

	for i in 1:3
		push!(xlabels, Label(fig[0, i], xlabelnames[i], tellwidth = false))
	end

	for i in 1:4
		push!(
			ylabels,
			Label(fig[i, 0], ylabelnames[i], rotation = pi / 2, tellheight = false),
		)
	end

	Colorbar(fig[:, 4], hm[1, 1], label = "Pearson coefficient", ticklabelsize = 15)
	display(fig)
	dir = mkpath(string(results_folder, "/plots/"))
	save(string(dir, "/", filename, ".png"), fig)
	readline()
end

function calc_extinction_frequencies()
	extinctions_array = Array{Integer, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)

	for j in eachindex(extinctions_array)
		extinctions_array[j]=0
	end

	for i in 1:3
		source_dir = "HZAM-J_2D_results_categorized/mmt_three_loci/Run3_three_loci_$i"
		files = readdir(source_dir)

		for file in files
			path = string(source_dir, "/", file)
			@load path outcome_array
			for j in eachindex(extinctions_array)
				if ismissing(outcome_array[j])
					continue
				end
				categorized_outcome = Integer(outcome_array[j])
				if outcome_array[j] == 6
					extinctions_array[j] += 1
				end
			end
		end
	end

	return extinctions_array
end

function calc_overlap_difference(folder1, folder2)
	overlap_outcomes1 = 0

	overlap_outcomes2 = 0

	for i in 1:3
		source_dir = "HZAM-J_2D_results_categorized/mmt_three_loci/Run3_three_loci_$i"
		files = readdir(source_dir)

		@load string(source_dir, "/", folder1) outcome_array
		outcome_array1 = outcome_array

		@load string(source_dir, "/", folder2) outcome_array
		outcome_array2 = outcome_array

		overlap_outcomes1 += count(x->x == 4, outcome_array1)
		overlap_outcomes2 += count(x->x == 4, outcome_array2)


		@load string(source_dir, "/", folder1) parameters_array

		for j in eachindex(parameters_array)
			if parameters_array[j].S_AM == 300
				if outcome_array1[j]==4
					overlap_outcomes1 += 1
				end

				if outcome_array2[j] == 4
					overlap_outcomes2 += 1
				end
			end
		end
	end

	return overlap_outcomes1 / overlap_outcomes2
end

function calc_outcome_frequencies(outcome; S_AM_range = 1:8, w_hyb_range = 1:13)
	folders = [
		"full_pleiotropy",
		"separate_hst",
		"separate_mmt",
		"separate_fmt",
		"low_reject_separate_fmt",
		"high_reject_separate_fmt",
		"no_pleiotropy",
		"low_reject_no_pleiotropy",
		"high_reject_no_pleiotropy",
	]

	function calc_outcome_frequency(folder)
		outcomes = 0
		extinction_outcomes = 0
		for i in 1:3
			source_dir = "HZAM-J_2D_results_categorized/mmt_three_loci/Run3_three_loci_$i"
			@load "$source_dir/$folder.JLD2" outcome_array
			@load "$source_dir/$folder.JLD2" parameters_array
			outcome_array = outcome_array[w_hyb_range, S_AM_range]
			parameters_array = parameters_array[w_hyb_range, S_AM_range]
			outcomes += length(outcome_array)
			extinction_outcomes += count(x->x==outcome, outcome_array)
		end

		return extinction_outcomes / outcomes
	end

	return Dict(zip(folders, map(x->calc_outcome_frequency(x), folders)))
end

"""
	extract_data(source_folder::String, source_setup::String, output_folder::String, categorization_loci::String)

Load all of the simulation outcomes in a folder and condense the data by removing anything not used in the paper.

# Arguments
- `source_folder::String`: the folder where the simulation outcomes are stored
- `source_setup::String`: the simulation setup of interest (i.e. "no_pleiotropy)
- `output_folder::String`: where the output is to be stored
- `categorization_loci::String`: whether mmt or fmt is used for categorizing outcomes
"""
function extract_data(source_folder::String, source_setup::String, output_folder::String, categorization_loci::String)
	source_dir = string(source_folder, "/", source_setup)
	files=[]
	try
		files = readdir(source_dir)
	catch
		return
	end
	outcome_array = Array{Integer, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	parameters_array = Array{Union{DataAnalysis.SimParams, Missing}, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	overlap_array = Array{Union{<:Real, Missing}, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	bimodality_array = Array{Union{<:Real, Missing}, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	for i in eachindex(outcome_array)
		outcome_array[i] = 7
		parameters_array[i] = missing
		overlap_array[i] = missing
		bimodality_array[i] = missing
	end
	for file in files
		path = string(source_dir, "/", file)
		if occursin(".jld2", path)
			@load path outcome

			if categorization_loci == "fmt"
				loci = outcome.sim_params.female_mating_trait_loci
			else
				loci = outcome.sim_params.male_mating_trait_loci
			end

			outcome_array[
				indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = categorize(outcome; categorization_loci)
			parameters_array[
				indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = outcome.sim_params
			overlap_array[
				indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = HZAM.Population.calc_species_overlap(
				outcome.population_data.population,
				0.03,
				0.01,
				loci,
			)[1]
			bimodality_array[
				indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = HZAM.DataAnalysis.calc_bimodality(outcome.population_data, loci)
		end
	end
	results_folder = mkpath(output_folder)
	@save string(results_folder, "/", source_setup, ".JLD2") outcome_array parameters_array overlap_array bimodality_array
end

"""
	extract_all(source_dir::String, output_dir::String; set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9, categorization_loci = "mmt")

Condense all of the simulation results data into a format easier to work with.

# Arguments
- `source_dir::String`: the directory storing the simulation outcomes
- `output_dir::String`: the output directory to be created
- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9`: the simulation setups used
- `categorization_loci = "mmt"`: the loci used to categorized simulation outcomes (either "fmt" or "mmt")
"""
function extract_all(source_dir::String, output_dir::String; set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9, categorization_loci = "mmt")
	input_folders = readdir(source_dir)
	for folder in input_folders
		for set in set_numbers
			try
				extract_data("$(source_dir)/$folder", set_names[set], "$(output_dir)/$folder", categorization_loci)
			catch
			end
		end
	end
end

"""
	first_to_three(source_dir::String, output_dir::String; set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9, categorization_loci = "mmt")

Check if any parameter setups have not resulted in three of the same outcome type and run another round of simulations for those parameter setups.

# Arguments
- `source_dir::String`: the directory storing the simulation outcomes
- `output_dir::String`: the output directory to be created
- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9`: the simulation setups used
- `categorization_loci = "mmt"`: the loci used to categorized simulation outcomes (either "fmt" or "mmt")
"""
function first_to_three(source_dir::String, output_dir::String; set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9, categorization_loci = "mmt")
	completed = true

	w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

	S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]

	function first_to_three(outcomes)
		for i in 1:6
			occurences = count(x->x==i, outcomes)
			if occurences ≥ 3
				return i
			end
		end
		return 7
	end

	for set_number in set_numbers
		println(set_names[set_number])
		outcome_arrays = HZAM.read_directory(source_dir, set_names[set_number])

		categorized_outcome_arrays = [categorize.(outcome_array; categorization_loci) for outcome_array in outcome_arrays]


		first_to_three_outcome_array = first_to_three.([getindex.(categorized_outcome_arrays, i) for i in CartesianIndices(outcome_arrays[1])])

		indices = findall(x->x==7, first_to_three_outcome_array)

		println(indices)
		if length(indices) > 0
			completed = false
		end

		results_folder = mkpath("$(output_dir)/$(set_names[set_number])")
		@sync @distributed for i in indices
			w_hyb = w_hyb_set[i[1]]
			S_AM = S_AM_set[i[2]]

			sample = outcome_arrays[1][i]

			run_name =
				string(
					"HZAM_simulation_run_gen", sample.sim_params.max_generations, "_SC",
					sample.sim_params.per_reject_cost, "_Whyb", sample.sim_params.w_hyb,
					"_SAM", sample.sim_params.S_AM,
				)


			println("Beginning: $run_name")
			# run one simulation by calling the function defined above:
			outcome = HZAM.run_one_HZAM_sim(
				w_hyb,
				S_AM,
				sample.sim_params.intrinsic_R;
				K_total = sample.sim_params.K_total,
				max_generations = sample.sim_params.max_generations,
				total_loci = sample.sim_params.total_loci,
				female_mating_trait_loci = sample.sim_params.female_mating_trait_loci,
				male_mating_trait_loci = sample.sim_params.male_mating_trait_loci,
				hybrid_survival_loci = sample.sim_params.hybrid_survival_loci,
				per_reject_cost = sample.sim_params.per_reject_cost,
				sigma_disp = sample.sim_params.sigma_disp,
				do_plot = false,
				cline_width_loci,
			)
			println("Ending: $run_name")

			if !isnothing(outcome)
				run_name_full_gen =
					string(
						"HZAM_simulation_run_gen", sample.sim_params.max_generations, "_SC",
						sample.sim_params.per_reject_cost, "_Whyb", sample.sim_params.w_hyb,
						"_SAM", sample.sim_params.S_AM,
					)
				filename = string(results_folder, "/", run_name_full_gen, ".jld2")
				@save filename outcome
			end
		end
	end

	return completed
end
