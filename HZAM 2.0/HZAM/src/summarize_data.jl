using Colors, ColorSchemes
using GLMakie
using JLD2 # needed for saving / loading data in Julia format
#using CSV # for saving in csv format
#using DataFrames # for converting data to save as CSV
using LsqFit: curve_fit

mean(itr) = sum(itr) / length(itr)

"The directory where all files are saved"
global results_folder = mkpath("HZAM-J_2D_results")
global working_dir = results_folder


"The set of hybrid fitnesses (w_hyb) values that will be run"
global w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]


"The set of assortative mating strengths (S_AM) values that will be run"
global S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

global spaced_locations = collect(Float32, 0:0.01:1)

struct OverlapData
	cline_width::Real
	overlap::Real
	bimodality::Real
end

function plot_num_loci(outcome_array)
	mmt_loci_count = Integer[]
	fmt_loci_count = Integer[]
	correlations = Real[]
	function categorize(outcome)

		pd = outcome.population_data

		genotypes = [
			vcat([d.genotypes_F for d in pd.population]...)
			vcat([d.genotypes_M for d in pd.population]...)
		]
		x_locations = [
			vcat([d.x_locations_F for d in pd.population]...)
			vcat([d.x_locations_M for d in pd.population]...)
		]
		y_locations = [
			vcat([d.y_locations_F for d in pd.population]...)
			vcat([d.y_locations_M for d in pd.population]...)
		]


		mmt_loci = outcome.sim_params.male_mating_trait_loci
		fmt_loci = outcome.sim_params.female_mating_trait_loci
		hst_loci = outcome.sim_params.hybrid_survival_loci

		sigmoid_curves, cline_widths = DataAnalysis.calc_sigmoid_curves(
			x_locations,
			y_locations,
			DataAnalysis.calc_traits_additive(genotypes, fmt_loci),
		)
		fmt_cline_width = mean(cline_widths)

		sigmoid_curves, cline_widths = DataAnalysis.calc_sigmoid_curves(
			x_locations,
			y_locations,
			DataAnalysis.calc_traits_additive(genotypes, mmt_loci),
		)
		mmt_cline_width = mean(cline_widths)
		sigmoid_curves, cline_widths = DataAnalysis.calc_sigmoid_curves(
			x_locations,
			y_locations,
			DataAnalysis.calc_traits_additive(genotypes, hst_loci),
		)
		hst_cline_width = mean(cline_widths)
		println("mmt_loci=$mmt_loci and fmt_loci=$fmt_loci")
		println((fmt_cline_width, mmt_cline_width, hst_cline_width))

		correlation = DataAnalysis.calc_trait_correlation(genotypes, mmt_loci, hst_loci)
		if fmt_cline_width < hst_cline_width * 0.8
			output = -1
		elseif hst_cline_width * 0.8 < fmt_cline_width < hst_cline_width * 1.2
			output = 0
		elseif correlation > 0.9
			output = 1
		else
			output = 2
		end


		return length(mmt_loci), length(fmt_loci), output
	end

	for o in outcome_array
		mmt_count, fmt_count, correlation = categorize(o)
		push!(mmt_loci_count, mmt_count)
		push!(fmt_loci_count, fmt_count)
		push!(correlations, correlation)
	end


	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colormap = :grayC)

	ax = Axis(
		fig[1, 1],
		xlabel = "mmt_loci",
		ylabel = "fmt_loci",
	)
	hm = heatmap!(ax, mmt_loci_count, fmt_loci_count, correlations)
	display(fig)
end


function run_HZAM_sets_complete(trial_name::String;
	set_numbers = 1:9,
	w_hyb_set_of_run = w_hyb_set,
	S_AM_set_of_run = S_AM_set,
	max_generations = 2000,
	K_total = 30000,
)
	set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
		"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
		"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
		"one_locus_full_pleiotropy", "one_locus_no_pleiotropy"]
	total_loci = [6, 12, 9, 9, 9, 6, 6, 9, 9, 1, 3]
	female_mating_trait_loci = [1:3, 4:6, 1:3, 4:6, 4:6, 1:3, 1:3, 4:6, 4:6, [1], [1]]
	male_mating_trait_loci = [1:3, 7:9, 4:6, 1:3, 4:6, 1:3, 1:3, 7:9, 7:9, [1], [2]]
	hybrid_survival_loci = [1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, [1], [3]]
	per_reject_cost = [0, 0, 0, 0, 0, 0.01, 0.05, 0.01, 0.05, 0, 0]

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
			S_AM_set_of_run;
			max_generations,
			K_total,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end

function run_HZAM_set_num_loci(
	set_name::String,
	per_reject_cost;
	intrinsic_R::Real = 1.1, K_total::Int = 30000, max_generations::Int = 2000,
	sigma_disp = 0.03, w_hyb = 0.95, S_AM = 10)
	dir = mkpath(string(results_folder, "/", set_name))
	# Loop through the different simulation sets
	for i in 1:4
		for j in 1:4
			run_name = string("HZAM_simulation_run", "_mmt_loci", i, "_fmt_loci", j, "_Whyb", w_hyb, "_SAM", S_AM)
			println("Beginning: $run_name")
			# run one simulation by calling the function defined above:
			outcome = run_one_HZAM_sim(
				w_hyb,
				S_AM,
				intrinsic_R;
				K_total,
				max_generations,
				total_loci = 11,
				female_mating_trait_loci = 8:(7+j),
				male_mating_trait_loci = 4:(3+i),
				hybrid_survival_loci = 1:3,
				per_reject_cost,
				sigma_disp,
				do_plot = true,
				plot_int = 500,
			)
			println("Ending: $run_name")

			filename = string(dir, "/", run_name, ".jld2")
			@save filename outcome
		end # of S_AM loop
	end # of w_hyb loop   
end


"""
	run_HZAM_set(
		set_name::String,
		total_loci::Int=3,
		female_mating_trait_loci=1:3,
		male_mating_trait_loci=1:3,
		hybrid_survival_loci=1:3,
		per_reject_cost::Real=0; 
		<keyword arguments>
	)

Run the simulation for 64 combinations of hybrid fitness and assortative mating strength and 
store the outcome of each simulation in a JLD2 file.

# Arguments
- `set_name::String`: the name assigned to the set of simulations.
- `intrinsic_R::Real`: the intrinsic growth rate.
- `K_total::Integer=40000`: the carrying capacity of the environment.
- `max_generations::Integer=1000`: the number of generations that the simulation will run for.
- `total_loci::Integer=6`: the total number of loci in the genome.
- `female_mating_trait_loci=1:3`: the loci specifying the female's mate preference.
- `male_mating_trait_loci=1:3`: the loci specifying the male's mating trait.
- `hybrid_survival_loci=1:3`: the loci specifying the probability of survival to adulthood.
- `survival_fitness_method:String="epistasis"`: the method used to calculate the probability of survival to adulthood.
- `per_reject_cost=0`: the fitness loss of female per male rejected (due to search time, etc.). Can take values of 0 to 1.
- `sigma_disp=0.05`: the standard deviation of the normal distribution determining how far offspring will disperse from their mothers.
"""
function run_HZAM_set(
	set_name::String,
	total_loci::Int,
	female_mating_trait_loci,
	male_mating_trait_loci,
	hybrid_survival_loci,
	per_reject_cost,
	w_hyb_set,
	S_AM_set;
	intrinsic_R::Real = 1.1, K_total::Int = 30000, max_generations::Int = 2000,
	survival_fitness_method::String = "epistasis", sigma_disp = 0.03,
)
	dir = mkpath(string(results_folder, "/", set_name))
	set_working_dir(dir)

	# Loop through the different simulation sets
	for i in eachindex(w_hyb_set)
		for j in eachindex(S_AM_set)
			w_hyb = w_hyb_set[i]
			S_AM = S_AM_set[j]

			run_name_half_gen =
				string(
					"HZAM_simulation_run_gen", trunc(max_generations / 2), "_SC",
					per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM,
				)

			println("Beginning: $run_name_half_gen")
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
				run_name = run_name_half_gen,
			)
			println("Ending: $run_name_half_gen")

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

function search_cost_clines()
	search_costs = [0.01, 0.05, 0.1]
	for s in search_costs
		run_HZAM_set(
			"search_cost_cline_sc_$s",
			1,
			[1],
			[1],
			[1],
			s,
			[1, 1, 1, 1, 1, 1],
			[3, 10, 30, 100, 300, 1000]
		)
	end
end

"""
	set_results_folder(dir::String)

Set the output parent directory for all methods in summarize_data.jl.
"""
function set_results_folder(dir::String)
	global results_folder = mkpath(dir)
end

"""
	set_working_dir(dir::String)

Set the output directory for all methods in summarize_data.jl.
"""
function set_working_dir(dir::String)
	global working_dir = mkpath(dir)
end

"""
	plot_output_field(
		outcomes::Array{<:Real},
		sim_params::Array{DataAnalysis.SimParams}
	)

Create a heatmap of an output variable vs hybrid fitness and assortative mating.

# Arguments
- `outcomes::Array{:Real}`: the output from the simulation to be displayed.
- `sim_params::Array{<:DataAnalysis.SimParams}`: the simulation parameters resulting in the outcomes.
"""
function plot_output_field(
	outcomes::Array{<:Real},
	sim_params::Array{DataAnalysis.SimParams},
)
	output = [outcomes...]
	w_hybs = [[s.w_hyb for s in sim_params]...]
	S_AMs = [[s.S_AM for s in sim_params]...]

	w_hyb_set = sort(union(w_hybs))
	S_AM_set = sort(union(S_AMs))
	xticks = collect(1:length(S_AM_set))
	yticks = collect(1:length(w_hyb_set))

	xs = map(s -> indexin(s, S_AM_set)[1], S_AMs)
	ys = map(w -> indexin(w, w_hyb_set)[1], w_hybs)


	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)  # this sets the standard font size

	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colormap = :grayC)

	ax = Axis(
		fig[1, 1],
		xlabel = "w_hyb",
		ylabel = "S_AM",
		xticks = (xticks, string.(S_AM_set)),
		yticks = (yticks, string.(w_hyb_set)),
		xticklabelsize = 20,
		yticklabelsize = 20,
		xlabelsize = 15,
		ylabelsize = 15,
	)
	hm = heatmap!(ax, xs, ys, output, colorrange = (0, 0.5))

	Colorbar(fig[:, 2], hm, ticklabelsize = 15)
	display(fig)
	readline()
	return fig
end

"""
	load_from_folder(dir::String)

Load the data from each simulation output file stored in the given directory into an array 
organized by hybrid fitness and assortative mating strength.
"""
function load_from_folder(dir::String)
	files = readdir(dir)
	outcome_array = Array{DataAnalysis.OutputData, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	for file in files
		path = string(dir, "/", file)
		if occursin(".jld2", path)
			@load path outcome
			println(outcome.sim_params)
			outcome_array[
				indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
				indexin(outcome.sim_params.S_AM, S_AM_set)[1],
			] = outcome
		end
	end
	return outcome_array
end


function load_from_folder_num_loci(dir::String)
	files = readdir(dir)
	outcome_array = Array{DataAnalysis.OutputData, 2}(
		undef, 4, 4,
	)
	for file in files
		path = string(dir, "/", file)
		if occursin(".jld2", path)
			@load path outcome
			println(outcome[1].sim_params)
			outcome_array[
				length(outcome[1].sim_params.male_mating_trait_loci),
				length(outcome[1].sim_params.female_mating_trait_loci),
			] = outcome[1]
		end
	end
	return outcome_array
end


function load_overlap_data_from_folder(dir::String)
	set_names = readdir(dir)
	outcome_arrays = Array{OverlapData, 2}[]

	for set in set_names
		set_dir = string(dir, "/", set)
		if isdir(set_dir)
			outcome_array = Array{OverlapData, 2}(
				undef, length(w_hyb_set), length(S_AM_set),
			)
			files = readdir(set_dir)

			for file in files
				path = string(set_dir, "/", file)
				if occursin(".jld2", path)
					@load path outcome
					outcome_array[
						indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
						indexin(outcome.sim_params.S_AM, S_AM_set)[1],
					] = OverlapData(outcome.hybrid_zone_width, outcome.population_overlap, outcome.bimodality)
				end
			end

			push!(outcome_arrays, outcome_array)
		end
	end

	return Dict(zip(set_names, outcome_arrays))
end








#=
"""
	load_from_csv(filepath::String)

Read the data from a CSV file and return vectors for hybrid fitness, assortative mating 
strength, and the output variable.
"""
function load_from_csv(filepath::String)
	df = DataFrame(CSV.File(filepath))

	return df[!, "w_hyb"], df[!, "S_AM"], df[:, 3]
end

"""
	convert_to_CSV(
		outcome_array::Array{<:Real},
		w_hyb_array::Array{<:Real},
		S_AM_array::Array{<:Real},
		name::String
	)

Convert an array of outcomes to a csv file for the given output field.

# Arguments
- `outcome_array::Array{<:Real}`: the array of the output data of interest.
- `w_hyb_array::Array{<:Real}`: the array of the w_hyb parameters used.
- `S_AM_array::Array{<:Real}`: the array of the S_AM parameters used.
- `field_name::Symbol`: the name of the output field of interest.
- `name::String`: the name for the CSV.
"""
function convert_to_CSV(
	outcome_array::Array{<:Real},
	w_hyb_array::Array{<:Real},
	S_AM_array::Array{<:Real},
	name::String
)
	output = [outcome_array...]
	w_hybs = [w_hyb_array...]
	S_AMs = [S_AM_array...]

	df = DataFrame(S_AM=S_AMs, w_hyb=w_hybs, output_field=output)

	dir = mkpath(string(results_folder, "/simulation_outcomes/CSV_data/"))

	CSV.write(string(dir, "/", name, ".csv"), df)
end
=#

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
				[fitnesses[j][phenotype] for j in i:min(length(fitnesses), i + n - 1)]
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

Create plots of the population size, hybrid zone width, hybrid index, and population overlap vs time.

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
search_cost.jld2, and no_magic.jld2. Each file should contain an outcome array from a simulation set.

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
		correlations[i, 1] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], fmt_loci[i]), outcome_array)...)
		correlations[i, 2] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], mmt_loci[i]), outcome_array)...)

		correlations[i, 3] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], neutral_loci[i]), outcome_array)...)

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

	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colorrange = (0, 1), colormap = :curl)

	println(correlations)

	for j in 1:4
		for i in 1:3
			ax[j, i] = Axis(fig[j, i], xlabel = "w_hyb", ylabel = "S_AM", xticks = (xticks, w_hyb_set), yticks = (yticks, S_AM_set), xticklabelsize = 20, yticklabelsize = 20, xlabelsize = 15, ylabelsize = 15)
			hm[j, i] = heatmap!(ax[j, i], xs, ys, correlations[j, i], colorrange = (-1, 1))
			println(correlations[j, i][16])
		end
		println("")
	end

	for i in 1:3
		push!(xlabels, Label(fig[0, i], xlabelnames[i], tellwidth = false))
	end

	for i in 1:4
		push!(ylabels, Label(fig[i, 0], ylabelnames[i], rotation = pi / 2, tellheight = false))
	end

	Colorbar(fig[:, 4], hm[1, 1], label = "Pearson coefficient", ticklabelsize = 15)
	display(fig)
	dir = mkpath(string(results_folder, "/plots/"))
	save(string(dir, "/", filename, ".png"), fig)
	readline()
end

function compare_phenotype_plots()
	f = Figure(resolution = (1800, 1200))

	ga = f[1, 1] = GridLayout()
	gb = f[1, 2] = GridLayout()
	gc = f[2, 1] = GridLayout()
	gd = f[2, 2] = GridLayout()

	make_phenotype_subplot(ga, "phenotype_counts/phenotype_counts_population_data_no_pleiotropy", "No Pleiotropy")
	f[1, 2] = make_phenotype_subplot(gb, "phenotype_counts/phenotype_counts_population_data_magic_cue", "Magic Cue")
	f[2, 1] = make_phenotype_subplot(gc, "phenotype_counts/phenotype_counts_population_data_magic_preference", "Magic Preference")
	f[2, 2] = make_phenotype_subplot(gd, "phenotype_counts/phenotype_counts_population_data_no_pleiotropy_high_search_cost", "No Pleiotropy, 5% Search Cost")

	display(f)
	save("phenotypes_over_time2.png", f)
end

function plot_transect(genotypes, x_locations, y_locations, loci, g)
	function sigmoid(x::Vector{<:Real}, p::Vector{<:Real})
		1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
	end
	hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, loci)
	sigma_comp = 0.01
	sorted_indices = DataAnalysis.sort_y(y_locations)[10]

	fit = curve_fit(sigmoid, x_locations[sorted_indices], hybrid_indices[sorted_indices], [0.0, 1.0])
	# compute the values of the sigmoid curve at evenly spaced x values
	sigmoid_curve = sigmoid(spaced_locations, fit.param)

	println(sigmoid_curve)
	readline()

	ax = Axis(g)

	lines!(ax, spaced_locations, sigmoid_curve)
end

function make_phenotype_subplot(g, filename, plotname)
	colgap!(g, 10)
	rowgap!(g, 10)

	function reshape_output(phenotypes)
		output = Matrix(undef, length(phenotypes), length(phenotypes[1]))

		for i in eachindex(phenotypes)
			for j in eachindex(phenotypes[1])
				output[i, j] = phenotypes[i][j]
			end
		end
		output
	end

	@load filename mmt_phenotype_counts fmt_phenotype_counts

	Label(g[1, 1, Top()], plotname, valign = :bottom,
		font = :bold,
		padding = (0, 0, 5, 0),
		fontsize = 32)

	ax1 = Axis(
		g[1, 1],
		ylabel = "Cue",
		yticks = ([1, 7], ["0", "1"]),
		ylabelsize = 24,
	)

	ax2 = Axis(
		g[2, 1],
		xlabel = "generation",
		ylabel = "Preference",
		yticks = ([1, 7], ["0", "1"]),
		xlabelsize = 24,
		ylabelsize = 24,
	)




	hidexdecorations!(ax1, ticks = false)

	linkxaxes!(ax1, ax2)

	heatmap!(ax1, reshape_output(mmt_phenotype_counts), colormap = Reverse(:grayC))

	heatmap!(ax2, reshape_output(fmt_phenotype_counts), colormap = Reverse(:grayC))

	g
end

function methods_plot()
	@load "phenotype_counts_population_data_magic_preference" outcome
	println("working2")

	pd = outcome.population_data

	genotypes = [
		vcat([d.genotypes_F for d in pd.population]...)
		vcat([d.genotypes_M for d in pd.population]...)
	]
	x_locations = [
		vcat([d.x_locations_F for d in pd.population]...)
		vcat([d.x_locations_M for d in pd.population]...)
	]
	y_locations = [
		vcat([d.y_locations_F for d in pd.population]...)
		vcat([d.y_locations_M for d in pd.population]...)
	]

	hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, outcome.sim_params.male_mating_trait_loci)

	fig = Figure(resolution = (1800, 1200))
	# create the axis and labels
	ax = Axis(
		fig[1, 1],
		xlabel = "x",
		ylabel = "y",
		xlabelsize = 25,
		ylabelsize = 25,
		xticklabelsize = 25,
		yticklabelsize = 25,
	)

	g = fig[1, 2] = GridLayout()
	# set the limits for the plotted area
	xlims!(-0.03, 1.03)
	ylims!(-0.03, 1.03)

	# add the location of every individual to the plot
	points = scatter!(
		ax,
		x_locations,
		y_locations,
		color = hybrid_indices,
		markersize = 5,
	)


	Colorbar(fig[:, 0], points, ticklabelsize = 25, flipaxis = false, label = "Male mating trait", labelsize = 25)

	y_coords = collect(0.1:0.2:0.9)
	hybrid_index_avgs, sigmoid_curves, widths = DataAnalysis.calc_transects(
		pd,
		outcome.sim_params.male_mating_trait_loci,
		y_coords,
	)

	sigmoid_axis = []

	for i in eachindex(sigmoid_curves)
		push!(
			sigmoid_axis,
			Axis(
				g[i, 1],
				xlabel = "x",
				xlabelsize = 25,
				xticklabelsize = 25,
			),
		)

		scatter!(
			sigmoid_axis[i],
			spaced_locations,
			hybrid_index_avgs[i],
			markersize = 4,
		)
		lines!(
			sigmoid_axis[i],
			spaced_locations,
			sigmoid_curves[i],
			linewidth = 10,
			color = (:gray, 0.5))
		text!(
			Point.(0.1, 0.8),
			text = "y = $(y_coords[i])",
			align = (:center, :center),
			color = :black,
			fontsize = 25,
		)
		hideydecorations!(sigmoid_axis[i])
		if y_coords[i] != y_coords[end]
			hidexdecorations!(sigmoid_axis[i])
		end
	end
	rowgap!(g, 0)

	display(fig)

	save(string("methods_fig.png"), fig)
end

function summarize_phenotypes(mmt_phenotypes, fmt_phenotypes)
	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colormap = :grayC)

	function reshape_output(phenotypes)
		output = Matrix(undef, length(phenotypes), length(phenotypes[1]))

		for i in eachindex(phenotypes)
			for j in eachindex(phenotypes[1])
				output[i, j] = phenotypes[i][j]
			end
		end
		output
	end

	ax1 = Axis(fig[1, :], xlabel = "generation", ylabel = "Mating cue")

	ax2 = Axis(fig[2, :], xlabel = "generation", ylabel = "Mating preference")
	heatmap!(ax1, reshape_output(mmt_phenotypes), colormap = Reverse(:grayC))

	display(fig)
	readline()

	heatmap!(ax2, reshape_output(fmt_phenotypes), colormap = Reverse(:grayC))
	display(fig)
	readline()

end

function plot_hybrid_index(hybrid_indices)
	ys = sort(hybrid_indices)
	xs = eachindex(ys)
	f = Figure()

	Axis(f[1, 1])

	lines!(xs, ys)

	display(f)

	readline()
end














function make_subplot(outcome_array, plot_title)
	sim_params = [o.sim_params for o in outcome_array]

	function assign_output_value(outcome)
		overlap_area, hybrid_area = Population.calc_species_overlap(
			outcome.population_data.population,
			0.03,
			0.01,
			outcome.sim_params.male_mating_trait_loci,
		)
		if overlap_area > hybrid_area && overlap_area > 0.1
			return 0 - overlap_area
		else
			return outcome.hybrid_zone_width
		end
	end

	output = [assign_output_value(outcome) for outcome in outcome_array]

	println(output)

	w_hyb_set = sort([1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])

	w_hybs = [[s.w_hyb for s in sim_params]...]
	S_AMs = [[s.S_AM for s in sim_params]...]
	xticks = collect(1:length(S_AM_set))
	yticks = collect(1:length(w_hyb_set))

	xs = map(s -> indexin(s, S_AM_set)[1], S_AMs)
	ys = map(w -> indexin(w, w_hyb_set)[1], w_hybs)
	output = [output...]

	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)  # this sets the standard font size

	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colormap = :balance)

	ax = Axis(
		fig[1, 1],
		title = plot_title,
		xlabel = "S_AM",
		ylabel = "w_hyb",
		xticks = (xticks, string.(S_AM_set)),
		yticks = (yticks, string.(w_hyb_set)),
		xticklabelsize = 20,
		yticklabelsize = 20,
		xlabelsize = 15,
		ylabelsize = 15,
		aspect = 1,
	)
	hm = heatmap!(ax, xs, ys, output, colorrange = (-0.5, 0.5))

	Colorbar(fig[:, 2], hm, ticklabelsize = 15)
	display(fig)
	save(string(plot_title, ".png"), fig)
	display(fig)
	return fig
end


function summarize_overlap_and_cline_width(source_dir::String; filename = "gene_correlations")
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
		correlations[i, 1] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], fmt_loci[i]), outcome_array)...)
		correlations[i, 2] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], mmt_loci[i]), outcome_array)...)

		correlations[i, 3] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], neutral_loci[i]), outcome_array)...)
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

	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colorrange = (0, 1), colormap = :curl)

	println(correlations)

	for j in 1:4
		for i in 1:3
			ax[j, i] = Axis(fig[j, i], xlabel = "w_hyb", ylabel = "S_AM", xticks = (xticks, w_hyb_set), yticks = (yticks, S_AM_set), xticklabelsize = 20, yticklabelsize = 20, xlabelsize = 15, ylabelsize = 15)
			hm[j, i] = heatmap!(ax[j, i], xs, ys, correlations[j, i], colorrange = (-1, 1))
			println(correlations[j, i][16])
		end
		println("")
	end

	for i in 1:3
		push!(xlabels, Label(fig[0, i], xlabelnames[i], tellwidth = false))
	end

	for i in 1:4
		push!(ylabels, Label(fig[i, 0], ylabelnames[i], rotation = pi / 2, tellheight = false))
	end

	Colorbar(fig[:, 4], hm[1, 1], label = "Pearson coefficient", ticklabelsize = 15)
	display(fig)
	dir = mkpath(string(results_folder, "/plots/"))
	save(string(dir, "/", filename, ".png"), fig)
	readline()
end



function plot_density(genotypes, locations, sigma_comp, title)
	hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, collect(1:3))

	locations_x = DataAnalysis.calc_distances_to_middle(genotypes, locations, 0.05, collect(1:3)) .+ Ref(0.5)

	function calc_density(focal_location, locations_x, sigma_comp)
		return sum(exp.(-(Ref(focal_location) .- locations_x) .^ 2 ./ Ref(2 * (sigma_comp^2))))
	end

	species_A_locations = locations_x[Bool[h == 0 for h in hybrid_indices]]
	species_B_locations = locations_x[Bool[h == 1 for h in hybrid_indices]]


	densities_A = calc_density.(spaced_locations, Ref(species_A_locations), Ref(sigma_comp))
	densities_B = calc_density.(spaced_locations, Ref(species_B_locations), Ref(sigma_comp))
	densities_other = calc_density.(spaced_locations, Ref(setdiff(locations_x, union(species_A_locations, species_B_locations))), Ref(sigma_comp))

	normalizing_factor = (0.001 * (sum(densities_A) + sum(densities_B)))^(-1)

	densities_A = densities_A .* Ref(normalizing_factor)
	densities_B = densities_B .* Ref(normalizing_factor)
	densities_other = densities_other .* Ref(normalizing_factor)
	f = Figure()

	Axis(f[1, 1])

	lines!(spaced_locations, densities_A)
	lines!(spaced_locations, densities_B)
	lines!(spaced_locations, densities_other)
	display(f)
	save(string(title, ".png"), f)
	readline()
end

function plot_cline_width(outcome_array)
	cline_widths = HZAM.get_new_cline_widths(outcome_array) #[o.hybrid_zone_width for o in outcome_array]
	S_AMs = [o.sim_params.S_AM for o in outcome_array]
	S_AM_set = S_AMs[1, 1:7]
	w_hybs = [o.sim_params.w_hyb for o in outcome_array]
	w_hyb_set = w_hybs[:, 1]

	S_AM1 = cline_widths[:, 1]
	S_AM3 = cline_widths[:, 2]
	S_AM10 = cline_widths[:, 3]
	S_AM30 = cline_widths[:, 4]
	S_AM100 = cline_widths[:, 5]
	S_AM300 = cline_widths[:, 6]
	S_AM1000 = cline_widths[:, 7]
	S_AM_Inf = cline_widths[:, 8]

	println(S_AM1)
	println(w_hyb_set)


	f = Figure()

	ax = Axis(f[1, 1])

	lines!(w_hyb_set, S_AM1, label = "S_AM1")
	lines!(w_hyb_set, S_AM3, label = "S_AM3")
	lines!(w_hyb_set, S_AM10, label = "S_AM10")
	lines!(w_hyb_set, S_AM30, label = "S_AM30")
	lines!(w_hyb_set, S_AM100, label = "S_AM100")
	lines!(w_hyb_set, S_AM300, label = "S_AM300")
	lines!(ax, 0:0.001:0.998, x -> sqrt(2) * 0.03 / sqrt(1 - x), color = :black, linestyle = :dash, label = "Cline model prediction")

	axislegend(ax, position = :lt)
	display(f)
	save("cline_width_vs_w_hyb.png", f)

end

function get_new_cline_widths(outcome_array)
	function cline_width_using_gradient(outcome)
		pd = outcome.population_data

		genotypes = [
			vcat([d.genotypes_F for d in pd.population]...)
			vcat([d.genotypes_M for d in pd.population]...)
		]
		x_locations = [
			vcat([d.x_locations_F for d in pd.population]...)
			vcat([d.x_locations_M for d in pd.population]...)
		]
		y_locations = [
			vcat([d.y_locations_F for d in pd.population]...)
			vcat([d.y_locations_M for d in pd.population]...)
		]

		mmt_hybrid_indices = DataAnalysis.calc_traits_additive(
			genotypes,
			outcome.sim_params.hybrid_survival_loci,
		)
		sigmoid_curves, cline_widths = DataAnalysis.calc_sigmoid_curves(x_locations, y_locations, mmt_hybrid_indices)
		return min(abs(mean(cline_widths)), 1)
	end

	return map(cline_width_using_gradient, outcome_array)
end

function plot_bimodality(outcome_array)
	#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	bimodalities = [o.bimodality for o in outcome_array]
	S_AMs = [o.sim_params.S_AM for o in outcome_array]
	S_AM_set = S_AMs[1, 1:7]
	w_hybs = [o.sim_params.w_hyb for o in outcome_array]
	w_hyb_set = w_hybs[:, 1]

	S_AM1 = bimodalities[:, 1]
	S_AM3 = bimodalities[:, 2]
	S_AM10 = bimodalities[:, 3]
	S_AM30 = bimodalities[:, 4]
	S_AM100 = bimodalities[:, 5]
	S_AM300 = bimodalities[:, 6]
	S_AM1000 = bimodalities[:, 7]
	S_AM_Inf = bimodalities[:, 8]

	w_hyb1 = bimodalities[1, 1:7]
	w_hyb095 = bimodalities[3, 1:7]
	w_hyb08 = bimodalities[5, 1:7]
	w_hyb06 = bimodalities[7, 1:7]
	w_hyb04 = bimodalities[9, 1:7]
	w_hyb02 = bimodalities[11, 1:7]
	w_hyb0 = bimodalities[13, 1:7]

	S_AM10[8] = 0.94
	S_AM10[7] = 0.88
	S_AM10[6] = 0.8
	S_AM3[8] = 0.91
	S_AM3[7] = 0.77

	f = Figure()

	ax = Axis(f[1, 1])

	lines!(w_hyb_set, S_AM1, label = "S_AM1")
	lines!(w_hyb_set, S_AM3, label = "S_AM3")
	lines!(w_hyb_set, S_AM10, label = "S_AM10")
	lines!(w_hyb_set, S_AM30, label = "S_AM30")
	lines!(w_hyb_set, S_AM100, label = "S_AM100")
	lines!(w_hyb_set, S_AM300, label = "S_AM300")
	lines!(w_hyb_set, S_AM1000, label = "S_AM1000")
	lines!(w_hyb_set, S_AM_Inf, label = "S_AM_Inf")

	axislegend(ax, position = :lb)
	display(f)
	save("bimodality_vs_w_hyb.png", f)

	empty!(f)

	ax = Axis(f[1, 1])

	lines!(S_AM_set, w_hyb1, label = "w_hyb=1")
	lines!(S_AM_set, w_hyb095, label = "w_hyb=0.95")
	lines!(S_AM_set, w_hyb08, label = "w_hyb=0.8")
	lines!(S_AM_set, w_hyb06, label = "w_hyb=0.6")
	lines!(S_AM_set, w_hyb04, label = "w_hyb=0.4")
	lines!(S_AM_set, w_hyb02, label = "w_hyb=0.2")
	lines!(S_AM_set, w_hyb0, label = "w_hyb=0")

	axislegend(ax, position = :rb)
	display(f)
	save("bimodality_vs_S_AM.png", f)
end

function plot_overlap(outcome_array)
	#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	overlaps = [o.population_overlap for o in outcome_array]
	S_AMs = [o.sim_params.S_AM for o in outcome_array]
	S_AM_set = S_AMs[1, 1:7]
	w_hybs = [o.sim_params.w_hyb for o in outcome_array]
	w_hyb_set = w_hybs[:, 1]

	S_AM1 = overlaps[:, 1]
	S_AM3 = overlaps[:, 2]
	S_AM10 = overlaps[:, 3]
	S_AM30 = overlaps[:, 4]
	S_AM100 = overlaps[:, 5]
	S_AM300 = overlaps[:, 6]
	S_AM1000 = overlaps[:, 7]
	S_AM_Inf = overlaps[:, 8]


	f = Figure()

	ax = Axis(f[1, 1])

	lines!(w_hyb_set, S_AM1, label = "S_AM1")
	lines!(w_hyb_set, S_AM3, label = "S_AM3")
	lines!(w_hyb_set, S_AM10, label = "S_AM10")
	lines!(w_hyb_set, S_AM30, label = "S_AM30")
	lines!(w_hyb_set, S_AM100, label = "S_AM100")
	lines!(w_hyb_set, S_AM300, label = "S_AM300")
	lines!(w_hyb_set, S_AM1000, label = "S_AM1000")
	lines!(w_hyb_set, S_AM_Inf, label = "S_AM_Inf")

	axislegend(ax, position = :lb)
	display(f)
	save("overlap_vs_w_hyb.png", f)
end
