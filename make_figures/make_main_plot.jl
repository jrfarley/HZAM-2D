include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using CairoMakie # used for plotting
using Colors
using DataFrames
using GLM
using HypothesisTests
using MannKendall
using StatsBase
using Distributed

CairoMakie.activate!()

ntz = []
zho = []

"Colours associated with the 6 outcome types."
colors = [
	RGB(246 / 255, 154 / 255, 153 / 255),
	RGB(228 / 255, 64 / 255, 32 / 255),
	RGB(166 / 255, 206 / 255, 227 / 255),
	RGB(31 / 255, 120 / 255, 180 / 255),
	RGB(177 / 255, 89 / 255, 40 / 255),
	RGB(0 / 255, 0 / 255, 0 / 255),
	RGB(255 / 255, 255 / 255, 255 / 255),
]

"The directory where the simulation run is stored."
dir = "$(dirname(@__DIR__))/HZAM-J_2D_results/Run3_one_locus_20250715/"

"Folders with all simulation outcome files from each set used for the main plot"
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

"Subplot titles"
names = [
	"Full Pleiotropy",
	"Matching",
	"Magic Preference",
	"Magic Cue (MC)",
	"Low Reject Cost, MC",
	"High Reject Cost, MC",
	"No Pleiotropy (NP)",
	"Low Reject Cost, NP",
	"High Reject Cost, NP",
]

"Colour squares for the legend"
group_color = [PolyElement(color = color, strokecolor = :transparent)
			   for color in colors]

"Outcome names for the legend"
labels = ["Bimodal hybrid zone", "Unimodal hybrid zone", "Narrow overlap zone", "Broad overlap zone", "Blended", "One species", "Uncategorizable"]


"""
	categorize(outcome::Union{HZAM.DataAnalysis.OutputData, Missing})

Classify the simulation outcome into one of six outcome types based on the cline width, bimodality, and population overlap.
"""

function categorize(outcome::Union{HZAM.DataAnalysis.OutputData, Missing}; categorization_loci = "mmt")

	# When the cline width is greater than 1, the simulations exit early. Thus missing simulations are classified as blended.
	if ismissing(outcome)
		return 5
	end

	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}} = [1]

	if categorization_loci == "fmt"
		loci = outcome.sim_params.female_mating_trait_loci
	else
		loci = outcome.sim_params.male_mating_trait_loci
	end

	if is_one_species_extinct(outcome; categorization_loci)
		return 6
	end

	if is_blended(outcome)
		return 5
	end

	overlap = HZAM.Population.calc_species_overlap(
		outcome.population_data.population,
		0.03,
		0.01,
		loci,
	)[1]

	bimodality = HZAM.DataAnalysis.calc_bimodality(outcome.population_data, loci)

	if bimodality > 0.95
		return overlap < 0.3 ? 3 : 4
	elseif bimodality < 0.5
		return 2
	end

	return 1
end

function is_one_species_extinct(outcome::HZAM.DataAnalysis.OutputData; categorization_loci = "mmt")

	phenotype_counts = outcome.population_tracking_data[2]

	if categorization_loci=="fmt"
		phenotype_counts = outcome.population_tracking_data[1]
	end

	final_phenotype_counts = phenotype_counts[end]

	num_mmt_loci = length(outcome.sim_params.male_mating_trait_loci)

	species_A_proportion = sum(final_phenotype_counts[1:(Integer(floor(0.2*num_mmt_loci))+1)])

	species_B_proportion = sum(final_phenotype_counts[(Integer(ceil(1.8*num_mmt_loci))+1):length(final_phenotype_counts)])

	return (species_A_proportion == 0 && species_B_proportion ≥ 0.5) || (species_B_proportion == 0 && species_A_proportion>0.5)
end

function is_cline_width_increasing(outcome::HZAM.DataAnalysis.OutputData)
	cline_widths = outcome.hybrid_zone_width[11:30]


	test = mk_original_test(cline_widths)


	return test.trend=="increasing"
end

function is_blended(outcome::HZAM.DataAnalysis.OutputData)

	if outcome.bimodality > 0.95
		return false
	end

	max_cline_width = maximum(outcome.hybrid_zone_width)

	if max_cline_width > 1
		return true
	end

	return is_cline_width_increasing(outcome)
end


"""
	save_outcomes()

Load all of the simulations from `directories`, classify into outcome types and save the results.
"""
function save_outcomes(; output_dir = dir, output_name = "categorical_outcomes.JLD2", set_numbers = 1:9)
	outcomes = []

	for i in set_numbers
		outcome_array = HZAM.load_from_folder(string(output_dir, folders[i]))
		output_array = categorize.(outcome_array)
		push!(outcomes, output_array)
	end
	@save output_name outcomes
end

"""
	make_subplot(output_array::Array{<:Real}, name::String, ax::Makie.Axis)

Make a subplot for the main figure showing the outcome types for different combinations of 
reduced hybrid fitness and assortative mating.

# Arguments
- `output_array::Array{<:Real}`: the outcome type for every simulation for the subplot .
- `name::String` : the subplot title.
- `ax::Makie.Axis` : the axis where the subplot will be placed.
"""
function make_subplot(output_array::Array{<:Real}, ax::Makie.Axis)
	println(output_array)
	output_array = output_array[end:-1:1, :]

	heatmap!(
		ax,
		output_array',
		colormap = colors,
		colorrange = (1, 7),
	)

	xs = collect(1:8)
	ys = collect(1:13)

	ax.xticks = (xs, ["1", "3", "10", "30", "100", "300", "1000", "Inf"])
	ax.yticks = (ys, string.(reverse([1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])))


end

"""
	make_main_plot()

Construct the main figure showing the outcomes of all the simulations.
"""
function make_main_plot()
	@load "categorical_outcomes.JLD2" outcomes

	f = Figure(size = (1000, 1150), pt_per_unit = 2)

	g = f[1, 1] = GridLayout()
	axes = Axis[]


	for i in 1:9
		ax = Axis(
			g[Int(ceil(i / 3)), (i-1)%3+1],
			xlabel = "S_AM",
			ylabel = "w_hyb",
			xticklabelrotation = π / 2,
			title = names[i],
			aspect = 1,
			xticklabelsize = 18,
			yticklabelsize = 18,
			xlabelsize = 20,
			ylabelsize = 20,
			titlesize = 22,
		)
		make_subplot(outcomes[i], ax)
		push!(axes, ax)
	end

	leg = Legend(
		f[2, :],
		group_color,
		labels,
		orientation = :horizontal,
		tellheight = true,
		labelsize = 20,
		nbanks = 2,
		labelvalign = :center,
	)
	resize_to_layout!(f)
	save("$(dirname(@__DIR__))/figures/main_plot_test_one_locus.png", f)
end

function test_categorization(dir)
	outcome_array = HZAM.load_from_folder(dir)
	output_array = categorize.(outcome_array; categorization_loci = "fmt")
	println(output_array)
end





function make_plot_first_to_three(outcome_folders, set_numbers, output_name, categorization_loci)
	first_to_three_outcome_arrays = Vector{<:Any}(undef, 11)

	function first_to_three(categorized_outcomes)
		for i in 1:6
			occurences = count(x->x==i, categorized_outcomes)
			if occurences ≥ 3
				return i
			end
		end
		return 7
	end

	for set in set_numbers
		categorized_outcome_arrays = []
		for i in 1:10
			outcome_array = HZAM.load_from_folder("$(outcome_folders[i])/$(folders[set])")
			push!(categorized_outcome_arrays, categorize.(outcome_array; categorization_loci))

			first_to_three_outcome_array = first_to_three.([getindex.(categorized_outcome_arrays, j) for j in CartesianIndices(outcome_array)])

			println(categorized_outcome_arrays[end])
			println(first_to_three_outcome_array)
			first_to_three_outcome_arrays[set] = first_to_three_outcome_array
			if 7 in first_to_three_outcome_array
				continue
			else
				break
			end
		end
	end

	f = Figure(size = (1000, 1150), pt_per_unit = 2)

	g = f[1, 1] = GridLayout()
	axes = Vector{<:Any}(undef, 9)


	for i in set_numbers
		ax = Axis(
			g[Int(ceil(i / 3)), (i-1)%3+1],
			xlabel = "S_AM",
			ylabel = "w_hyb",
			xticklabelrotation = π / 2,
			title = names[i],
			aspect = 1,
			xticklabelsize = 18,
			yticklabelsize = 18,
			xlabelsize = 20,
			ylabelsize = 20,
			titlesize = 22,
		)
		make_subplot(first_to_three_outcome_arrays[i], ax)
		push!(axes, ax)
	end
	leg = Legend(
		f[2, :],
		group_color,
		labels,
		orientation = :horizontal,
		tellheight = true,
		labelsize = 20,
		nbanks = 2,
		labelvalign = :center,
	)
	resize_to_layout!(f)

	display(f)
	save("$(dirname(@__DIR__))/figures/$output_name.png", f)
	readline()
end

function save_outcomes_first_to_three(folders)
	for i in eachindex(folders)
		save_outcomes(output_dir = folders[i], output_name = "categorical_outcomes$(i).JLD2", set_numbers = 4:6)
	end
end


function first_to_three(input_folders, output_folder, set_numbers, cline_width_loci)

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

		outcome_arrays = [HZAM.load_from_folder("$(input_folder)/$(folders[set_number])") for input_folder in input_folders]

		categorized_outcome_arrays = [categorize.(outcome_array) for outcome_array in outcome_arrays]


		first_to_three_outcome_array = first_to_three.([getindex.(categorized_outcome_arrays, i) for i in CartesianIndices(outcome_arrays[1])])

		indices = findall(x->x==7, first_to_three_outcome_array)
		results_folder = mkpath("$(output_folder)/$(folders[set_number])")
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
end
#save_outcomes()
#make_main_plot()
#=
make_plot_first_to_three(
	["HZAM-J_2D_results/Run4_fmt_cline_$(i)_nine_loci" for i in 1:8],
	4:6,
	"nine_loci_fmt",
	"fmt",
)
	=#

