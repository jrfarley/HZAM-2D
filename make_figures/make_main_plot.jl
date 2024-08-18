include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using CairoMakie # used for plotting
using Colors
using DataFrames
using GLM

ntz = []
zho = []

"Colours associated with the 6 outcome types."
colors = [
	RGB(246 / 255, 154 / 255, 153 / 255),
	RGB(228 / 255, 64 / 255, 32 / 255),
	RGB(166 / 255, 206 / 255, 227 / 255),
	RGB(31 / 255, 120 / 255, 180 / 255),
	RGB(177 / 255, 89 / 255, 40 / 255),
	#RGB(106 / 255, 66 / 255, 154 / 255),
]

"The directory where the simulation run is stored."
dir = "$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/Run2_2024July22/"

"Folders with all simulation outcome files from each set used for the main plot"
folders = [
	"full_pleiotropy",
	"low_reject_full_pleiotropy",
	"high_reject_full_pleiotropy",
	"no_pleiotropy",
	"low_reject_no_pleiotropy",
	"high_reject_no_pleiotropy",
	"separate_fmt",
	"separate_mmt",
	"separate_hst",
]

"Subplot titles"
names = [
	"Full Pleiotropy (FP)",
	"Low Reject Cost, FP",
	"High Reject Cost, FP",
	"No Pleiotropy (NP)",
	"Low Reject Cost, NP",
	"High Reject Cost, NP",
	"Magic Cue",
	"Magic Preference",
	"Matching",
]

"Colour squares for the legend"
group_color = [PolyElement(color = color, strokecolor = :transparent)
			   for color in colors]

"Outcome names for the legend"
labels = ["Bimodal hybrid zone", "Unimodal hybrid zone", "Narrow overlap zone", "Broad overlap zone", "Blended",]


"""
	categorize(outcome::Union{HZAM.DataAnalysis.OutputData, Missing})

Classify the simulation outcome into one of six outcome types based on the cline width, bimodality, and population overlap.
"""
function categorize(outcome::Union{HZAM.DataAnalysis.OutputData, Missing})
	
	# When the cline width is greater than 1, the simulations exit early. Thus missing simulations are classified as blended.
	if ismissing(outcome)
		return 5
	end

	widths = outcome.hybrid_zone_width
	last_widths = widths[length(widths)-19:end]
	overlaps = outcome.population_overlap
	last_overlaps = overlaps[length(overlaps)-19:end]
	bimodality = HZAM.DataAnalysis.calc_bimodality(outcome.population_data, outcome.sim_params.male_mating_trait_loci)

	width = sum(last_widths) / length(last_widths)
	overlap = sum(last_overlaps) / length(last_overlaps)

	generations = collect(1:20)

	df = DataFrame(X=generations, Y=last_widths)

	# Linear regression of cline width vs generation
	ols = lm(@formula(Y ~ X), df)

	# If the slope exceeds 0.03 units per 1000 generations (95% confidence interval) then the simulation is categorized as blended
	if confint(ols)[2,1] > 0.00075 && bimodality ≤ 0.95
		return 5
	end
	
	if bimodality > 0.95
		return overlap < 0.3 ? 3 : 4
	elseif bimodality < 0.5
		return 2
	end

	return 1
end

"""
	save_outcomes()

Load all of the simulations from `directories`, classify into outcome types and save the results.
"""
function save_outcomes()
	outcomes = []

	for i in 1:9
		outcome_array = HZAM.load_from_folder(string(dir, folders[i]))
		output_array = categorize.(outcome_array)
		push!(outcomes, output_array)
	end
	@save "MainPlot_GitIgnore/categorical_outcomes.JLD2" outcomes
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
function make_subplot(output_array::Array{<:Real}, name::String, ax::Makie.Axis)
	println(output_array)
	output_array = output_array[end:-1:1, :]

	heatmap!(
		ax,
		output_array',
		colormap = colors,
		colorrange = (1,5),
	)

	xs = collect(1:8)
	ys = collect(1:13)

	ax.xticks = (xs, ["1", "3", "10", "30", "100", "300", "1000", "Inf"])
	ax.yticks = (ys, string.(reverse([1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])))
	ax.xlabel = "S_AM"
	ax.ylabel = "w_hyb"
	ax.xticklabelsize = 18
	ax.yticklabelsize = 18
	ax.xlabelsize = 20
	ax.ylabelsize = 20
	ax.titlesize = 22
	ax.aspect = 1
	ax.title = name
	ax.xticklabelrotation = π / 2
end

"""
	make_main_plot()

Construct the main figure showing the outcomes of all the simulations.
"""
function make_main_plot()
	f = Figure(size = (1000, 1250))

	@load "MainPlot_GitIgnore/categorical_outcomes.JLD2" outcomes
	axes = Axis[]

	for i in 1:9
		ax = Axis(f[Int(ceil(i / 3)), (i-1)%3+1])
		make_subplot(outcomes[i], names[i], ax)
		push!(axes, ax)
	end
	
	leg = Legend(
		f[4, :],
		group_color,
		labels,
		orientation = :horizontal,
		tellheight = true,
		labelsize = 20,
		nbanks = 2,
		labelvalign = :center
	)
	trim!(f.layout)
	display(f)
	save("$(dirname(@__DIR__))/figures/main_plot3.png", f)
end

#save_outcomes()
#=
f = Figure()
hist(f[1, 1], ntz)
hist(f[1, 2], zho)
display(f)
readline()
=#
make_main_plot()
