include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using CairoMakie # used for plotting
using Colors
using DataFrames
using GLM
#using HypothesisTests
using MannKendall
#using StatsBase
using Distributed

CairoMakie.activate!()
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
dir = "$(dirname(@__DIR__))/HZAM-J_2D_results/"

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
	"Low Search Cost, MC",
	"High Search Cost, MC",
	"No Pleiotropy (NP)",
	"Low Search Cost, NP",
	"High Search Cost, NP",
]

"Colour squares for the legend"
group_color = [PolyElement(color = color, strokecolor = :transparent)
			   for color in colors]

"Outcome names for the legend"
labels = ["Bimodal hybrid zone", "Unimodal hybrid zone", "Narrow overlap zone", "Broad overlap zone", "Blended", "One species", "Uncategorizable"]

w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]


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
function make_main_plot(dir, set_numbers)
	outcomes = Vector{<:Any}(undef, 9)

	for i in set_numbers
		@load string(dir, "/$(folders[i]).JLD2") outcome_array
		outcomes[i] = outcome_array
	end

	f = Figure(size = (1000, 1150), pt_per_unit = 2)

	g = f[1, 1] = GridLayout()
	axes = Axis[]


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
		make_subplot(outcomes[i], ax)
		push!(axes, ax)
	end

	leg = Legend(
		f[2, :],
		group_color[1:6],
		labels[1:6],
		orientation = :horizontal,
		tellheight = true,
		labelsize = 20,
		nbanks = 2,
		labelvalign = :center,
	)
	resize_to_layout!(f)
	save("$(dirname(@__DIR__))/figures/main_plot_test.png", f)
end

function make_plot_first_to_three(outcome_folders, set_numbers, output_name; single_row = false)
	first_to_three_outcome_arrays = Vector{<:Any}(undef, 11)

	function first_to_three(categorized_outcomes)
		for i in 1:6
			occurences = count(x->!ismissing(x) && x==i, categorized_outcomes)
			if occurences ≥ 3
				return i
			end
		end
		return 7
	end

	for set in set_numbers
		categorized_outcome_arrays = []
		for folder in outcome_folders
			try
				@load "$(folder)/$(folders[set]).JLD2" outcome_array

				push!(categorized_outcome_arrays, outcome_array)
			catch
			end
		end

		first_to_three_outcome_array = first_to_three.([getindex.(categorized_outcome_arrays, j) for j in CartesianIndices(categorized_outcome_arrays[end])])

		println(first_to_three_outcome_array)
		first_to_three_outcome_arrays[set] = first_to_three_outcome_array
	end

	if !single_row
		f = Figure(size = (1000, 1150), pt_per_unit = 2)
	else
		f=Figure(size = (1000, 500), pt_per_unit = 2)
	end

	g = f[1, 1] = GridLayout()
	axes = Vector{<:Any}(undef, 9)


	for i in set_numbers
		row = single_row ? 1 : Int(ceil(i / 3))
		ax = Axis(
			g[row, (i-1)%3+1],
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
		group_color[1:6],
		labels[1:6],
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

#make_plot_first_to_three(["HZAM-J_2D_results_categorized/mmt_three_loci/Run3_three_loci_$i" for i in 1:11], 1:9, "three_loci_mmt")
#make_plot_first_to_three(["HZAM-J_2D_results_categorized/mmt_nine_loci/Run3_nine_loci_$i" for i in 1:10], 1:9, "nine_loci_mmt")
#make_plot_first_to_three(["HZAM-J_2D_results_categorized/mmt_one_locus/Run3_one_locus_$i" for i in 1:9], 1:9, "one_locus_mmt")
make_plot_first_to_three(["HZAM-J_2D_results_categorized/fmt_three_loci/fmt_clines_three_loci_$i" for i in 1:11], 4:6, "three_loci_fmt"; single_row = true)

#make_main_plot()


