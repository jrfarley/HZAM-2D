include("HZAM 2.0/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using GLMakie
using Statistics
using MannKendall
using Colors

colors = [
	RGB(177 / 255, 89 / 255, 40 / 255),
	RGB(246 / 255, 154 / 255, 153 / 255),
	RGB(228 / 255, 64 / 255, 32 / 255),
	RGB(166 / 255, 206 / 255, 227 / 255),
	RGB(31 / 255, 120 / 255, 180 / 255),
	RGB(106 / 255, 66 / 255, 154 / 255),
]

function categorize(outcome)#=
	if missing(outcome)
		return -1
	end

	widths = outcome.hybrid_zone_width
	widths = widths[length(widths) - 19 : end]
	overlaps = outcome.population_overlap
	overlaps = overlaps[length(overlaps) - 19 : end]

	mk_test_result = mk_original_test(widths)

	if pvalue(mk_test_result) < 0.05
		return -1
	end

	width = mean(widths)
	overlap = mean(overlaps)
	=#

	width = outcome.hybrid_zone_width
	overlap = outcome.population_overlap
	bimodality = HZAM.DataAnalysis.calc_bimodality(outcome.population_data, outcome.sim_params.male_mating_trait_loci)

	if width > 1
		return 1
	end

	if bimodality > 0.95
		return overlap < 0.3 ? 4 : 5
	elseif overlap < 0.1
		return width < 0.3 ? 2 : 3
	end

	return 6
end

function make_subplot(dir)
	outcome_array = HZAM.load_from_folder(dir)

	output_array = categorize.(outcome_array)

	println(output_array)
    output_array = output_array[end:-1:1, :]

    fig, ax, hm = heatmap(
        output_array',
        colormap = colors,
    )
    
    xs = collect(1:8)
    ys = collect(1:13)
    
    ax.xticks = (xs, ["1", "3", "10", "30", "100", "300", "1000", "Inf"])
    ax.yticks=(ys, string.(reverse([1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])))
    ax.xlabel = "S_AM"
    ax.ylabel = "w_hyb"
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
    ax.xlabelsize = 15
    ax.ylabelsize = 15
    ax.aspect = 1
    
    display(fig)
    readline()

end

make_subplot("HZAM_Sym_Julia_results_GitIgnore/HZAM_simulation_outcomes_Feb2024_GitIgnore/full_pleiotropy")
