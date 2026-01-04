using JLD2
using Distributed
include("../HZAM/src/HZAM.jl")
include("make_main_plot.jl")
import .HZAM

"Colours associated with the 6 outcome types."
colors = [
	RGB(246 / 255, 154 / 255, 153 / 255),
	RGB(228 / 255, 64 / 255, 32 / 255),
	RGB(166 / 255, 206 / 255, 227 / 255),
	RGB(31 / 255, 120 / 255, 180 / 255),
	RGB(177 / 255, 89 / 255, 40 / 255),
	RGB(0 / 255, 0 / 255, 0 / 255)#,
	#RGB(255 / 255, 255 / 255, 255 / 255),
]

LABELSIZE = 20



S_AM_set = [3, 10, 30, 100, 300, 1000]
sc_set = [0.001, 0.003, 0.01, 0.02, 0.05, 0.1]

function load_four_outcomes(index, path)
	j = 1
	outcomes = Array{Integer}(undef, 4)
	parameters = []

	for i in 1:8
		filename = string(path, "_", i, ".JLD2")
		try
			@load filename outcome_array parameters_array
			parameters = parameters_array[index]
			outcome = outcome_array[index]
			if outcome != 7
				outcomes[j] = outcome
				j+=1
			end
		catch
		end
		i += 1
		if j>4
			break
		end
	end

	if undef in outcomes
		println("ERROR")
		println(parameters)
	end

	return outcomes
end

function load_set(path)
	outcome_array = Array{Missing, 2}(
		undef, length(sc_set), length(S_AM_set),
	)

	load_four_outcomes.(CartesianIndices(outcome_array), Ref(path))
end

function subplot(path, gl)
	grid_data = load_set(path)
	println(grid_data)

	for i in 1:length(sc_set), j in 1:length(S_AM_set)
		ax = Axis(gl[i,j], aspect = 1)

		data = grid_data[i,j]

		# Count categories
		counts = [count(==(k), data) for k in 1:6]

		# Normalize for pie
		fracs = counts ./ sum(counts)
		println(fracs)

		pie!(ax, fracs, color = colors)
		hidedecorations!(ax)  # hides ticks, grid and lables
		hidespines!(ax)  # hide the frame
	end

	for i in 1:length(sc_set)
		Label(gl[i,0], "$(sc_set[i]*100)", tellwidth=false)
	end
	for i in 1:length(S_AM_set)
		Label(gl[7, i], "$(S_AM_set[i])", tellwidth=false)
	end

	Label(gl[1:6,-1], "Search cost (%)", rotation=pi/2, tellwidth=false)
	colsize!(gl, -1, 15)
	Label(gl[8,1:6], "Strength of conspecific mate preference", tellwidth=false)
	#Box(gl[1,0], color = :transparent, strokecolor = :black, strokewidth = 2)
end

function make_plot(path)
	paths = [
		"HZAM-J_2D_results_categorized/supplement/multivariate_replicate",
		"HZAM-J_2D_results_categorized/supplement/multivariate_w_hyb_095",
		"HZAM-J_2D_results_categorized/supplement/one_mmt_one_fmt_replicate",
		"HZAM-J_2D_results_categorized/supplement/one_mmt_one_fmt_w_hyb_095",
		"HZAM-J_2D_results_categorized/supplement/one_mmt_three_fmt_replicate",
		"HZAM-J_2D_results_categorized/supplement/one_mmt_three_fmt_w_hyb_095"
	]

	fig = Figure(resolution = (900,1100))
	g = fig[1,1] = GridLayout()

	for i in 1:6
		gl = g[2*((i+1)÷2)-1, 2+(-1)^(i%2)] = GridLayout()
		subplot(paths[i], gl)
		#Box(g[(i-1)÷2 + 1, (i-1) % 2 + 1], color = :transparent, strokecolor = :black) 
	end

	#Box(fig[1,1], color = :transparent, strokecolor = :black) 

	#Label(g[:,2], "")
	Label(g[0,1], "100% Hybrid fitness", tellwidth=false, fontsize=LABELSIZE)
	Label(g[0,3], "95% Hybrid fitness", tellwidth=false, fontsize = LABELSIZE)
	colsize!(g, 2, 10)

	Label(g[1,0], "Multivariate preference", rotation = π/2, tellwidth=false, fontsize=LABELSIZE)
	Label(g[3,0], "Single locus preference", rotation = π/2, tellwidth=false, fontsize=LABELSIZE)
	Label(g[5,0], "Multilocus preference", rotation = π/2, tellwidth=false, fontsize=LABELSIZE)
	colsize!(g, 0, LABELSIZE)
	rowsize!(g, 2, 10)
	rowsize!(g, 4, 10)
	
	#rowgap!(g, 50)
	#colgap!(g, 50)
	group_color = [PolyElement(color = color, strokecolor = :transparent)
			   for color in colors]
	labels = ["Bimodal hybrid zone", "Unimodal hybrid zone", "Narrow overlap zone", "Broad overlap zone", "Blended", "One species", "Uncategorizable"]


	leg = Legend(
		fig[2, :],
		group_color[1:6],
		labels[1:6],
		orientation = :horizontal,
		tellheight = true,
		labelsize = 20,
		nbanks = 2,
		labelvalign = :center,
	)


	display(fig)
	save("$(dirname(@__DIR__))/figures/test_supplement.png", fig)
	#readline()
end

make_plot("HZAM-J_2D_results_categorized/supplement/multivariate_w_hyb_095")