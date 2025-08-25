# July 2025 simulation run

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
end

# run main simulation sets
#=
for i in 1:3
	HZAM.run_HZAM_sets_complete_nine_loci("Run4_fmt_cline_$i", set_numbers = [4, 10, 11], cline_width_loci="fmt")
end
=#
first_to_three(
		["HZAM-J_2D_results/Run4_fmt_cline_$(i)_nine_loci" for i in 1:3],
		"HZAM-J_2D_results/Run4_fmt_cline_4_nine_loci",
		[5,6],
		"fmt",
	)

for j in 5:8
first_to_three(
		["HZAM-J_2D_results/Run4_fmt_cline_$(i)_nine_loci" for i in 1:(j-1)],
		"HZAM-J_2D_results/Run4_fmt_cline_$(j)_nine_loci",
		[4,5,6],
		"fmt",
	)
end
