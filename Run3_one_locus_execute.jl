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

	import .HZAM
	using Dates
	HZAM.set_results_folder("HZAM-J_2D_results/Run3_one_and_nine_loci_$(Dates.format(today(), "yyyymmdd"))")
end

# run main simulation sets
HZAM.run_HZAM_sets_complete_one_locus(set_numbers=[3])
HZAM.run_HZAM_sets_complete_nine_loci(set_numbers=[1,2,4,5,6,7,8,9])
