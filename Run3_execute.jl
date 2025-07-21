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
end

# run main simulation sets
HZAM.run_HZAM_sets_complete_nine_loci(set_numbers=[3])

for i in 1:2
	HZAM.run_HZAM_sets_complete_one_locus()
end
