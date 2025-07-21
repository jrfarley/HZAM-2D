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

# run the three loci sets 4 more times
for i in 1:4
	HZAM.run_HZAM_sets_complete_three_loci("Run3_replicate_$i")
end
