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
end

# run main simulation sets
HZAM.run_HZAM_sets_complete_nine_loci("Run3_replicate_1", set_numbers=2:9)
HZAM.run_HZAM_sets_complete_one_locus("Run3_replicate_1", set_numbers=3:9)

for i in 2:4
	HZAM.run_HZAM_sets_complete_one_locus("Run3_replicate_$i")
	HZAM.run_HZAM_sets_complete_nine_loci("Run3_replicate_$i")
end
