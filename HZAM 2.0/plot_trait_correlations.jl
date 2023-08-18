include("HZAM/src/HZAM.jl")
import .HZAM
using JLD2

struct SimulationParameters
    intrinsic_R
    ecolDiff
    w_hyb
    S_AM
    K_total
    max_generations
    sigma_disp
    total_loci
    female_mating_trait_loci
    male_mating_trait_loci
    competition_trait_loci
    hybrid_survival_loci
    num_neutral_loci
    survival_fitness_method
    per_reject_cost
end


dir = mkpath(string("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/", "trial1"))
#=
filename = string(dir, "/", "no_magic", ".jld2")

@load filename sim_params outcome_array

println(HZAM.find_extinct_alleles(outcome_array[1,2]))
println(sim_params[1,2])
readline()

outcome = HZAM.run_one_HZAM_sim(0.95, 10, 0, 1.1; # these values are 
# hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
K_total=20000, max_generations=1000,
sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=20,
total_loci=10,
female_mating_trait_loci=1:2,
male_mating_trait_loci=3:4,
competition_trait_loci=1:2,
hybrid_survival_loci=5:6,
per_reject_cost=0,
gene_plot=true,
save_plot=false
)
outcome_array[1,2] = outcome
@save filename sim_params outcome_array
=#
HZAM.summarize_gene_correlations(dir)
