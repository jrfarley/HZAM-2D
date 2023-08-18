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

# the set of hybrid fitnesses (w_hyb) values that will be run
w_hyb_set = [0.95, 0.9, 0.7, 0.5]#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
# the set of assortative mating strengths (S_AM) that will be run
S_AM_set = [1, 10, 100, 1000]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

#=
parameter space for next time:

w_hyb_set = [1,0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5]

S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]

=#

intrinsic_R = 1.1
ecolDiff = 0
K_total = 20000
max_generations = 1000
sigma_disp = 0.05
total_loci = 10
female_mating_trait_loci = 1:2
male_mating_trait_loci = 3:4
competition_trait_loci = 1:2
hybrid_survival_loci = 5:6
num_neutral_loci = 4
survival_fitness_method = "epistasis"
per_reject_cost = 0.05


# set up array of strings to record outcomes
outcome_array = Array{Any,2}(undef, length(w_hyb_set), length(S_AM_set))

sim_params = Array{SimulationParameters,2}(undef, length(w_hyb_set), length(S_AM_set))


dir = mkpath(string("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/", "trial1"))

# Loop through the different simulation sets
Threads.@threads for i in eachindex(w_hyb_set)
    for j in eachindex(S_AM_set)
        w_hyb = w_hyb_set[i]
        S_AM = S_AM_set[j]
        println("ecolDiff = ", ecolDiff, "; w_hyb = ", w_hyb, "; S_AM = ", S_AM)

        # set up initial values for one simulation
        extinction = false
        K = K_total

        #K = (1+ecolDiff) * K_total

        # run one simulation by calling the function defined above:
        outcome = HZAM.run_one_HZAM_sim(w_hyb, S_AM, 0, 1.1; # these values are 
            # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
            K_total=K, max_generations=1000,
            sigma_disp=0.05, sigma_comp=0.01, do_plot=false, plot_int=20,
            total_loci=10,
            female_mating_trait_loci=1:2,
            male_mating_trait_loci=3:4,
            competition_trait_loci=1:2,
            hybrid_survival_loci=5:6,
            per_reject_cost=0.05,
            gene_plot=false,
            save_plot=false
        )

        parameters = SimulationParameters(
            intrinsic_R,
            ecolDiff,
            w_hyb,
            S_AM,
            K_total,
            max_generations,
            sigma_disp,
            total_loci,
            female_mating_trait_loci,
            male_mating_trait_loci,
            competition_trait_loci,
            hybrid_survival_loci,
            num_neutral_loci,
            survival_fitness_method,
            per_reject_cost
        )
        println(string("S_AM: ", S_AM, ", w_hyb: ", w_hyb))
        println("Outcome was: ", length(outcome))
        outcome_array[i, j] = outcome
        sim_params[i, j] = parameters
    end # of S_AM loop
end # of w_hyb loop

filename = string(dir, "/", "search_cost", ".jld2")
@save filename sim_params outcome_array
