using BenchmarkTools
using CategoricalArrays
using Colors, ColorSchemes
import ColorSchemes.plasma
using Plots

include("simulation.jl")

function run_HZAM_set(set_name::String, ecolDiff, intrinsic_R, replications;  # the semicolon makes the following optional keyword arguments 
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3, neutral_loci=4:6, total_functional_loci=1:3, num_neutral_loci=3,
    survival_fitness_method::String="epistasis", per_reject_cost=0,
    starting_pop_ratio=1.0, optimize=true)
    # ecolDiff should be from 0 to 1 (parameter "E" in the paper)
    # intrinsic_R is called "R" in the paper
    # replications should be somehting like "1:10" or just "1" for 1 replicate, or something like "2:5" to add replicates after 1 is done
    # K_total (default 1000) is the total carrying capacity, divided equally between K_alpha and K_beta
    # max_generations (default 1000) is the number of generations each simulation will run 
    # total_loci (default 6) is the total number of loci, regardless of their role (with indices 1:total_loci, referred to below)
    # Need to specify indices (columns) of four types of functional loci (can be the same). At least one should begin with index 1:
    # female_mating_trait_loci (default 1:3) is indices of the loci that determine the female mating trait
    # male_mating_trait_loci (default 1:3) is indices of the loci that determine the male mating trait
    # competition_trait_loci (default 1:3) is indices of the loci that determine the ecological trait (used in fitness related to resource use)
    # hybrid_survival_loci (default 1:3) is indices of the loci that determine survival probability of offspring to adulthood (can be viewed as incompatibilities and/or fitness valley based on ecology)
    # Specify indices (columns) of neutral loci (which have no effect on anything, just along for the ride):
    # neutral_loci (default 4:6) is indices of neutral loci (used for neutral measure of hybrid index; not used in the HZAM-sym paper)
    # per_reject_cost (default 0, can take values of 0 to 1) is fitness loss of female per male rejected (due to search time, etc.)
    # starting_pop_ratio (default 1) is the starting ratio of pop A to pop B

    save_outcomes_JL = true
    save_outcomes_csv = true  # whether to save the whole outcome array as csv files (with each rep as separate file)
    save_each_sim = false  # whether to save detailed data for each simulation

    # the set of hybrid fitnesses (w_hyb) values that will be run
    w_hyb_set = [1.0]#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    # the set of assortative mating strengths (S_AM) that will be run
    S_AM_set = [300]#[1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

    if survival_fitness_method == "epistasis"
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        short_survFitnessMethod = "Het"
    end

    # set up array of strings to record outcomes
    outcome_array = Array{String,3}(undef, length(w_hyb_set), length(S_AM_set), length(replications))

    for k in 1:length(replications)  # loop through the replicate runs
        replicate_ID = replications[k]
        run_set_name = string(set_name, "_rep", replicate_ID)

        # Loop through the different simulation sets
        Threads.@threads for i in 1:length(w_hyb_set)
            for j in 1:length(S_AM_set)
                w_hyb = w_hyb_set[i]
                S_AM = S_AM_set[j]
                println("w_hyb = ", w_hyb, "; S_AM = ", S_AM)

                run_name = string("HZAM_animation_run", run_set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM)

                # set up initial values for one simulation
                extinction = false
                outcome = []
                final_distribution = []

                # run one simulation by calling the function defined above:
                extinction, functional_HI_all_inds, HI_NL_all_inds = run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;
                    K_total, max_generations,
                    total_loci, female_mating_trait_loci, male_mating_trait_loci,
                    competition_trait_loci, hybrid_survival_loci, neutral_loci,
                    survival_fitness_method, per_reject_cost,
                    starting_pop_ratio, do_plot=false, optimize)


                if extinction  # whole simulation went extinct
                    outcome = "extinction"
                    if save_each_sim
                        @save string("HZAM_Sym_Julia_results_GitIgnore/simulation_data.", run_name, ".jld2") outcome
                    end
                else  # no complete extinction
                    species0_proportion = sum(functional_HI_all_inds .< 0.1) / length(functional_HI_all_inds)
                    species1_proportion = sum(functional_HI_all_inds .> 0.9) / length(functional_HI_all_inds)
                    species0_proportion_NL = sum(HI_NL_all_inds .== 0) / length(HI_NL_all_inds)
                    species1_proportion_NL = sum(HI_NL_all_inds .== 1) / length(HI_NL_all_inds)

                    if species0_proportion >= 0.85 || species1_proportion >= 0.85
                        outcome = "one_species"
                    elseif (species0_proportion + species1_proportion >= 0.85) && (species0_proportion >= 0.15) && (species1_proportion >= 0.15)
                        outcome = "two_species"
                    else
                        outcome = "blended"
                    end

                    if save_each_sim
                        @save string("HZAM_Sym_Julia_results_GitIgnore/simulation_data.", run_name, ".jld2") outcome functional_HI_all_inds species0_proportion species1_proportion HI_NL_all_inds species0_proportion_NL species1_proportion_NL
                    end
                end
                println(run_name, "  outcome was: ", outcome)
                outcome_array[i, j, k] = outcome
            end # of S_AM loop
        end # of w_hyb loop   
    end # of replicate loop

    if save_outcomes_JL
        filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set", set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, ".jld2")
        @save filename outcome_array
    end

    if save_outcomes_csv
        for i in 1:size(outcome_array, 3)
            filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set", set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, "_rep", replications[i])
            CSV.write(filename, Tables.table(outcome_array[:, :, i]), writeheader=false)
        end
    end
    return outcome_array
end

#### functions for summarizing and plotting results

# convert outcome array (an array of strings) to categorical array
function convert_to_cat_array(outcome_array)
    cat_outcome_array = compress(CategoricalArray(outcome_array))
    levels!(cat_outcome_array, ["extinction", "blended", "one_species", "two_species"])
    return cat_outcome_array
end 

# function for plotting grid of pie charts showing distribution of four outcomes, using categorical array as input
function plot_all_outcomes(cat_outcome_array)
    num_outcome_types = length(levels(cat_outcome_array))
    outcome_counts = Array{Int16, 2}(undef, num_outcome_types, (size(cat_outcome_array, 1)*size(cat_outcome_array, 2) )) 
    for i in 1:size(cat_outcome_array, 1) 
        for j in 1:size(cat_outcome_array, 2) 
            for outcome_num in 1:num_outcome_types
                outcome_counts[outcome_num, j + (i-1)*size(cat_outcome_array, 2)] = sum(cat_outcome_array[i,j,:] .== levels(cat_outcome_array)[outcome_num])
            end
        end
    end
    colors_of_outcomes = [RGB(0,0,0), plasma[0.525], plasma[0.2], plasma[0.9]] # colors for 4 outcome categories
    pie(outcome_counts, layout = grid(size(cat_outcome_array, 1), size(cat_outcome_array, 2)), legend = false, palette = colors_of_outcomes#=, margin = -2.0mm=#)
    plot!(size=(800,1300))
end

function get_most_common_outcomes(cat_outcome_array)
    levels_of_outcomes = levels(cat_outcome_array)
    num_outcome_types = length(levels_of_outcomes)
    most_common_outcomes = CategoricalArray{String, 2}(undef, size(cat_outcome_array, 1), size(cat_outcome_array, 2))
    levels!(most_common_outcomes, ["extinction", "blended", "one_species", "two_species"])
    for i in 1:size(cat_outcome_array, 1) 
        for j in 1:size(cat_outcome_array, 2)
            outcome_counts = Vector{Int}(undef, num_outcome_types) 
            for outcome_num in 1:num_outcome_types
                outcome_counts[outcome_num] = sum(cat_outcome_array[i,j,:] .== levels(cat_outcome_array)[outcome_num])
            end
            outcomes_with_max_count = findall(outcome_counts .== maximum(outcome_counts))
            if length(outcomes_with_max_count) == 1
                most_common_outcomes[i,j] = levels_of_outcomes[outcomes_with_max_count[1]]
            elseif length(outcomes_with_max_count) >= 2
                # if tie in outcome count, choose one randomly
                most_common_outcomes[i,j] = levels_of_outcomes[sample(outcomes_with_max_count)]
            end
        end
    end
    return most_common_outcomes
end

function plot_common_outcomes(common_outcome_array)
    # make heat map of outcomes
    colors_of_outcomes = [RGB(0,0,0), plasma[0.525], plasma[0.2], plasma[0.9]] # colors for 4 outcome categories
    one_outcome_array = reverse(common_outcome_array, dims=1)
    x_midpoints = log10.([1, 3, 10, 30, 100, 300, 1000, 5000])  # the S_AM values, with Inf convert to 3000 for graphing 
    w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    y_midpoints = reverse(w_hyb_set)
    min_color = minimum(levelcode.(one_outcome_array)) # this and next line needed to choose the proper colors for the figure
    max_color = maximum(levelcode.(one_outcome_array))
    heatmap(x_midpoints, y_midpoints, one_outcome_array, c = colors_of_outcomes[min_color:max_color], yflip = false, tick_direction = :out, colorbar = false, size = (440,310), framestyle = :box)
    xaxis!("Strength of conspecific mate preference")
    xticklabels = ["1", "3", "10", "30", "100", "300", "1000", "complete"]
    plot!(xticks=(x_midpoints, xticklabels))
    yaxis!("Hybrid fitness")
    yticklabels = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "", "", "1.0"]
    plot!(yticks=(y_midpoints, yticklabels), tick_direction = :out)
    # add white lines
    plot!([xlims()[1], xlims()[2]], [0.05, 0.05], linecolor = :white, widen = false, legend = false, linewidth=3)
    x_for_line = mean(x_midpoints[[length(x_midpoints)-1 length(x_midpoints)]])
    plot!([x_for_line, x_for_line], [ylims()[1], ylims()[2]], linecolor = :white, widen = false, legend = false, linewidth=3)
end

function make_and_save_figs(ResultsFolder, RunName, RunOutcomes)
    cat_RunOutcomes = convert_to_cat_array(RunOutcomes)
    display(plot_all_outcomes(cat_RunOutcomes))
    savefig(string(ResultsFolder,"/",RunName,"_AllOutcomes.png"))
    #=savefig(string(ResultsFolder,"/",RunName,"_AllOutcomes.pdf"))
    most_common_outcomes = get_most_common_outcomes(cat_RunOutcomes)
    display(plot_common_outcomes(most_common_outcomes))
    savefig(string(ResultsFolder,"/",RunName,"_MostCommonOutcomes.png"))
    savefig(string(ResultsFolder,"/",RunName,"_MostCommonOutcomes.pdf"))=#
end
