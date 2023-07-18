using CategoricalArrays
using Colors, ColorSchemes
import ColorSchemes.plasma
#using Plots.PlotMeasures  # needed for plot margin adjustment
using GLMakie
using Statistics: mean  # needed for "mean" function
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as CSV

function run_HZAM_set(set_name::String, intrinsic_R, ecolDiff;  # the semicolon makes the following optional keyword arguments 
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3, neutral_loci=4:6, total_functional_loci=1:3, num_neutral_loci=3,
    survival_fitness_method::String="epistasis", per_reject_cost=0, optimize=true)
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
    w_hyb_set = [1, 0.95, 0.9, 0.7, 0.5, 0.3, 0.1, 0]#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    # the set of assortative mating strengths (S_AM) that will be run
    S_AM_set = [1, 3, 10, 30, 100, 300, 1000]#[1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

    #ecolDiff_set = [1]

    if survival_fitness_method == "epistasis"
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        short_survFitnessMethod = "Het"
    end

    # set up array of strings to record outcomes
    outcome_array = Array{OutputData,3}(undef, length(w_hyb_set), length(S_AM_set), length(ecolDiff_set))

    outcome_w_hyb = Array{Float64,3}(undef, length(w_hyb_set), length(S_AM_set), length(ecolDiff_set))

    outcome_S_AM = Array{Integer,3}(undef, length(w_hyb_set), length(S_AM_set), length(ecolDiff_set))

    for k in eachindex(ecolDiff_set)  # loop through the replicate runs
        #ecolDiff = ecolDiff_set[k]

        # Loop through the different simulation sets
        Threads.@threads for i in eachindex(w_hyb_set)
            for j in eachindex(S_AM_set)
                w_hyb = w_hyb_set[i]
                S_AM = S_AM_set[j]
                println("ecolDiff = ", ecolDiff, "; w_hyb = ", w_hyb, "; S_AM = ", S_AM)

                run_name = string("HZAM_animation_run", set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM)

                # set up initial values for one simulation
                extinction = false
                K = K_total

                #K = (1+ecolDiff) * K_total

                # run one simulation by calling the function defined above:
                outcome = run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;
                    K_total=Int(trunc(K)), max_generations,
                    total_loci, female_mating_trait_loci, male_mating_trait_loci,
                    competition_trait_loci, hybrid_survival_loci,
                    survival_fitness_method, per_reject_cost, do_plot=false)

                if extinction  # whole simulation went extinct
                #= outcome = "extinction"
                 if save_each_sim
                     @save string("HZAM_Sym_Julia_results_GitIgnore/simulation_data.", run_name, ".jld2") outcome
                 end=#
                else  # no complete extinction

                    if save_each_sim
                        @save string("HZAM_Sym_Julia_results_GitIgnore/simulation_data.", run_name, ".jld2") String(string(outcome))
                    end
                end
                println(run_name, "  outcome was: ", outcome)
                outcome_array[i, j, k] = outcome
                outcome_w_hyb[i, j, k] = w_hyb
                outcome_S_AM[i, j, k] = S_AM
            end # of S_AM loop
        end # of w_hyb loop   
    end # of replicate loop

    if save_outcomes_JL
        filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set", set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff_set, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, ".jld2")
        @save filename outcome_array
    end

    if save_outcomes_csv
        for i in 1:size(outcome_array, 3)
            filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set", set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff_set[i], "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost)
            CSV.write(filename, Tables.table(outcome_array[:, :, i]), writeheader=false)
        end
    end
    return outcome_array, outcome_w_hyb, outcome_S_AM
end


#### functions for summarizing and plotting results

function plot_bimodality(outcomes, w_hybs, S_AMs, run_name)
    bimodalities = [map(o -> o.bimodality, outcomes)...]
    w_hybs = vcat(w_hybs...)
    S_AMs = vcat(S_AMs...)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string("Bimodality--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, w_hybs, S_AMs, color=bimodalities, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="bimodality", height=Relative(0.5))

    fig
end

function plot_gene_flow(outcomes, w_hybs, S_AMs, run_name)
    gene_flow = [map(o -> o.gene_flow, outcomes)...]
    w_hybs = vcat(w_hybs...)
    S_AMs = vcat(S_AMs...)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string("Gene Flow--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, w_hybs, S_AMs, color=gene_flow, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="gene flow", height=Relative(0.5))

    fig
end

function plot_width(outcomes, w_hybs, S_AMs, run_name)
    hybrid_zone_width = [map(o -> o.hybrid_zone_width, outcomes)...]
    w_hybs = vcat(w_hybs...)
    S_AMs = vcat(S_AMs...)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string("Hybrid Zone Width--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, w_hybs, S_AMs, color=hybrid_zone_width, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="cline width", height=Relative(0.5))

    fig
end


function plot_overlap(outcomes, w_hybs, S_AMs, run_name)
    overlap = [map(o -> o.overlap, outcomes)...]
    w_hybs = vcat(w_hybs...)
    S_AMs = vcat(S_AMs...)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string("Overlap--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, w_hybs, S_AMs, color=overlap, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="overlap", height=Relative(0.5))

    fig
end

function make_and_save_figs(ResultsFolder, run_name, outcomes, w_hybs, S_AMs)
    save(string(ResultsFolder, "/", run_name, "_bimodality.png"), plot_bimodality(outcomes, w_hybs, S_AMs, run_name))
    save(string(ResultsFolder, "/", run_name, "_gene_flow.png"), plot_gene_flow(outcomes, w_hybs, S_AMs, run_name))
    save(string(ResultsFolder, "/", run_name, "_zone_width.png"), plot_width(outcomes, w_hybs, S_AMs, run_name))
    save(string(ResultsFolder, "/", run_name, "_overlap.png"), plot_overlap(outcomes, w_hybs, S_AMs, run_name))
end