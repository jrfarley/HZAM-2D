using CategoricalArrays
using Colors, ColorSchemes
import ColorSchemes.plasma
#using Plots.PlotMeasures  # needed for plot margin adjustment
using GLMakie
using Statistics: mean  # needed for "mean" function
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as CSV

# stores the parameters for each simulation
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

# runs the simulation for various combinations of hybrid fitness and assortative mating strength
# stores the outcome of each simulation in a JLD2 file
function run_HZAM_set(set_name::String, intrinsic_R, ecolDiff;  # the semicolon makes the following optional keyword arguments 
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3,
    survival_fitness_method::String="epistasis", per_reject_cost=0, sigma_disp=0.05)
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

    # the set of hybrid fitnesses (w_hyb) values that will be run
    w_hyb_set = [1, 0.95, 0.9, 0.7, 0.5, 0.3, 0.1, 0]#[1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    # the set of assortative mating strengths (S_AM) that will be run
    S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

    total_functional_loci = union(female_mating_trait_loci, male_mating_trait_loci, competition_trait_loci, hybrid_survival_loci)
    neutral_loci = setdiff((1:total_loci), total_functional_loci)
    num_neutral_loci = length(neutral_loci)
    num_functional_loci = length(total_functional_loci)

    if survival_fitness_method == "epistasis"
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        short_survFitnessMethod = "Het"
    end

    # set up array of strings to record outcomes
    outcome_array = Array{OutputData,2}(undef, length(w_hyb_set), length(S_AM_set))

    sim_params = Array{SimulationParameters,2}(undef, length(w_hyb_set), length(S_AM_set))


    dir = mkpath(string("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/", set_name))

    # Loop through the different simulation sets
    Threads.@threads for i in eachindex(w_hyb_set)
        for j in eachindex(S_AM_set)
            w_hyb = w_hyb_set[i]
            S_AM = S_AM_set[j]
            println("ecolDiff = ", ecolDiff, "; w_hyb = ", w_hyb, "; S_AM = ", S_AM)

            run_name = string("HZAM_animation_run", "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", num_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM)

            # set up initial values for one simulation
            extinction = false
            K = K_total

            #K = (1+ecolDiff) * K_total

            # run one simulation by calling the function defined above:
            outcome = run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;
                K_total=Int(trunc(K)), max_generations,
                total_loci, female_mating_trait_loci, male_mating_trait_loci,
                competition_trait_loci, hybrid_survival_loci,
                survival_fitness_method, per_reject_cost, sigma_disp, do_plot=false)

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

            filename = string(dir, "/", run_name, ".jld2")
            @save filename parameters outcome

            println(run_name, "  outcome was: ", outcome)
            outcome_array[i, j] = outcome
            sim_params[i, j] = parameters
        end # of S_AM loop
    end # of w_hyb loop   

    return outcome_array, sim_params
end


#### functions for summarizing and plotting results

# creates a plot of all the outcomes in a folder with hybrid fitness on the x axis, strength of assortative mating on the y axis
# and the colour of each point representing the output variable of interest
function plot_output_field(outcomes, sim_params, run_name, fieldname)
    output = [[getfield(outcome, fieldname) for outcome in outcomes]...]
    w_hybs = [[s.w_hyb for s in sim_params]...]
    S_AMs = [[s.S_AM for s in sim_params]...]

    fieldname = String(fieldname)

    cr = (0, 1) # color range
    if cmp(fieldname, "variance") == 0
        cr = (0, 0.05)
    end

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60, colorrange=cr)
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string(fieldname, "--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    xlims!(-0.1, 1.1)
    points = scatter!(ax, w_hybs, S_AMs, color=output, markersize=50, colorrange=cr) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="overlap", height=Relative(0.5))

    fig
end

# creates and saves a set of plots based on the outcomes stored in a given folder
function make_and_save_figs(ResultsFolder, run_name, outcomes, sim_params)
    dir = mkpath(string(ResultsFolder, "/plots/", run_name))

    for field in fieldnames(OutputData)
        if cmp(String(field), "gene_flows") == 0
            for loci_type in fieldnames(GeneFlows)
                save(string(dir, "/", run_name, "_", loci_type, ".png"), plot_gene_flow(outcomes, sim_params, run_name, loci_type))
            end
        elseif cmp(String(field), "cline_widths") != 0 && !(occursin("position", String(field)))
            save(string(dir, "/", run_name, "_", String(field), ".png"), plot_output_field(outcomes, sim_params, run_name, field))
        end
    end
end

# creates a plot of gene flow (for a given trait) vs hybrid fitness on the x axis and assortative mating strength on the y axis
function plot_gene_flow(outcomes, sim_params, run_name, loci_type)
    all_gene_flow = [map(o -> o.gene_flows, outcomes)...]
    gene_flow_at_loci = [getfield(gf, loci_type) for gf in all_gene_flow]
    w_hybs = [[s.w_hyb for s in sim_params]...]
    S_AMs = [[s.S_AM for s in sim_params]...]

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60, colorrange=(0, 0.3))
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel="S_AM", title=string(loci_type, "--", run_name), yscale=log10, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    xlims!(-0.1, 1.1)

    points = scatter!(ax, w_hybs, S_AMs, color=gene_flow_at_loci, markersize=50, colorrange=(0, 0.3)) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="gene flow", height=Relative(0.5))

    fig
end

# creates a plot of all the outcomes in a csv file where the x axis is the strength of assortative mating
# and the colour of the points is the hybrid fitness
function plot_trait_vs_S_AM(filepath, name)
    w_hybs, S_AMs, output = load_from_file_csv(filepath)

    w_hyb_values = sort(union(w_hybs))
    S_AM_values = sort(union(S_AMs))

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60, colorrange=(0, 0.3))
    ax = Axis(fig[1, 1], xlabel="S_AM", ylabel=name, title=name, xticklabelsize=45, yticklabelsize=45, titlegap=30, xscale=log10) # creates the axes and labels

    ylims!(-0.04, 0.3)
    points = []
    xs = S_AM_values

    for w_hyb in w_hyb_values
        indices = filter(i -> w_hybs[i] == w_hyb, eachindex(output))
        xs = S_AMs[indices]
        ys = output[indices]
        push!(points, scatter!(xs, ys, markersize=30))
    end
    Legend(fig[1, 2],
        points,
        string.(w_hyb_values))

    display(fig)
    readline()
end

# plots all of the outcomes in a csv file vs w_hyb (ignores S_AM)
function plot_trait_vs_w_hyb(filepath, name)
    w_hybs, S_AMs, output = load_from_file_csv(filepath)

    w_hyb_values = sort(union(w_hybs))

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60, colorrange=(0, 0.3))
    ax = Axis(fig[1, 1], xlabel="w_hyb", ylabel=name, title=name, xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    ylims!(-0.04, 0.3)

    points = scatter!(w_hybs, output, markersize=30)
    display(fig)
    readline()
end

# This function adds jitter (small random shifts in position, to better visualize overlapping points)
function jitter(n, factor=0.005)
    n .+ (0.5 .- rand(length(n))) .* factor
end

# creates and saves to file a series of plots showing how hybrid fitness and assortative mating strength influence gene flow for each trait
function make_and_save_figs_gene_flow(ResultsFolder, run_name, outcomes, sim_params)
    dir = mkpath(string(ResultsFolder, "/plots/", run_name))

    for loci_type in fieldnames(GeneFlows)
        save(string(dir, "/", run_name, "_", loci_type, ".png"), plot_gene_flow(outcomes, sim_params, run_name, loci_type))
    end
end

# loads the data from each simulation output file stored in the given folder into an array of outcomes and a matching array of
# the simulation parameters
function load_from_folder(OutcomesFolder)
    simulation_outcomes = OutputData[]
    simulation_parameters = SimulationParameters[]
    files = readdir(OutcomesFolder)
    for file in files
        parameters = missing
        outcome = missing
        path = string(OutcomesFolder, "/", file)
        if occursin(".jld2", path)
            @load path parameters outcome
            push!(simulation_outcomes, outcome)
            push!(simulation_parameters, parameters)
        end
    end

    simulation_parameters, simulation_outcomes
end

# reads a csv and returns the data in a format that's easy to plot
# (a vector of the hybrid fitnesses, a vector of the assortative mating strengths, and a vector of the output data)
function load_from_file_csv(filepath)
    df = DataFrame(CSV.File(filepath))

    return df[!, "w_hyb"], df[!, "S_AM"], df[:, 3]
end

# converts a list of outcomes to a csv file for the given trait
function convert_to_CSV(sim_params, outcomes, output_field, output_folder)
    output = [[getfield(outcome, output_field) for outcome in outcomes]...]
    w_hybs = [[s.w_hyb for s in sim_params]...]
    S_AMs = [[s.S_AM for s in sim_params]...]

    df = DataFrame(S_AM=S_AMs, w_hyb=w_hybs, output_field=output)

    dir = mkpath(string("HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/CSV_data/", output_folder))

    CSV.write(string(dir, "/", String(output_field), ".csv"), df)
end

# creates csv files of the gene flow for each trait
function convert_to_CSV_gene_flows(sim_params, outcomes, output_folder)
    gene_flows = [[o.gene_flows for o in outcomes]...]
    for loci_type in fieldnames(GeneFlows)
        convert_to_CSV(sim_params, gene_flows, loci_type, output_folder)
    end
end

# creates csv files of the cline width for each trait
function convert_to_CSV_cline_widths(sim_params, outcomes, output_folder)
    cline_widths = [[o.cline_widths for o in outcomes]...]
    for loci_type in fieldnames(DataAnalysis.ClineWidths)
        convert_to_CSV(sim_params, cline_widths, loci_type, output_folder)
    end
end

# saves the genotypes to a file given the population data at the end of a simulation
function save_genotypes(pd, filepath)
    genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
    #filepath = string("genotypes.jld2")
    @save filepath genotypes
end

# calculates the Pearson coefficient between two loci
function calc_linkage_diseq(genotypes, l1, l2)
    if l2 != l1
        genotypes = [g[:, [l1, l2]] for g in genotypes]

        haplotypes = vcat([g[1, :] for g in genotypes], [g[2, :] for g in genotypes])

        p_A = count(h -> h[1] == 0, haplotypes) / length(haplotypes)
        p_B = count(h -> h[2] == 0, haplotypes) / length(haplotypes)

        p_AB = count(h -> h == [0, 0], haplotypes) / length(haplotypes)
        D = (p_AB - (p_A * p_B))
        pearson_coefficient = (D^2) / (p_A * (1 - p_A) * p_B * (1 - p_B))

        if isnan(pearson_coefficient)
            pearson_coefficient = 1
        end
        return pearson_coefficient
    else
        return 1
    end
end


# calculates the linkage disequilibrium between loci using Pearson coefficients
# and returns a table of values representing the average correlation between two traits
function calc_linkage_diseq_all(path, plot_title)

    @load path genotypes


    num_loci = size(genotypes[1], 2)

    rows = (1:num_loci)
    cols = (1:num_loci)'

    calc_linkage_diseq(genotypes, 5, 6)

    linkage_diseq = calc_linkage_diseq.(Ref(genotypes), rows, cols)

    traits = ["f mating", "m mating", "competition", "hybrid survival", "neutral"]

    loci = [1:4, 5:8, 9:12, 13:16, 17:20]

    function average_linkage_diseq(l1, l2)
        return l1, l2, mean(linkage_diseq[loci[l1], loci[l2]])
    end

    avg_linkage_diseq = [average_linkage_diseq.((1:5), (1:5)')...]

    xs = [a[1] for a in avg_linkage_diseq]
    ys = [a[2] for a in avg_linkage_diseq]
    output = [a[3] for a in avg_linkage_diseq]

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="trait", ylabel="trait", title=plot_title, xticklabelsize=45, yticklabelsize=45, xticks = (1:5, traits), yticks = (1:5, traits), xticklabelrotation = pi/2, titlegap=30) # creates the axes and labels


    points = scatter!(ax, xs, ys, color=output, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="Pearson coefficient", height=Relative(1.0))

    display(fig)
    dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/gene_linkages")
    filepath = string(dir, "/", plot_title,".png")

    save(filepath, fig)
end

# creates and saves a graph showing the competition trait cline
function check_competition_trait(path)
    @load path genotypes 
    hybrid_indices = calc_traits_additive(genotypes, 9:12)
    sort!(hybrid_indices)

    
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="competition_trait", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices,markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/competition_trait-", path,".png"), fig)
end

# creates and saves a graph showing the male mating trait cline
function check_male_mating_trait(path)
    @load path genotypes 
    hybrid_indices = calc_traits_additive(genotypes, 5:8)
    sort!(hybrid_indices)

    
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="male_mating_trait", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices,markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/male_mating_trait-", path,".png"), fig)
end