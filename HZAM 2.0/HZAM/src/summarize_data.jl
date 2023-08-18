using Colors, ColorSchemes
import ColorSchemes.plasma
using GLMakie
using Statistics: mean
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

function t_test(genotypes)
    function calc_heterozygosity(locus)
        count(x -> x == [0; 1] || x == [1; 0], [g[:, locus] for g in genotypes]) / length(genotypes)
    end

    heterozygosities = map(calc_heterozygosity, (1:20))

    neutral_heterozygosities = heterozygosities[17:20]

    loci = [1:4, 5:8, 9:12, 13:16]

    for range in loci
        println(EqualVarianceTTest(heterozygosities[range], neutral_heterozygosities))
    end
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
            println(string("loci 1: ", l1, " loci 2: ", l2))
            println((p_A * (1 - p_A) * p_B * (1 - p_B)))
            pearson_coefficient = 1
        end
        return pearson_coefficient
    else
        return 1
    end
end

function calc_chi_squared(genotypes, locus)
    genotypes = [g[:, locus] for g in genotypes]
    alleles = vcat([g[1] for g in genotypes], [g[2] for g in genotypes])

    N = length(genotypes)
    p_A = count(x -> x == 0, alleles) / (2 * N)
    p_B = 1 - p_A
    n_AB = count(x -> x == [0; 1] || x == [1; 0], genotypes)
    n_AA = count(x -> x == [0; 0], genotypes)
    n_BB = count(x -> x == [1; 1], genotypes)


    return (((n_AA - N * p_A^2)^2) / (N * p_A^2)) +
           (((n_AB - 2 * N * p_A * p_B)^2) / (2 * N * p_A * p_B)) +
           (((n_BB - N * p_B^2)^2) / (N * p_B^2))
end

function load_genotypes(path)
    @load path genotypes

    return genotypes
end

# calculates the linkage disequilibrium between loci using Pearson coefficients
# and returns a table of values representing the average correlation between two traits
function plot_trait_correlations_all(genotypes, loci, plot_title)

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)

    xs = [t[1] for t in trait_correlations]
    ys = [t[2] for t in trait_correlations]
    output = [t[3] for t in trait_correlations]

    xs = map(x -> findall(z -> z == x, keys(loci))[1], xs)
    ys = map(y -> findall(z -> z == y, keys(loci))[1], ys)

    ticks = (1:length(loci))
    labels = [string.(keys(loci))...]

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="trait", ylabel="trait", title=plot_title, xticklabelsize=45, yticklabelsize=45, xticks=(ticks, labels), yticks=(ticks, labels), xticklabelrotation=pi / 2, titlegap=30) # creates the axes and labels


    points = scatter!(ax, xs, ys, color=output, markersize=50) # adds the location of every individual to the plot

    Colorbar(fig[1, 4], points, label="Pearson coefficient", height=Relative(1.0))

    display(fig)
    dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/trait_correlations_ecolDiff1")
    filepath = string(dir, "/", plot_title, ".png")

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
    ax = Axis(fig[1, 1], title="competition_trait, ecolDiff=1", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices, markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/competition_trait-ecolDiff1_w_hyb_0.8_S_AM100.png"), fig)
end

# creates and saves a graph showing the competition trait cline
function check_neutral_trait(path)
    @load path genotypes
    hybrid_indices = calc_traits_additive(genotypes, 17:20)
    sort!(hybrid_indices)


    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="neutral trait, ecolDiff=1", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices, markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/neutral_trait_ecolDiff_1_w_hyb_0.8_S_AM100.png"), fig)
end

# creates and saves a graph showing the male mating trait cline
function check_male_mating_trait(path)
    @load path genotypes
    hybrid_indices = calc_traits_additive(genotypes, 5:8)
    sort!(hybrid_indices)


    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="male_mating_trait, ecolDiff=1", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices, markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/male_mating_trait_ecolDiff_1_w_hyb_0.8_S_AM100.png"), fig)
end


# creates and saves a graph showing the male mating trait cline
function check_female_mating_trait(path)
    @load path genotypes
    hybrid_indices = calc_traits_additive(genotypes, 5:8)
    sort!(hybrid_indices)


    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="female_mating_trait, ecolDiff=1", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices, markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/female_mating_trait_ecolDiff_1_w_hyb_0.8_S_AM100.png"), fig)
end


# creates and saves a graph showing the male mating trait cline
function check_hybrid_survival(path)
    @load path genotypes
    hybrid_indices = calc_traits_additive(genotypes, 5:8)
    sort!(hybrid_indices)


    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], title="hybrid_survival, ecolDiff=1", xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels

    points = scatter!(ax, collect(eachindex(hybrid_indices)), hybrid_indices, markersize=20) # adds the location of every individual to the plot

    display(fig)
    save(string("HZAM_Sym_Julia_results_GitIgnore/gene_linkages/hybrid_survival_ecolDiff_1_w_hyb_0.8_S_AM100.png"), fig)
end

function find_extinct_alleles(genotypes)

    extinct = []
    for i in 1:size(genotypes[1], 2)
        for j in 0:1
            n = count(x -> x[1, i] == j || x[2, i] == j, genotypes)
            if n == 0
                push!(extinct, i)
            end
        end
    end
    return extinct
end


global ax, points, fig

function create_gene_plot(
    genotypes,
    loci,
    generation,
    save_plot
)
    global ax = Vector(undef, 6)
    global points = Vector(undef, 6)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # set the standard font size
    global fig = Figure(resolution=(1800, 1200), figure_padding=60)
    # create the axis and labels
    global ax[2] = Axis(
        fig[1, 1],
        xlabel="trait value",
        ylabel="# individuals",
        title="Female mating trait",
        limits=(-0.05, 1.05, 0, 6000),
        yticks=(0:1000:6000))
    global ax[3] = Axis(
        fig[1, 2],
        xlabel="trait value",
        ylabel="# individuals",
        title="Male mating trait",
        limits=(-0.05, 1.05, 0, 12000),
        yticks=(0:3000:12000, ["0", "3000", "6000", "9000", "12000"])
    )
    global ax[4] = Axis(
        fig[1, 3],
        xlabel="trait value",
        ylabel="# individuals",
        title="Competition trait",
        limits=(-0.05, 1.05, 0, 5000),
        yticks=(0:1000:5000)
    )
    global ax[5] = Axis(
        fig[2, 1],
        xlabel="trait value",
        ylabel="# individuals",
        title="Hybrid survival trait",
        limits=(-0.05, 1.05, 0, 12000),
        yticks=(0:3000:12000, ["0", "3000", "6000", "9000", "12000"])
    )
    global ax[1] = Axis(
        fig[2, 2],
        xlabel="trait value",
        ylabel="# individuals",
        title="Neutral trait",
        limits=(-0.05, 1.05, 0, 5000),
        yticks=(0:1000:5000)
    )

    global ax[6] = Axis(
        fig[2, 3],
        xlabel="trait value",
        ylabel="# individuals",
        title="Expected neutral",
        limits=(-0.05, 1.05, 0, 5000),
        yticks=(0:1000:5000)
    )

    function get_hybrid_indices(range, genotypes)
        hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, range)
        return hybrid_indices
    end

    function get_points(range, genotypes)
        num_loci = length(range)
        indices_range = (0:(1/(2*num_loci)):1)
        hybrid_indices = get_hybrid_indices(range, genotypes)
        return map(i -> count(x -> x == i, hybrid_indices), indices_range)
    end

    function get_expected(i)
        n = length(genotypes)
        return (n * binomial(8, Int(8 * i))) / 2^8
    end

    hybrid_survival_indices = get_hybrid_indices(loci.hybrid_survival, genotypes)

    species_A_indices = findall(h -> h == 0, hybrid_survival_indices)
    species_B_indices = findall(h -> h == 1, hybrid_survival_indices)
    other_indices = setdiff(eachindex(genotypes), union(species_A_indices, species_B_indices))

    species_A_output = get_points.(values(loci), Ref(genotypes[species_A_indices]))
    species_B_output = get_points.(values(loci), Ref(genotypes[species_B_indices]))
    other_output = get_points.(values(loci), Ref(genotypes[other_indices]))

    for i in 3:7
        num_loci = length(values(loci)[i])
        range = (0:(1/(2*num_loci)):1)
        xs = [range; range; range]
        stk = [fill(1, length(species_A_output[i])); fill(3, length(species_B_output[i])); fill(2, length(other_output[i]))]
        ys = [species_A_output[i]; species_B_output[i]; other_output[i]]

        global points[i-2] = barplot!(
            ax[i-2],
            xs,
            ys,
            color=stk,
            stack=stk,
            colorrange=(1.1, 2.9),
            highclip=:orange,
            lowclip=:blue)
    end
    global points[6] = barplot!(ax[6], (0:0.125:1), get_expected.((0:0.125:1)), color=:blue)

    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse5")
        save(string(dir, "/", generation, ".png"), fig)
    end
    display(fig)
end

function update_gene_plot(
    genotypes,
    loci,
    generation,
    save_plot
)

    [delete!(ax[i], points[i]) for i in 1:6] # remove the old points from the plot


    function get_hybrid_indices(range, genotypes)
        hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, range)
        return hybrid_indices
    end

    function get_points(range, genotypes)
        num_loci = length(range)
        indices_range = (0:(1/(2*num_loci)):1)
        hybrid_indices = get_hybrid_indices(range, genotypes)
        return map(i -> count(x -> x == i, hybrid_indices), indices_range)
    end

    function get_expected(i)
        n = length(genotypes)
        return (n * binomial(8, Int(8 * i))) / 2^8
    end

    hybrid_survival_indices = get_hybrid_indices(loci.hybrid_survival, genotypes)


    species_A_indices = findall(h -> h == 0, hybrid_survival_indices)
    species_B_indices = findall(h -> h == 1, hybrid_survival_indices)
    other_indices = setdiff(eachindex(genotypes), union(species_A_indices, species_B_indices))

    species_A_output = get_points.(values(loci), Ref(genotypes[species_A_indices]))
    species_B_output = get_points.(values(loci), Ref(genotypes[species_B_indices]))
    other_output = get_points.(values(loci), Ref(genotypes[other_indices]))

    for i in 3:7
        num_loci = length(values(loci)[i])
        range = (0:(1/(2*num_loci)):1)
        xs = [range; range; range]
        stk = [fill(1, length(species_A_output[i])); fill(3, length(species_B_output[i])); fill(2, length(other_output[i]))]
        ys = [species_A_output[i]; species_B_output[i]; other_output[i]]

        global points[i-2] = barplot!(
            ax[i-2],
            xs,
            ys,
            color=stk,
            stack=stk,
            colorrange=(1.1, 2.9),
            highclip=:orange,
            lowclip=:blue)
    end
    global points[6] = barplot!(ax[6], (0:0.125:1), get_expected.((0:0.125:1)), color=:blue)
    println(string("generation: ", generation))

    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse5")
        save(string(dir, "/", generation, ".png"), fig)
    end
end

function chi_squared_traits(genotypes, loci)
    n = length(genotypes)

    function get_hybrid_indices(range)
        hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, range)
        return sort(hybrid_indices)
    end

    neutral_hybrid_indices = get_hybrid_indices(loci.neutral)

    function get_q(range)
        haplotypes1 = [g[1, range] for g in genotypes]
        haplotypes2 = [g[2, range] for g in genotypes]
        (sum(map(h -> count(x -> x == 1, h), haplotypes1)) + sum(map(h -> count(x -> x == 1, h), haplotypes2))) / (n * 8)
    end

    function calc_probability(trait_index, q)
        return (q^trait_index) * ((1 - q)^(8 - trait_index)) * binomial(8, Int(trait_index))
        binomial(8, Int(8 * trait_index)) / (2^8)
    end

    function count_indv(i, hybrid_indices)
        count(x -> x == i, hybrid_indices)
    end

    function calc_expected(i, q)
        # return n * calc_probability(8*i, q)
        return count(x -> x == i, neutral_hybrid_indices)
    end

    function calc_phenotype(i, q, indices)
        m = calc_expected(i, q)
        x = count_indv(i, indices)
        ((x - m)^2) / m
    end

    function check_trait(trait)
        q = get_q(loci[trait])
        indices = get_hybrid_indices(loci[trait])
        intermediate = map(x -> calc_phenotype(x, q, indices), collect(0:0.125:1))
        output = 0
        return sum(intermediate)
    end
    println(keys(loci)[3:7])
    println(map(check_trait, keys(loci)[3:7]))
end

function plot_fitnesses(fitnesses)
    @save "fitness.jld2" fitnesses
    generations = collect(eachindex(fitnesses))
    num_generations = length(generations)
    phenotypes = collect(keys(fitnesses[1]))
    num_phenotypes = length(fitnesses[1])
    xs = vcat([fill(gen, num_phenotypes) for gen in eachindex(fitnesses)]...)
    ys = vcat([collect(keys(f)) for f in fitnesses]...)
    zs = vcat([collect(values(f)) for f in fitnesses]...)


    function maximum_fitness()
        return map(
            f -> reduce((x, y) -> f[x] ≥ f[y] ? x : y, keys(f)),
            fitnesses
        )
    end

    function ratio()
        return map(
            f -> f[0.5f0] / f[0.0f0],
            fitnesses
        )
    end

    function average_fitness(phenotype)
        average_fitnesses = []
        for i in 1:10:length(fitnesses)
            fitness_values = []
            for j in i:min(i + 9, length(fitnesses))
                push!(fitness_values, fitnesses[j][phenotype])
            end
            push!(average_fitnesses, mean(fitness_values))
        end
        return average_fitnesses
    end

    #=
    xs = vcat([(1:10:num_generations) for phenotype in phenotypes]...)
    ys = vcat([fill(phenotype, Int(trunc(num_generations/10))) for phenotype in phenotypes]...)
    zs = vcat(map(average_fitness, phenotypes)...)
    zs = Float32.(zs)
    =#
    #=
    xs = generations
    ys = ratio()
    =#

    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="generation", ylabel="phenotype")
    scatter!(ax, xs, ys)
    display(fig)
    display(heatmap(xs, ys, zs, colormap=:grayC))
    readline()
end

function plot_tracking_data(filepath)
    @load filepath tracking_data
    overlaps = [t.overlap for t in tracking_data]
    hybridnesses = [t.hybridness for t in tracking_data]
    widths = [t.width for t in tracking_data]
    populations = [t.population for t in tracking_data]

    xs = collect(eachindex(tracking_data))
    ax = Vector(undef, 4)
    points = Vector(undef, 4)

    fontsize_theme = Theme(fontsize=40)
    set_theme!(fontsize_theme)  # set the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    # create the axis and labels
    ax[1] = Axis(
        fig[1, 1],
        ylabel="# individuals",
        yticklabelsize=30.0f0
    )
    ax[2] = Axis(
        fig[2, 1],
        ylabel="cline width",
        limits=((nothing, nothing), (-0.05, 1.05))
    )
    ax[3] = Axis(
        fig[3, 1],
        ylabel="hybridness",
        yticklabelsize=30.0f0
    )
    ax[4] = Axis(
        fig[4, 1],
        ylabel="overlap",
        xlabel="generation",
        limits=((nothing, nothing), (-0.05, 1.05))
    )
    points[1] = scatter!(ax[1], xs, populations)
    points[2] = scatter!(ax[2], xs, widths)
    points[3] = scatter!(ax[3], xs, hybridnesses)
    points[4] = scatter!(ax[4], xs, overlaps)

    yspace = 1.5 * maximum(tight_yticklabel_spacing!, [ax[1], ax[2]])

    ax[1].yticklabelspace = yspace
    ax[2].yticklabelspace = yspace
    ax[3].yticklabelspace = yspace
    ax[4].yticklabelspace = yspace

    display(fig)
    readline()

    save("hz_formation_ecolDiff1_S_AM_350_w_hyb_0.87.png", fig)
end

function summarize_gene_correlations(dir)
    names = ["magic_preference", "magic_cue", "search_cost", "no_magic"]
    correlations = Matrix(undef, 4, 3)

    magic_loci = [[1, 5], [3, 5], [5, 6], [5, 6]]
    fmt_loci = [[2], 1:2, 1:2, 1:2]
    mmt_loci = [3:4, [4], 3:4, 3:4]
    neutral_loci = [6:9, 6:9, 7:10, 7:10]

    for i in 1:4
        filename = string(dir, "/", names[i], ".jld2")
        @load filename sim_params outcome_array
        correlations[i, 1] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], fmt_loci[i]), outcome_array)...)
        correlations[i, 2] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], mmt_loci[i]), outcome_array)...)

        correlations[i, 3] = vcat(map(genotypes -> DataAnalysis.calc_trait_correlation(genotypes, magic_loci[i], neutral_loci[i]), outcome_array)...)

    end

    w_hyb_set = string.([0.95, 0.9, 0.7, 0.5])
    S_AM_set = string.([1, 10, 100, 1000])
    xticks = [4, 3, 2, 1]
    yticks = [1, 2, 3, 4]

    xs = vcat(fill(xticks, 4)...)
    ys = vcat([fill(s, 4) for s in yticks]...)

    ax = Matrix(undef, 4, 4)
    hm = Matrix(undef, 4, 4)


    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)  # this sets the standard font size

    fig = Figure(resolution=(1800, 1200), figure_padding=60, colorrange=(0, 1), colormap=:curl)

    println(correlations)

    for j in 1:4
        for i in 1:3
            ax[j, i] = Axis(fig[j, i], xlabel="w_hyb", ylabel="S_AM", xticks=(xticks, w_hyb_set), yticks=(yticks, S_AM_set), xticklabelsize=20, yticklabelsize=20, xlabelsize=15, ylabelsize=15)
            hm[j, i] = heatmap!(ax[j, i], xs, ys, correlations[j, i], colorrange=(-1, 1))
            println(correlations[j, i][16])
        end
        println("")
    end

    xlabel1 = Label(fig[0, 1], "Female mating trait", tellwidth=false)
    xlabel2 = Label(fig[0, 2], "Male mating trait", tellwidth=false)
    xlabel3 = Label(fig[0, 3], "Neutral trait", tellwidth=false)

    ylabel1 = Label(fig[1, 0], "Magic preference", rotation=pi / 2, tellheight=false)
    ylabel2 = Label(fig[2, 0], "Magic cue", rotation=pi / 2, tellheight=false)
    ylabel3 = Label(fig[3, 0], "Search cost", rotation=pi / 2, tellheight=false)
    ylabel4 = Label(fig[4, 0], "No pleiotropy", rotation=pi / 2, tellheight=false)
    Colorbar(fig[:, 4], hm[1, 1], label="Pearson coefficient", ticklabelsize=15)
    display(fig)
    readline()
    save("trial1.png", fig)
end

function plot_phenotypes(phenotypes)
    ylabels = ["Mating preference", "Mating cue", "Hybrid survival trait"]
    function calc_proportion(dict)
        return_dict = Dict(collect(keys(dict)).=>Float32.(collect(values(dict))))
        counts = collect(values(dict))
        sum_counts = sum(counts)
        map!(x->x/sum_counts, values(return_dict))
        return return_dict
    end
    function plot_phenotypes_at_loci(loci, ax)
        proportions = [calc_proportion(p[loci]) for p in phenotypes]
        generations = collect(eachindex(proportions))
        num_phenotypes = length(proportions[1])
        xs = vcat([fill(gen, num_phenotypes) for gen in generations]...)
        ys = vcat([collect(keys(p)) for p in proportions]...)
        zs = vcat([collect(values(p)) for p in proportions]...)
        return heatmap!(ax, xs, ys, zs)
    end
    
    fig = Figure(resolution=(1800, 1200), figure_padding=60, colormap=:grayC)
    ax = Vector(undef, 3)
    hm = Vector(undef, 3)

    for i in 1:3
        ax[i] = Axis(fig[i,1], xlabel="generation", ylabel=ylabels[i])
    end

    hm[1] = plot_phenotypes_at_loci(:female_mating_trait, ax[1])
    hm[2] = plot_phenotypes_at_loci(:male_mating_trait, ax[2])
    hm[3] = plot_phenotypes_at_loci(:hybrid_survival, ax[3])
    display(fig)
    readline()
    
    save("phenotypes2.png", fig)
end