using Colors, ColorSchemes
import ColorSchemes.plasma
using GLMakie
using Statistics: mean
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as CSV

"The directory where all files are saved"
global results_folder = "HZAM_Sym_Julia_results_GitIgnore"

"""
    set_results_folder(dir::String)

Set the output directory for all methods in summarize_data.jl.
"""
function set_results_folder(dir::String)
    global results_folder = dir
end

"""
    run_HZAM_set(
        set_name::String, 
        intrinsic_R::Real,
        ecolDiff::Real; 
        <keyword arguments>
    )

Run the simulation for 64 combinations of hybrid fitness and assortative mating strength and 
store the outcome of each simulation in a JLD2 file.

# Arguments
- `set_name::String`: the name assigned to the set of simulations.
- `intrinsic_R::Real`: the intrinsic growth rate.
- `ecolDiff::Real`: the ecological difference between the two species. Between 0 and 1.
- `K_total::Integer=20000`: the carrying capacity of the environment.
- `max_generations::Integer=1000`: the number of generations that the simulation will run for.
- `total_loci::Integer=6`: the total number of loci in the genome.
- `female_mating_trait_loci=1:3`: the loci specifying the female's mate preference.
- `male_mating_trait_loci=1:3`: the loci specifying the male's mating trait.
- `competition_trait_loci=1:3`: the loci specifying the ecological trait (used in fitness related to resource use).
- `hybrid_survival_loci=1:3`: the loci specifying the probability of survival to adulthood.
- `survival_fitness_method:String="epistasis"`: the method used to calculate the probability of survival to adulthood.
- `per_reject_cost=0`: the fitness loss of female per male rejected (due to search time, etc.). Can take values of 0 to 1.
- `sigma_disp=0.05`: the standard deviation of the normal distribution determining how far offspring will disperse from their mothers.
- `track_spatial_data=false`: if true, keeps track of bimodality, gene flow, variance, cline widths, and cline positions.
- `track_population_data=true`: if true, keeps track of population size, hybridness, overlap, and hybrid zone width.
- `track_fitness=false`: if true, keeps track of the average number of offspring for each phenotype.
- `track_mating_success=false`: if true, keeps track of the average number of mates per male for each phenotype.
- `track_phenotypes=false`: if true, keeps track of the number of individuals with each phenotype.
"""
function run_HZAM_set(set_name::String, intrinsic_R::Real, ecolDiff::Real;  # the semicolon makes the following optional keyword arguments 
    K_total::Int=20000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3,
    survival_fitness_method::String="epistasis", per_reject_cost=0, sigma_disp=0.05,
    track_spatial_data=false, track_population_data=true, track_fitness=false,
    track_mating_success=false, track_phenotypes=false)

    # the set of hybrid fitnesses (w_hyb) values that will be run
    w_hyb_set = [1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5]
    S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific


    # set up array of strings to record outcomes
    outcome_array = Array{DataAnalysis.OutputData,2}(undef, length(w_hyb_set), length(S_AM_set))

    dir = mkpath(string(results_folder, "/simulation_outcomes/", set_name))

    # Loop through the different simulation sets
    Threads.@threads for i in eachindex(w_hyb_set)
        for j in eachindex(S_AM_set)
            w_hyb = w_hyb_set[i]
            S_AM = S_AM_set[j]
            println("ecolDiff = ", ecolDiff, "; w_hyb = ", w_hyb, "; S_AM = ", S_AM)

            run_name = string("HZAM_animation_run", "_ecolDiff", ecolDiff, "_gen", max_generations, "_SC", per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM)

            # run one simulation by calling the function defined above:
            outcome = run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;
                K_total, max_generations,
                total_loci, female_mating_trait_loci, male_mating_trait_loci,
                competition_trait_loci, hybrid_survival_loci,
                survival_fitness_method, per_reject_cost, sigma_disp, do_plot=false,
                track_spatial_data, track_population_data, track_fitness,
                track_mating_success, track_phenotypes)



            filename = string(dir, "/", run_name, ".jld2")
            @save filename outcome

            println(run_name, "  completed.")
            outcome_array[i, j] = outcome
        end # of S_AM loop
    end # of w_hyb loop   

    return outcome_array
end

"""
    plot_output_field(
        outcomes::Array{:Real},
        sim_params::Array{<:DataAnalysis.SimParams}
    )

Create a heatmap of an output variable vs hybrid fitness and assortative mating.

# Arguments
- `outcomes::Array{:Real}`: the output from the simulation to be displayed.
- `sim_params::Array{<:DataAnalysis.SimParams}`: the simulation parameters resulting in the outcomes.
"""
function plot_output_field(
    outcomes::Array{:Real},
    sim_params::Array{<:DataAnalysis.SimParams}
)
    output = [outcomes...]
    w_hybs = [[s.w_hyb for s in sim_params]...]
    S_AMs = [[s.S_AM for s in sim_params]...]

    w_hyb_set = sort(union(w_hybs))
    S_AM_set = sort(union(S_AMs))
    xticks = collect(1:length(w_hyb_set))
    yticks = collect(1:length(S_AM_set))

    xs = map(w -> indexin(w, w_hyb_set)[1], w_hybs)
    ys = map(s -> indexin(s, S_AM_set)[1], S_AMs)


    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)  # this sets the standard font size

    fig = Figure(resolution=(1800, 1200), figure_padding=60, colormap=:grayC)

    ax = Axis(
        fig[1, 1],
        xlabel="w_hyb",
        ylabel="S_AM",
        xticks=(xticks, string.(w_hyb_set)),
        yticks=(yticks, string.(S_AM_set)),
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=15,
        ylabelsize=15
    )
    hm = heatmap!(ax, xs, ys, output)

    Colorbar(fig[:, 2], hm, ticklabelsize=15)
    display(fig)
    readline()
    return fig
end

"""
    load_from_folder(dir::String)

Load the data from each simulation output file stored in the given directory into an array 
organized by hybrid fitness and assortative mating strength.
"""
function load_from_folder(dir::String)
    files = readdir(dir)
    w_hyb_set = [1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5]
    S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]
    outcome_array = Array{DataAnalysis.OutputData,2}(
        undef, length(w_hyb_set), length(S_AM_set)
    )
    for file in files
        path = string(dir, "/", file)
        if occursin(".jld2", path)
            @load path outcome
            outcome_array[
                indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
                indexin(outcome.sim_params.S_AM, S_AM_set)[1]
            ] = outcome
        end
    end
    return outcome_array
end

"""
    load_from_csv(filepath::String)

Read the data from a CSV file and return vectors for hybrid fitness, assortative mating 
strength, and the output variable.
"""
function load_from_csv(filepath::String)
    df = DataFrame(CSV.File(filepath))

    return df[!, "w_hyb"], df[!, "S_AM"], df[:, 3]
end

"""
    convert_to_CSV(
        outcome_array::Array{<:Real},
        w_hyb_array::Array{<:Real},
        S_AM_array::Array{<:Real},
        name::String
    )

Convert an array of outcomes to a csv file for the given output field.

# Arguments
- `outcome_array::Array{<:Real}`: the array of the output data of interest.
- `w_hyb_array::Array{<:Real}`: the array of the w_hyb parameters used.
- `S_AM_array::Array{<:Real}`: the array of the S_AM parameters used.
- `field_name::Symbol`: the name of the output field of interest.
- `name::String`: the name for the CSV.
"""
function convert_to_CSV(
    outcome_array::Array{<:Real},
    w_hyb_array::Array{<:Real},
    S_AM_array::Array{<:Real},
    name::String
)
    output = [outcome_array...]
    w_hybs = [w_hyb_array...]
    S_AMs = [S_AM_array...]

    df = DataFrame(S_AM=S_AMs, w_hyb=w_hybs, output_field=output)

    dir = mkpath(string(results_folder, "/simulation_outcomes/CSV_data/"))

    CSV.write(string(dir, "/", name, ".csv"), df)
end

"""
    plot_fitnesses(fitnesses::Vector{<:Dict})

Produce a plot of fitnesses per phenotype over time.

# Arguments
- `fitnesses::Vector{<:Dict}`: number of offspring per phenotype each generation.
"""
function plot_fitnesses(fitnesses::Vector{<:Dict})
    dir = mkpath(string(results_folder, "/plots/"))
    @save string(dir, "fitness.jld2") fitnesses
    generations = collect(eachindex(fitnesses))
    num_generations = length(generations)
    phenotypes = collect(keys(fitnesses[1]))
    num_phenotypes = length(fitnesses[1])
    xs = vcat([fill(gen, num_phenotypes) for gen in eachindex(fitnesses)]...)
    ys = vcat([collect(keys(f)) for f in fitnesses]...)
    zs = vcat([collect(values(f)) for f in fitnesses]...)

    """
    Determines which phenotype has the highest fitness each generation.
    """
    function maximum_fitness()
        return map(
            f -> reduce((x, y) -> f[x] ≥ f[y] ? x : y, keys(f)),
            fitnesses
        )
    end

    """
    Compute the running average of the fitnesses for a given phenotype.
    """
    function average_fitness(phenotype)
        n = 10
        average_fitnesses = []
        for i in eachindex(fitnesses)
            fitness_values =
                [fitnesses[j][phenotype] for j in i:min(length(fitnesses), i + n - 1)]
            push!(average_fitness, mean(fitness_values))
        end
        return average_fitnesses
    end

    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    ax = Axis(fig[1, 1], xlabel="generation", ylabel="phenotype")
    scatter!(ax, xs, ys)
    display(fig)
    display(heatmap(xs, ys, zs, colormap=:grayC))
    readline()
end

"""
    plot_population_tracking_data(filepath::String)

Create plots of the population size, hybrid zone width, hybrid index, and population overlap vs time.

# Arguments
- `filepath::String`:: the filepath for the file containing the outcome of the simulation.
"""
function plot_population_tracking_data(filepath::String)
    @load filepath population_tracking_data
    overlaps = [t.overlap for t in population_tracking_data]
    hybridnesses = [t.hybridness for t in population_tracking_data]
    widths = [t.width for t in tracking_data]
    populations = [t.population for t in population_tracking_data]

    xs = collect(eachindex(population_tracking_data))
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
    
    dir = mkpath(string(results_folder, "/plots/"))

    save(string(
        dir,
        "hz_formation_ecolDiff",
        sim_params.ecolDiff,
        "_S_AM",
        sim_params.S_AM,
        "_w_hyb",
        sim_params.w_hyb,
        ".png",
        fig
    ))

    readline()
end

"""
    summarize_gene_correlations(source_dir::String; filename="gene_correlations")

Create plots of the gene correlations between different traits vs S_AM and w_hyb for 
different combinations of the same loci controlling different traits.

source_dir should point to a folder containing magic_preference.jld2, magic_cue.jld2, 
search_cost.jld2, and no_magic.jld2. Each file should contain an outcome array from a simulation set.

# Arguments
- `source_dir::String`: the folder where the outcomes of all the simulations are stored.
- `file_name="gene_correlations`: the name of the image file for the plot. 
"""
function summarize_gene_correlations(source_dir::String; filename="gene_correlations")
    names = ["magic_preference", "magic_cue", "search_cost", "no_magic"]
    xlabelnames = ["Female mating trait", "Male mating trait", "Neutral trait"]
    ylabelnames = ["Magic preference", "Magic cue", "Search cost", "No pleiotropy"]
    xlabels = []
    ylabels = []
    correlations = Matrix(undef, 4, 3)

    magic_loci = [[1, 5], [3, 5], [5, 6], [5, 6]]
    fmt_loci = [[2], 1:2, 1:2, 1:2]
    mmt_loci = [3:4, [4], 3:4, 3:4]
    neutral_loci = [6:9, 6:9, 7:10, 7:10]

    for i in 1:4
        filename = string(source_dir, "/", names[i], ".jld2")
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

    for i in 1:3
        push!(xlabels, Label(fig[0, i], xlabelnames[i], tellwidth=false))
    end

    for i in 1:4
        push!(ylabels, Label(fig[i, 0], ylabelnames[i], rotation=pi / 2, tellheight=false))
    end

    Colorbar(fig[:, 4], hm[1, 1], label="Pearson coefficient", ticklabelsize=15)
    display(fig)
    dir = mkpath(string(results_folder, "/plots/"))
    save(string(dir, "/", filename, ".png"), fig)
    readline()
end

"""
    plot_phenotypes(phenotypes; filename="phenotype_frequencies")

Create plots of mating preference, mating cue, and hybrid survival trait phenotype 
frequencies over time.

# Arguments
- `filename="phenotype_frequencies"`: the name for the saved plot image.
"""
function plot_phenotypes(phenotypes; filename="phenotype_frequencies")
    ylabels = ["Mating preference", "Mating cue", "Hybrid survival trait"]
    function calc_proportion(dict)
        return_dict = Dict(collect(keys(dict)) .=> Float32.(collect(values(dict))))
        counts = collect(values(dict))
        sum_counts = sum(counts)
        map!(x -> x / sum_counts, values(return_dict))
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
        ax[i] = Axis(fig[i, 1], xlabel="generation", ylabel=ylabels[i])
    end

    hm[1] = plot_phenotypes_at_loci(:female_mating_trait, ax[1])
    hm[2] = plot_phenotypes_at_loci(:male_mating_trait, ax[2])
    hm[3] = plot_phenotypes_at_loci(:hybrid_survival, ax[3])
    display(fig)
    
    dir = mkpath(string(results_folder, "/plots/"))
    save(string(dir, filename, ".png"), fig)
    readline()
end