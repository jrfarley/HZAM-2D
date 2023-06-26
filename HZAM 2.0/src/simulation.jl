include("population.jl")

using .Population

using Distributions # needed for "Poisson" function
using Statistics  # needed for "mean" function
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as csv
using LsqFit


# This is the function to run a single HZAM simulation
function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3, neutral_loci=4:6,
    survival_fitness_method::String="epistasis", per_reject_cost=0,
    geographic_limits::Vector{Location}=[Location(0.0f0, 0.0f0), Location(0.999f0, 0.999f0)],
    starting_range_pop0=[Location(0.0f0, 0.0f0), Location(0.48f0, 1.0f0)], starting_range_pop1=[Location(0.52f0, 0.0f0), Location(1.0f0, 1.0f0)],
    sigma_disp=0.01, sigma_comp=0.01,
    do_plot=true, plot_int=10, optimize=true)

    functional_loci_range = setdiff(collect(1:total_loci), collect(neutral_loci))

    pop = []

    # get the chosen survival fitness function
    if survival_fitness_method == "epistasis"
        get_survival_fitness = get_survival_fitness_epistasis
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        get_survival_fitness = get_survival_fitness_hetdisadvantage
        short_survFitnessMethod = "Het"
    else
        println("ERROR--no survival fitness method chosen--should be either epistasis or hetdisadvantage")
    end

    # convert S_AM to pref_ratio (for use in math below)
    if S_AM == 1
        S_AM_for_math = 1 + 10^(-15)
    elseif S_AM == Inf
        S_AM_for_math = 10^(15)
    else
        S_AM_for_math = S_AM
    end
    pref_SD = sqrt(-1 / (2 * log(1 / S_AM_for_math)))  # width of female acceptance curve for male trait

    # generates the starting genotypes/locations
    # and calculates the growth rates based on individual resource use
    pd = PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp)

    extinction = false  # if extinction happens later this will be set true

    if do_plot
        plot_population(pd, functional_loci_range)
    end


    # loop throught the generations
    for generation in 1:max_generations
        genotypes_daughters_all = Matrix{Vector{Matrix{Int8}}}(undef, NUM_DEMES, NUM_DEMES)
        genotypes_sons_all = Matrix{Vector{Matrix{Int8}}}(undef, NUM_DEMES, NUM_DEMES)
        locations_daughters_all = Matrix{Vector{Location}}(undef, NUM_DEMES, NUM_DEMES)
        locations_sons_all = Matrix{Vector{Location}}(undef, NUM_DEMES, NUM_DEMES)
        mitochondria_daughters_all = Matrix{Vector{Int8}}(undef, NUM_DEMES, NUM_DEMES)
        mitochondria_sons_all = Matrix{Vector{Int8}}(undef, NUM_DEMES, NUM_DEMES)
        new_active_demes = CartesianIndex[]

        for i in 1:NUM_DEMES
            for j in 1:NUM_DEMES
                genotypes_daughters_all[i, j] = []
                genotypes_sons_all[i, j] = []
                locations_daughters_all[i, j] = []
                locations_sons_all[i, j] = []
                mitochondria_daughters_all[i, j] = []
                mitochondria_sons_all[i, j] = []
            end
        end

        if false
            active_demes = pd.active_demes
        else
            active_demes = eachindex(IndexCartesian(), pd.population)
        end

        for zone_index in active_demes
            # Prepare for mating and reproduction
            N_F = length(pd.population[zone_index].locations_F)
            father_deme = missing

            # loop through mothers, mating and reproducing
            for mother in 1:N_F
                # initialize tracking variables
                mate = false # becomes true when female is paired with male
                rejects = 0 # will track number of rejected males (in cases there is cost--which there isn't in main HZAM-Sym paper)
                father = [] # will contain the index of the male mate
                # make vector of indices of eligible males
                elig_M = Dict{CartesianIndex, Vector{Int64}}()

                unused_deme_indices = eachindex(IndexCartesian(), pd.population)

                neighbourhood_size = 0

                while mate==false
                    lower_left = max(zone_index - CartesianIndex(neighbourhood_size, neighbourhood_size), CartesianIndex(1, 1))
                    upper_right = min(zone_index + CartesianIndex(neighbourhood_size, neighbourhood_size), CartesianIndex(NUM_DEMES, NUM_DEMES))
                    neighbourhood = filter(e -> e ∈ unused_deme_indices, lower_left:upper_right) # remove used demes

                    for i in neighbourhood
                        elig_M[i] = 1:length(pd.population[i].locations_M)
                    end

                    if sum(map(length, values(elig_M)))==0
                        continue
                    end

                    unused_deme_indices = filter(e -> e ∉ neighbourhood, unused_deme_indices) # remove used demes from list

                    # determine male mate of female
                    while mate == false
                        # present female with random male (sympatric case) or closest male (spatial case), and remove him from list:
                        focal_male, elig_M, male_deme = choose_closest_male(pd.population, neighbourhood, elig_M, pd.population[zone_index].locations_F[mother])
                        # compare male trait with female's trait (preference), and determine
                        # whether she accepts; note that match_strength is determined by a
                        # Gaussian, with a maximum of 1 and minimum of zero.
                        match_strength = calc_match_strength(pd.population[zone_index].genotypes_F[mother], pd.population[male_deme].genotypes_M[focal_male], pref_SD, female_mating_trait_loci, male_mating_trait_loci)
                        if rand() < match_strength
                            # she accepts male, and mates
                            father = focal_male
                            father_deme = male_deme
                            mate = true
                        else
                            # she rejects male
                            rejects += 1
                            if sum(map(length, values(elig_M)))==0
                                break
                            end
                        end
                    end
                    neighbourhood_size += 1
                end

                # now draw the number of offspring from a poisson distribution with a mean of reproductive_fitness
                if !isempty(father)
                    # determine fitness cost due to mate search (number of rejected males)
                    search_fitness = (1 - per_reject_cost)^rejects    # (in most of HZAM-Sym paper, per_reject_cost = 0)

                    #combine for total fitness:   
                    reproductive_fitness = 2 * pd.growth_rates_F[zone_index][mother] * search_fitness  # the 2 is because only females, not males, produce offspring

                    offspring = rand(Poisson(reproductive_fitness))

                    # if offspring, generate their genotypes and sexes
                    if offspring >= 1
                        for kid in 1:offspring
                            kid_genotype = generate_offspring_genotype(pd.population[zone_index].genotypes_F[mother], pd.population[father_deme].genotypes_M[father])
                            kid_mitochondria = pd.population[zone_index].mitochondria_F[mother]
                            survival_fitness = get_survival_fitness(kid_genotype[:, hybrid_survival_loci], w_hyb)#######

                            if survival_fitness > rand()
                                # determine sex and location of kid

                                #genotype_average = sum(kid_genotype) / (2 * total_loci)
                                #hybrid = genotype_average != kid_mitochondria

                                new_location = Location(pd.population[zone_index].locations_F[mother], sigma_disp, geographic_limits)

                                zone = assign_zone(new_location)

                                if rand() > 0.5 # kid is daughter
                                    push!(genotypes_daughters_all[zone], kid_genotype)
                                    push!(mitochondria_daughters_all[zone], kid_mitochondria)
                                    push!(locations_daughters_all[zone], new_location)
                                else # kid is son
                                    push!(genotypes_sons_all[zone], kid_genotype)
                                    push!(mitochondria_sons_all[zone], kid_mitochondria)
                                    push!(locations_sons_all[zone], new_location)
                                end
                            end
                        end
                    end
                end # of loop through mothers
            end # of loop through the rows of zones
        end # of loop through the columns of zones

        #= check if either no daughters or no sons, and end the simulation if so
        if (size(genotypes_daughters, 3) == 0) || (size(genotypes_sons, 3) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end=#

        # assign surviving offspring to new adult population
        pd = update_population(pd,
            genotypes_daughters_all,
            genotypes_sons_all,
            mitochondria_daughters_all,
            mitochondria_sons_all,
            locations_daughters_all,
            locations_sons_all,
            new_active_demes,
            competition_trait_loci,
            K_total,
            sigma_comp,
            intrinsic_R,
            ecolDiff,
            false)

        # update the plot
        if (do_plot && (generation % plot_int == 0))
            update_plot(pd, generation, functional_loci_range)
        end

        push!(pop, sum([length(d.locations_F) for d in pd.population]) + sum([length(d.locations_M) for d in pd.population]))

    end # of loop through generations

    genotypes = vcat([[d.genotypes_F; d.genotypes_M] for d in pd.population]...)
    functional_HI_all_inds = calc_traits_additive(genotypes, functional_loci_range)
    HI_NL_all_inds = calc_traits_additive(genotypes, neutral_loci)

    println(string("average_pop: ", sum(pop) / length(pop)))

    return extinction, functional_HI_all_inds, HI_NL_all_inds
end

# This function calculates survival fitness of each individual according to heterozygosity.
function get_survival_fitness_hetdisadvantage(genotype::Array{Int8,2}, w_hyb::Real)::Vector{Float32}
    num_loci = size(genotype, 2)
    s_per_locus = 1 - w_hyb^(1 / num_loci)  # loss in fitness due to each heterozygous locus 
    num_hetloci = sum(genotype[1, :] .!= genotype[2, :])

    hetdisadvantage_fitness = (1 - s_per_locus)^num_hetloci
    return hetdisadvantage_fitness
end

# This function calculates survival fitness of each individual according to epistasis,
# with the beta parameter set to one as a default.
function get_survival_fitness_epistasis(genotype::Array{Int8,2}, w_hyb::Real, beta=1::Real)::Float32
    survival_HI = mean(genotype[:, :])
    epistasis_fitnesses = 1 - (1 - w_hyb) * (4 * survival_HI * (1 - survival_HI))^beta
    return epistasis_fitnesses
end

