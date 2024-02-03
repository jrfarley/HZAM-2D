using Distributions: Poisson # needed for "Poisson" functional_HI_all_inds

"""
    run_one_HZAM_sim(
        w_hyb::Real, 
        S_AM::Real,
        intrinsic_R::Real; 
        <keyword arguments>
    )

Run a single HZAM simulation.

# Arguments
- `w_hyb::Real`: the hybrid fitness.
- `S_AM::Real`: the strength of assortative mating.
- `intrinsic_R::Real`: the intrinsic growth rate.
- `K_total::Integer=20000`: the carrying capacity of the environment.
- `max_generations::Integer=1000`: the number of generations that the simulation will run for.
- `total_loci::Integer=6`: the total number of loci in the genome.
- `female_mating_trait_loci=1:3`: the loci specifying the female's mate preference.
- `male_mating_trait_loci=1:3`: the loci specifying the male's mating trait.
- `hybrid_survival_loci=1:3`: the loci specifying the probability of survival to adulthood.
- `survival_fitness_method:String="epistasis"`: the method used to calculate the probability of survival to adulthood.
- `per_reject_cost=0`: the fitness loss of female per male rejected (due to search time, etc.). Can take values of 0 to 1.
- `sigma_disp::Float32=0.05f0`: the standard deviation of the normal distribution determining how far offspring will disperse from their mothers.
- `sigma_comp::Float32=0.01f0`: the standard deviation for the normal curve used in calculating local density.
- `do_plot=true`: the program will display a plot of the locations and hybrid indices of every individual while the simulation is running if this is true.
- `plot_int=10`: the interval (measured in generations) between updating the plot.
- `gene_plot=false`: if true, generates phenotype plots.
- `save_plot=false`: if true, saves each plot to a PNG file.
- `track_population_data=false`: if true, stores the population size, hybridity, overlap, and cline width for each generation
"""
function run_one_HZAM_sim(w_hyb::Real, S_AM::Real, intrinsic_R::Real;
    # the semicolon makes the following optional keyword arguments  
    K_total::Integer=40000, max_generations::Integer=1000,
    total_loci::Integer=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    hybrid_survival_loci=1:3, survival_fitness_method::String="epistasis",
    per_reject_cost=0, sigma_disp=0.03f0,
    sigma_comp=0.01f0, do_plot=true, plot_int=10, gene_plot=false, save_plot=false,
    track_population_data=false)

    population_tracking_data = DataAnalysis.PopulationTrackingData[]

    overall_loci_range = collect(1:total_loci)

    # list of loci specifying the phenotype
    functional_loci_range = union(
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci
    )

    # list of loci that do not affect the phenotype
    neutral_loci_range = setdiff(overall_loci_range, functional_loci_range)

    # put all the loci together in an easy to access data format
    loci = @NamedTuple{
        overall,
        functional,
        neutral,
        female_mating_trait,
        male_mating_trait,
        hybrid_survival
    }((
        overall_loci_range,
        functional_loci_range,
        neutral_loci_range,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci
    ))

    # get the chosen survival fitness function
    if survival_fitness_method == "epistasis"
        calc_survival_fitness = Population.calc_survival_fitness_epistasis
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        calc_survival_fitness = Population.calc_survival_fitness_hetdisadvantage
        short_survFitnessMethod = "Het"
    else
        println("ERROR--no survival fitness method chosen--should be either epistasis or 
        hetdisadvantage")
    end


    # width of female acceptance curve for male trait
    pref_SD = 1

    # convert S_AM to pref_SD (for use in math below)
    if S_AM == 1
        pref_SD = Inf
    elseif S_AM == Inf
        pref_SD = 10^(-15)
    else
        pref_SD = sqrt(-1 / (2 * log(1 / S_AM)))
    end


    #=
    Generate the starting genotypes/locations and calculate the growth rates 
    for all individuals in the simulation.
    =#
    pd = Population.PopulationData(
        K_total,
        total_loci,
        intrinsic_R,
        sigma_comp
    )


    if do_plot
        # display plot of individual locations and genotypes
        genotypes = [
            vcat([d.genotypes_F for d in pd.population]...)
            vcat([d.genotypes_M for d in pd.population]...)
        ]
        x_locations = [
            vcat([d.x_locations_F for d in pd.population]...)
            vcat([d.x_locations_M for d in pd.population]...)
        ]
        y_locations = [
            vcat([d.y_locations_F for d in pd.population]...)
            vcat([d.y_locations_M for d in pd.population]...)
        ]

        if gene_plot
            PlotData.create_gene_plot(genotypes, loci, save_plot)
        else
            PlotData.create_population_plot(
                DataAnalysis.calc_traits_additive(genotypes, functional_loci_range),
                x_locations,
                y_locations,
                save_plot
            )
        end
    end

    zone_cartesian_indices = CartesianIndex(1, 1):CartesianIndex(NUM_ZONES, NUM_ZONES)

    # loop through the generations
    for generation in 1:max_generations

        # set up empty matrices to store the offspring data
        genotypes_daughters_all = [Matrix{Int8}[] for i in zone_cartesian_indices]
        genotypes_sons_all = [Matrix{Int8}[] for i in zone_cartesian_indices]
        x_locations_daughters_all = [Float32[] for i in zone_cartesian_indices]
        y_locations_daughters_all = [Float32[] for i in zone_cartesian_indices]
        x_locations_sons_all = [Float32[] for i in zone_cartesian_indices]
        y_locations_sons_all = [Float32[] for i in zone_cartesian_indices]
        #offspring_per_parent_F = Matrix{Vector{Integer}}(undef, NUM_ZONES, NUM_ZONES)
        #offspring_per_parent_M = Matrix{Vector{Integer}}(undef, NUM_ZONES, NUM_ZONES)

        for zone_index in zone_cartesian_indices
            # prepare for mating and reproduction

            zone = pd.population[zone_index]
            father_zone = missing

            # loop through mothers, mating and reproducing
            for mother in eachindex(zone.x_locations_F)
                # initialize tracking variables
                mate = false # becomes true when female is paired with male
                rejects = 0 # will track number of rejected males 
                father = missing # will contain the index of the male mate


                #=
                Construct a dictionary where the keys are the zone indices and the values 
                are vectors of the indices of the eligible males in each zone.
                =#
                elig_M = Dict{CartesianIndex,Vector{Int64}}()

                for i in eachindex(IndexCartesian(), pd.population)
                    elig_M[i] = 1:length(pd.population[i].x_locations_M)
                end

                # set a distance cutoff on the search for a mate
                neighbourhood_size::Float32 = sigma_comp


                while mate == false && neighbourhood_size <= Float32(1/(2*Population.NUM_ZONES))

                    #=
                    Find the index for the zone that's the neighbourhood size away 
                    towards the bottom left.
                    =#
                    lower_left = max(
                        Population.assign_zone(
                            zone.x_locations_F[mother] - neighbourhood_size,
                            zone.y_locations_F[mother] - neighbourhood_size
                        ),
                        CartesianIndex(1, 1)
                    )


                    #=
                    Find the index for the zone that's the neighbourhood size away 
                    towards the top right.
                    =#
                    upper_right = min(
                        Population.assign_zone(
                            zone.x_locations_F[mother] + neighbourhood_size,
                            zone.y_locations_F[mother] + neighbourhood_size
                        ),
                        CartesianIndex(NUM_ZONES, NUM_ZONES)
                    )

                    # collect the indices of the nearby zones
                    neighbourhood = [lower_left:upper_right...]

                    # remove zones with no eligible males left
                    filter!(i -> length(elig_M[i]) > 0, neighbourhood)

                    #=
                    If none of the zones in the neighbourhood have eligible males remaining, 
                    increase the neighbourhood size by 0.01 to look for males further away.
                    =#
                    if length(neighbourhood) == 0
                        neighbourhood_size += sigma_comp
                        continue
                    end

                    # determine male mate of female
                    while mate == false
                        # present female with closest male, and remove him from list:
                        focal_male, male_zone = Mating.choose_closest_male(
                            pd.population,
                            neighbourhood,
                            elig_M,
                            zone.x_locations_F[mother],
                            zone.y_locations_F[mother],
                            neighbourhood_size
                        )

                        # check if choose_closest_male returned a valid male index
                        if focal_male ≠ -1
                            #=
                            Compare the male trait with female's trait (preference), and 
                            determine whether she accepts; note that match_strength is 
                            determined by a Gaussian, with a maximum of 1 and minimum of 
                            zero.
                            =#
                            match_strength = Mating.calc_match_strength(
                                zone.genotypes_F[mother],
                                pd.population[male_zone].genotypes_M[focal_male],
                                pref_SD,
                                female_mating_trait_loci,
                                male_mating_trait_loci
                            )
                            if rand() < match_strength
                                # she accepts male, and mates
                                father = focal_male
                                father_zone = male_zone
                                mate = true
                            else
                                # she rejects male
                                rejects += 1
                                # remove rejected male from list of eligible males
                                deleteat!(
                                    elig_M[male_zone],
                                    findall(x -> x == focal_male, elig_M[male_zone])
                                )

                                # exit the loop if there are no remaining eligible males
                                if sum(map(length, values(elig_M))) == 0
                                    break
                                end
                            end
                        else
                            # there are no remaining eligible males in the neighbourhood
                            break
                        end
                    end
                    # increase the neighbourhood size to look for males further away
                    neighbourhood_size += sigma_comp
                end

                #=
                Now draw the number of offspring from a poisson distribution with a mean 
                value of the reproductive_fitness.
                =#
                if !ismissing(father)
                    # determine fitness cost due to mate search (number of rejected males)
                    search_fitness = (1 - per_reject_cost)^rejects # typically equals 1

                    #combine for total fitness:   
                    # the 2 is because only females, not males, produce offspring
                    reproductive_fitness = 2 *
                                           pd.growth_rates_F[zone_index][mother] *
                                           search_fitness

                    offspring = rand(Poisson(reproductive_fitness))

                    # if there are offspring, generate their genotypes and sexes
                    if offspring >= 1
                        for kid in 1:offspring
                            kid_genotype = Mating.generate_offspring_genotype(
                                zone.genotypes_F[mother],
                                pd.population[father_zone].genotypes_M[father]
                            )
                            survival_fitness = calc_survival_fitness(
                                kid_genotype,
                                hybrid_survival_loci,
                                w_hyb
                            )

                            if survival_fitness > rand()
                                #offspring_per_parent_F[zone_index][mother] += 1
                                #offspring_per_parent_M[father_zone][father] += 1

                                # generate offspring location
                                x, y = Population.new_location(
                                    zone.x_locations_F[mother],
                                    zone.y_locations_F[mother],
                                    sigma_disp
                                )

                                # determine which zone the offspring dispersed into
                                kid_zone = Population.assign_zone(x, y)

                                # update variables for tracking offspring data
                                if rand() > 0.5 # kid is daughter
                                    push!(genotypes_daughters_all[kid_zone], kid_genotype)
                                    push!(x_locations_daughters_all[kid_zone], x)
                                    push!(y_locations_daughters_all[kid_zone], y)
                                else # kid is son
                                    push!(genotypes_sons_all[kid_zone], kid_genotype)
                                    push!(x_locations_sons_all[kid_zone], x)
                                    push!(y_locations_sons_all[kid_zone], y)
                                end
                                #end
                            end
                        end
                    end
                end
            end # of loop through mothers
        end # of loop through the zones

        # assign surviving offspring to new adult population
        pd = Population.PopulationData(
            genotypes_daughters_all,
            genotypes_sons_all,
            x_locations_daughters_all,
            y_locations_daughters_all,
            x_locations_sons_all,
            y_locations_sons_all,
            K_total,
            sigma_comp,
            intrinsic_R
        )

        # update the plot of locations and hybrid indices
        if (do_plot && (generation % plot_int == 0))
            genotypes = [
                vcat([d.genotypes_F for d in pd.population]...)
                vcat([d.genotypes_M for d in pd.population]...)
            ]
            x_locations = [
                vcat([d.x_locations_F for d in pd.population]...)
                vcat([d.x_locations_M for d in pd.population]...)
            ]
            y_locations = [
                vcat([d.y_locations_F for d in pd.population]...)
                vcat([d.y_locations_M for d in pd.population]...)
            ]


            if gene_plot
                PlotData.update_gene_plot(genotypes, loci, generation, save_plot)
            else
                PlotData.update_population_plot(
                    DataAnalysis.calc_traits_additive(genotypes, functional_loci_range),
                    x_locations,
                    y_locations,
                    generation,
                    save_plot
                )
            end
        end

        if generation % 100 == 0 
            println("generation: $generation")
        end

        if track_population_data
            genotypes = [
                vcat([d.genotypes_F for d in pd.population]...)
                vcat([d.genotypes_M for d in pd.population]...)
            ]

            x_locations = [
                vcat([d.x_locations_F for d in pd.population]...)
                vcat([d.x_locations_M for d in pd.population]...)
            ]
            y_locations = [
                vcat([d.y_locations_F for d in pd.population]...)
                vcat([d.y_locations_M for d in pd.population]...)
            ]

            overlap = Population.calc_species_overlap(pd.population, 0.03, sigma_comp, male_mating_trait_loci)

            push!(
                population_tracking_data,
                DataAnalysis.PopulationTrackingData(genotypes, locations, male_mating_trait_loci, overlap)
            )
        end
    end # of loop through generations

    parameters = DataAnalysis.SimParams(
        intrinsic_R,
        w_hyb,
        S_AM,
        K_total,
        max_generations,
        sigma_disp,
        total_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        per_reject_cost
    )

    genotypes = [
        vcat([d.genotypes_F for d in pd.population]...)
        vcat([d.genotypes_M for d in pd.population]...)
    ]
    x_locations = [
        vcat([d.x_locations_F for d in pd.population]...)
        vcat([d.x_locations_M for d in pd.population]...)
    ]
    y_locations = [
        vcat([d.y_locations_F for d in pd.population]...)
        vcat([d.y_locations_M for d in pd.population]...)
    ]

    overlap = Population.calc_species_overlap(pd.population, 0.03, sigma_comp, male_mating_trait_loci)

    output = DataAnalysis.OutputData(
        genotypes,
        x_locations,
        y_locations,
        male_mating_trait_loci,
        parameters,
        population_tracking_data,
        overlap,
        pd
    )
    fig = PlotData.create_population_plot(
        DataAnalysis.calc_traits_additive(genotypes, functional_loci_range),
        x_locations,
        y_locations,
        false
    )

    return output
end # of module

