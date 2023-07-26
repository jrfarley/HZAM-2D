using Distributions: Poisson # needed for "Poisson" functional_HI_all_inds


# This is the function to run a single HZAM simulation
function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3,
    survival_fitness_method::String="epistasis", per_reject_cost=0,
    geographic_limits::Vector{Location}=[Location(0.0f0, 0.0f0), Location(0.999f0, 0.999f0)],
    sigma_disp=0.05, sigma_comp=0.01,
    do_plot=true, plot_int=10)

    output_data = []
    hybrid_zone_positions = []
    functional_loci_range = union(female_mating_trait_loci, male_mating_trait_loci, competition_trait_loci, hybrid_survival_loci)# list of loci responsible for mating trait, competition trait, and hybrid survival

    #functional_loci_range = collect(1:total_loci)
    # get the chosen survival fitness function
    if survival_fitness_method == "epistasis"
        calc_survival_fitness = calc_survival_fitness_epistasis
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        calc_survival_fitness = calc_survival_fitness_hetdisadvantage
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
        # displays plot of individual locations and genotypes
        genotypes = [
            vcat([d.genotypes_F for d in pd.population]...)
            vcat([d.genotypes_M for d in pd.population]...)
        ]
        locations = [
            vcat([d.locations_F for d in pd.population]...)
            vcat([d.locations_M for d in pd.population]...)
        ]

        create_new_plot(
            calc_traits_additive(genotypes, functional_loci_range),
            locations
        )
    end


    # loop throught the generations
    for generation in 1:max_generations

        # sets up empty matrices to store the offspring genotypes, locations, and mitochondria
        genotypes_daughters_all = Matrix{Vector{Matrix{Int8}}}(undef, NUM_ZONES, NUM_ZONES)
        genotypes_sons_all = Matrix{Vector{Matrix{Int8}}}(undef, NUM_ZONES, NUM_ZONES)
        locations_daughters_all = Matrix{Vector{Location}}(undef, NUM_ZONES, NUM_ZONES)
        locations_sons_all = Matrix{Vector{Location}}(undef, NUM_ZONES, NUM_ZONES)
        mitochondria_daughters_all = Matrix{Vector{Int8}}(undef, NUM_ZONES, NUM_ZONES)
        mitochondria_sons_all = Matrix{Vector{Int8}}(undef, NUM_ZONES, NUM_ZONES)

        # fills all the matrices with empty vectors
        for i in 1:NUM_ZONES
            for j in 1:NUM_ZONES
                genotypes_daughters_all[i, j] = Matrix{Int8}[]
                genotypes_sons_all[i, j] = Matrix{Int8}[]
                locations_daughters_all[i, j] = Location[]
                locations_sons_all[i, j] = Location[]
                mitochondria_daughters_all[i, j] = Int8[]
                mitochondria_sons_all[i, j] = Int8[]
            end
        end

        for zone_index in eachindex(IndexCartesian(), pd.population)
            # Prepare for mating and reproduction
            N_F = length(pd.population[zone_index].locations_F)
            father_zone = missing

            # loop through mothers, mating and reproducing
            for mother in 1:N_F
                # initialize tracking variables
                mate = false # becomes true when female is paired with male
                rejects = 0 # will track number of rejected males (in cases there is cost--which there isn't in main HZAM-Sym paper)
                father = [] # will contain the index of the male mate


                # make dict where the keys are the zone indices and the values are vectors of indices of eligible males
                elig_M = Dict{CartesianIndex,Vector{Int64}}()

                for i in eachindex(IndexCartesian(), pd.population)
                    elig_M[i] = 1:length(pd.population[i].locations_M)
                end

                neighbourhood_size = 0.02f0 # how far away the simulation is checking for eligible mates (0 means only the current zone)


                while mate == false && neighbourhood_size < 0.1

                    # finds the coordinates for the zone that's the neighbourhood size away towards the bottom left
                    lower_left = max(assign_zone(Location(pd.population[zone_index].locations_F[mother].x - neighbourhood_size,
                            pd.population[zone_index].locations_F[mother].y - neighbourhood_size)),
                        CartesianIndex(1, 1))


                    # finds the coordinates for the zone that's the neighbourhood size away towards the top right
                    upper_right = min(assign_zone(Location(pd.population[zone_index].locations_F[mother].x + neighbourhood_size,
                            pd.population[zone_index].locations_F[mother].y + neighbourhood_size)),
                        CartesianIndex(NUM_ZONES, NUM_ZONES))

                    neighbourhood = filter(e -> length(elig_M[e]) > 0, lower_left:upper_right) # remove zones with no eligible males left

                    if length(neighbourhood) == 0 # if there are no eligible males remaining increase the neighbourhood size by 0.01 to look for males slightly further away
                        neighbourhood_size += 0.04f0
                        continue
                    end

                    # determine male mate of female
                    while mate == false
                        # present female with closest male, and remove him from list:
                        focal_male, male_zone = choose_closest_male(pd.population, neighbourhood, elig_M, pd.population[zone_index].locations_F[mother], neighbourhood_size)
                        if focal_male ≠ -1 # check if choose_closest_male returned a valid male index
                            # compare male trait with female's trait (preference), and determine
                            # whether she accepts; note that match_strength is determined by a
                            # Gaussian, with a maximum of 1 and minimum of zero.
                            match_strength = calc_match_strength(pd.population[zone_index].genotypes_F[mother], pd.population[male_zone].genotypes_M[focal_male], pref_SD, female_mating_trait_loci, male_mating_trait_loci)
                            if rand() < match_strength
                                # she accepts male, and mates
                                father = focal_male
                                father_zone = male_zone
                                mate = true
                            else
                                # she rejects male
                                rejects += 1
                                deleteat!(elig_M[male_zone], findall(x -> x == focal_male, elig_M[male_zone])) # remove rejected male from list of eligible males
                                if sum(map(length, values(elig_M))) == 0
                                    break
                                end
                            end
                        else
                            break # exit the loop since there are no eligible males left within the neighbourhood size distance
                        end
                    end
                    neighbourhood_size += 0.01f0 # increase the neighbourhood size to look for males slightly further away
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
                            kid_genotype = generate_offspring_genotype(pd.population[zone_index].genotypes_F[mother], pd.population[father_zone].genotypes_M[father])
                            kid_mitochondria = pd.population[zone_index].mitochondria_F[mother]
                            survival_fitness = calc_survival_fitness(kid_genotype, hybrid_survival_loci, w_hyb)

                            if survival_fitness > rand()
                                new_location = Location(pd.population[zone_index].locations_F[mother], sigma_disp) # generates offspring location based on dispersal distance, range, and location of mother

                                zone = assign_zone(new_location) # determines which zone the offspring dispersed into

                                #if zone in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
                                # add offspring data to the variables that track what the next generation's population will be
                                if rand() > 0.5 # kid is daughter
                                    push!(genotypes_daughters_all[zone], kid_genotype)
                                    push!(mitochondria_daughters_all[zone], kid_mitochondria)
                                    push!(locations_daughters_all[zone], new_location)
                                else # kid is son
                                    push!(genotypes_sons_all[zone], kid_genotype)
                                    push!(mitochondria_sons_all[zone], kid_mitochondria)
                                    push!(locations_sons_all[zone], new_location)
                                end
                                #end
                            end
                        end
                    end
                end
            end # of loop through mothers
        end # of loop through the zones

        # assign surviving offspring to new adult population
        pd = PopulationData(genotypes_daughters_all,
            genotypes_sons_all,
            mitochondria_daughters_all,
            mitochondria_sons_all,
            locations_daughters_all,
            locations_sons_all,
            competition_trait_loci,
            K_total,
            sigma_comp,
            intrinsic_R,
            ecolDiff)

        # updates the plot of locations and hybrid indices (plotting handled by the Plot_Data module)
        if (do_plot && (generation % plot_int == 0))
            genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
            locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]
            update_population_plot(
                calc_traits_additive(genotypes, functional_loci_range),
                locations,
                generation
            )
        end

        if generation > max_generations - 20

            neutral_loci = setdiff(collect(1:total_loci), functional_loci_range)
            genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
            locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]
            mitochondria = [vcat([d.mitochondria_F for d in pd.population]...); vcat([d.mitochondria_M for d in pd.population]...)]

            hybrid_indices_functional = calc_traits_additive(genotypes, functional_loci_range)

            gene_flows = calc_all_gene_flow(genotypes, hybrid_indices_functional, female_mating_trait_loci, male_mating_trait_loci, competition_trait_loci, hybrid_survival_loci, neutral_loci)
            cline_widths, cline_positions = calc_all_cline_widths(genotypes, locations, female_mating_trait_loci, male_mating_trait_loci, competition_trait_loci, hybrid_survival_loci, total_loci)
            average_linkage_diseq = calc_average_linkage_diseq(genotypes, female_mating_trait_loci, male_mating_trait_loci, competition_trait_loci, hybrid_survival_loci, neutral_loci)

            push!(output_data, calc_output_data(hybrid_indices_functional, locations, sigma_disp, gene_flows, genotypes, cline_widths, cline_positions, mitochondria, average_linkage_diseq))
        end

    end # of loop through generations

    return average_output_data(output_data), pd
end

