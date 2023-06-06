include("population.jl")

using .Population

using Distributions # needed for "Poisson" function
using Statistics  # needed for "mean" function
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as csv
using LsqFit

global distances = Vector{Float64}() # tracks the distances between paired individuals


# This is the function to run a single HZAM simulation
function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3, neutral_loci=4:6,
    survival_fitness_method::String="epistasis", per_reject_cost=0,
    starting_pop_ratio=1.0, sympatry=false, geographic_limits::Vector{Float64}=[0.0, 1.0],
    starting_range_pop0=[0.0, 0.48], starting_range_pop1=[0.52, 1.0],
    sigma_disp=0.01, sigma_comp=0.01,
    do_plot=true, plot_int=10)

    # WHEN ECOLDIFF = 1, THEN SHOULD NOT HAVE THE 2 (1 INSTEAD)
    # WHEN ECOLDIFF = 0, THEN SHOULD HAVE THE 2
    # SO 2-ecolDiff ?


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
    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp)

    extinction = false  # if extinction happens later this will be set true

    # Pick function for choosing potential male mate (these two functions defined above) 
    if sympatry # if sympatry, then choose random male
        pick_potential_mate = choose_random_male
    else # if space matters, then choose closest male
        pick_potential_mate = choose_closest_male
    end
    if do_plot
        plot_population()
    end


    # loop throught the generations
    for generation in 1:max_generations

        # Prepare for mating and reproduction
        N_F, N_M = get_population()

        #println("generation: ", generation, "; individuals: ", N_F + N_M)

        #plot(spaced_locations, local_growth_rates_resourceA)
        #plot(spaced_locations, local_growth_rates_resourceB)  

        #plot(spaced_locations, local_growth_rates)  

        ### NEED TO CAREFULLY PROOF THE ABOVE, PARTICULARLY IF ECOLDIF > 0

        # make empty arrays for storing genotypes of daughters and sons
        genotypes_daughters = [zeros(Int8, 2, 3) for i in 1:0] # a little sketchy FIX LATER
        genotypes_sons = [zeros(Int8, 2, 3) for i in 1:0] # likewise

        # make empty arrays for storing locations of daughters and sons
        locations_daughters = Array{Float32,1}(undef, 0)
        locations_sons = Array{Float32,1}(undef, 0)

        # create structures for recording indices of mother and father (for error checking, and potentially for tracking genealogies)
        # first column for mother index (3rd dim of genotypes_F) and second column for father index (3rd dim of genotypes_M) 
        daughter_parent_IDs = Array{Int}(undef, 0, 2)
        son_parent_IDs = Array{Int}(undef, 0, 2)

        # flags for when it is necessary to expand the active region of the simulation
        expand_left = false
        expand_right = false

        # loop through mothers, mating and reproducing
        for mother in active_F
            # initialize tracking variables
            mate = false # becomes true when female is paired with male
            rejects = 0 # will track number of rejected males (in cases there is cost--which there isn't in main HZAM-Sym paper)
            father = [] # will contain the index of the male mate
            # make vector of indices of eligible males
            if (length(active_M)>0)
                elig_M = copy(active_M)
            else
                break
            end
            # determine male mate of female
            while mate == false
                # present female with random male (sympatric case) or closest male (spatial case), and remove him from list:
                focal_male, elig_M = pick_potential_mate(elig_M, mother)
                # compare male trait with female's trait (preference), and determine
                # whether she accepts; note that match_strength is determined by a
                # Gaussian, with a maximum of 1 and minimum of zero.
                match_strength = calc_match_strength(mother, focal_male, pref_SD)
                if rand() < match_strength
                    # she accepts male, and mates
                    father = focal_male
                    mate = true
                else
                    # she rejects male
                    rejects += 1
                    if length(elig_M) == 0
                        break
                    end
                end
            end

            # now draw the number of offspring from a poisson distribution with a mean of reproductive_fitness
            if !isempty(father)
                push!(distances, calc_distance(mother, father))
                # determine fitness cost due to mate search (number of rejected males)
                search_fitness = (1 - per_reject_cost)^rejects    # (in most of HZAM-Sym paper, per_reject_cost = 0)

                # determine fitness due to female use of available resources
                growth_rate_of_focal_female = calc_female_growth_rate(mother)

                #combine for total fitness:   
                reproductive_fitness = 2 * growth_rate_of_focal_female * search_fitness  # the 2 is because only females, not males, produce offspring

                offspring = rand(Poisson(reproductive_fitness))

                # if offspring, generate their genotypes and sexes
                if offspring >= 1
                    for kid in 1:offspring
                        kid_info = generate_offspring_genotype(mother, father, total_loci)
                        survival_fitness = get_survival_fitness(kid_info[:, hybrid_survival_loci], w_hyb)#######
                        if survival_fitness > rand()
                            # determine sex and location of kid
                            new_location = disperse_individual(mother, sigma_disp, geographic_limits)
                            genotype_sum = sum(kid_info)
                            if genotype_sum > 0 && new_location < left_boundary / 10
                                expand_left = true
                            elseif genotype_sum < total_loci * 2 && new_location > (right_boundary - 1) / 10
                                expand_right = true
                            end

                            if rand() > 0.5 # kid is daughter
                                push!(genotypes_daughters, kid_info)
                                locations_daughters = [locations_daughters; disperse_individual(mother, sigma_disp, geographic_limits)]
                            else # kid is son
                                push!(genotypes_sons, kid_info)
                                locations_sons = [locations_sons; disperse_individual(mother, sigma_disp, geographic_limits)]
                            end
                        end
                    end
                end
            end
        end # of loop through mothers

        # check if either no daughters or no sons, and end the simulation if so
        if (size(genotypes_daughters, 3) == 0) || (size(genotypes_sons, 3) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end

        # assign surviving offspring to new adult population
        update_population(genotypes_daughters, genotypes_sons, locations_daughters, locations_sons, expand_left, expand_right, generation)

        # check if there are no remaining females in the hybrid zone
        if (length(active_F) == 0)
            println("NO FEMALES REMAINING IN ACTIVE ZONE")
            break
        end

        # update the plot
        if (do_plot && (generation % plot_int == 0))
            update_plot(generation)
        end
    end # of loop through generations

    #return genotypes_F, locations_F, genotypes_M, locations_M, extinction
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

