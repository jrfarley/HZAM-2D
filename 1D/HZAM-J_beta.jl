#= HZAM-J_beta.jl
Hybrid Zone with Assortative Mating, Julia implementation,
by Darren Irwin.
Starting this file on 30Jan2022 to develop an integrated HZAM model that 
combines features of HZAM_Sym_Julia_V1.jl and HZAM_release_2.0.R into 

If you use this code, please cite the paper below, and possibly my GitHub or Dryad repository 
where you obtained this code:

Irwin, D., and Schluter, D. Hybridization and the coexistence of species. bioRxiv 2021.04.04.438369; doi: https://doi.org/10.1101/2021.04.04.438369 
=#

#= You need to make sure these packages below are added to your Julia environment.
To add e.g. the package Distributions, 
type "]" to get the "pkg>" prompt, then type e.g. "add Distributions";
or run the commands below: 
=#
# import Pkg; 
# Pkg.add("Distributions") 
# Pkg.add("Statistics") 
# Pkg.add("JLD2")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("LsqFit") 
### Pkg.add("Plots")
# Pkg.add("CategoricalArrays")
# Pkg.add("Colors")
# Pkg.add("ColorSchemes")
### Pkg.add("Makie") 
### Pkg.add("CairoMakie")
# Pkg.add("GLMakie")


using Distributions # needed for "Poisson" function
using Statistics  # needed for "mean" function

# for plotting:
using Plots

# to start Julia with multiple threads, type in terminal e.g.:
# julia --threads 4
# To check, type in Julia: Threads.nthreads()
# 4

# set up functions (should not define functions repeatedly in loop, as causes re-compilation, slows things)

# This function sets up the genotypes of the starting population
# in a 3D array, where:
# rows (D1) are alleles (row 1 from mother, row 2 from father),
# columns (D2) are loci, 
# pages (D3) are individuals.  
function generate_genotype_array(N_pop0::Integer, N_pop1::Integer, loci::Integer)::Array{Int8,3}
    total_N = N_pop0 + N_pop1
    genotypes = Array{Int8,3}(undef, 2, loci, total_N) # The "Int8" is the type (8-bit integer), and "undef" means an unitialized array, so values are meaningless
    genotypes[:, :, 1:N_pop0] .= 0  # assigns genotypes of pop01
    genotypes[:, :, (N_pop0+1):total_N] .= 1  # assigns genotypes of pop1
    return genotypes
end

# This function calculates each mean values of the genotypes passed to it (for each individual).
# Used to determine trait values in an additive way.
# Only those loci that are additive trait loci should be passed to this function.
function calc_traits_additive(genotypes::Array{Int8,3})::Vector{Float32}
    return [mean(genotypes[:, :, i]) for i in 1:size(genotypes, 3)]
end

# This function calculates survival fitness of each individual according to epistasis,
# with the beta parameter set to one as a default.
function get_survival_fitnesses_epistasis(genotypes::Array{Int8,3}, w_hyb::Real, beta=1::Real)::Vector{Float32}
    survival_HI = calc_traits_additive(genotypes)
    epistasis_fitnesses = 1 .- (1 - w_hyb) .* (4 .* survival_HI .* (1 .- survival_HI)) .^ beta
    return epistasis_fitnesses
end

# This function calculates survival fitness of each individual according to heterozygosity.
function get_survival_fitnesses_hetdisadvantage(genotypes::Array{Int8,3}, w_hyb::Real)::Vector{Float32}
    N = size(genotypes, 3)
    num_loci = size(genotypes, 2)
    s_per_locus = 1 - w_hyb^(1 / num_loci)  # loss in fitness due to each heterozygous locus 
    num_hetloci = Vector{Integer}(undef, N)
    for ind in 1:N  # count number of het loci per individual
        num_hetloci[ind] = sum(genotypes[1, :, ind] .!= genotypes[2, :, ind])
    end
    hetdisadvantage_fitnesses = (1 - s_per_locus) .^ num_hetloci
    return hetdisadvantage_fitnesses
end

# This function determines breeding location of one individual,
# based on a normal distribution with width sigma_disp,
# centred on birth location. Constrained to be within range. 
# geographic_limits should be a vector with two numbers.
function disperse_individual(start_location::Real, sigma_disp::Real, geographic_limits::Vector{Float64})::Float32
    while true
        new_location = start_location + (sigma_disp * randn())
        if (new_location >= geographic_limits[1]) & (new_location <= geographic_limits[2]) # checks if location is in range
            return new_location
        end
    end
end

# These next two functions define the way potential male mates are chosen.
# elig_M is a vector of indices of possible male mates.
# When in sympatric model, choose random male (note the second and third arguments not used but allows function to be called in same way as in spatial model):
function choose_random_male(elig_M::Vector{UInt32}, locations_M::Vector{Float32}, focal_location::Real)
    focal_male = splice!(elig_M, rand(eachindex(elig_M))) # this gets the index of a random male, and removes that male from the list in elig_M
    return focal_male, elig_M
end
# When in spatial model, choose closest male:
function choose_closest_male(elig_M::Vector{UInt32}, locations_M::Vector{Float32}, focal_location::Real)
    focal_male = splice!(elig_M, argmin(abs.(locations_M[elig_M] .- focal_location))) # this gets the index of a closest male, and removes that male from the list in elig_M
    return focal_male, elig_M
end
# This is the function to run a single HZAM simulation
function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int=1000, max_generations::Int=1000,
    total_loci::Int=6, female_mating_trait_loci=1:3, male_mating_trait_loci=1:3,
    competition_trait_loci=1:3, hybrid_survival_loci=1:3, neutral_loci=4:6,
    survival_fitness_method::String="epistasis", per_reject_cost=0,
    starting_pop_ratio=1.0, sympatry=false, geographic_limits::Vector{Float64}=[0.0, 1.0],
    starting_range_pop0=[0.0, 0.48], starting_range_pop1=[0.52, 1.0],
    sigma_disp=0.02, sigma_comp=0.01,
    do_plot=true, plot_int=10)

    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    competAbility_useResourceA_pop0 = (1 + ecolDiff) / 2    # equals 1 when ecolDiff = 1   
    competAbility_useResourceB_pop0 = 1 - competAbility_useResourceA_pop0
    competAbility_useResourceA_pop1 = (1 - ecolDiff) / 2   # equals 0 when ecolDiff = 1
    competAbility_useResourceB_pop1 = 1 - competAbility_useResourceA_pop1
    # set up carying capacities on each resource, and starting pop sizes of each species
    K_A = K_total / 2  # EVEN NUMBER; carrying capacity (on resource alpha) of entire range (for two sexes combined), regardless of species 
    K_B = K_total / 2   # EVEN NUMBER; carrying capacity (on resource beta) of entire range (for two sexes combined), regardless of species

    position = []

    if sympatry  # if sympatry is true, then set both pop sizes according to full range 
        starting_range_pop0 = [0.0, 1.0]
        starting_range_pop1 = [0.0, 1.0]
        pop0_starting_N = K_A   # starting N of species 0
        pop0_starting_N_half = Int(pop0_starting_N / 2)  # The "Int" is to ensure no decimal
        pop1_starting_N = Int(round(starting_pop_ratio * K_B))   # starting N of species 1, which can be lower if starting_pop_ratio is below 1)
        pop1_starting_N_half = Int(pop1_starting_N / 2)
    else  # sympatry is false, then set pop sizes assuming no range overlap--this is why the "2" is in the formulae below  
        pop0_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop0) + (K_B * competAbility_useResourceB_pop0)) * (starting_range_pop0[2] - starting_range_pop0[1]) / (geographic_limits[2] - geographic_limits[1]))
        pop0_starting_N_half = Int(pop0_starting_N / 2)  # starting number for each sex
        pop1_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop1) + (K_B * competAbility_useResourceB_pop1)) * (starting_range_pop1[2] - starting_range_pop1[1]) / (geographic_limits[2] - geographic_limits[1]))
        pop1_starting_N_half = Int(pop1_starting_N / 2)
    end

    # WHEN ECOLDIFF = 1, THEN SHOULD NOT HAVE THE 2 (1 INSTEAD)
    # WHEN ECOLDIFF = 0, THEN SHOULD HAVE THE 2
    # SO 2-ecolDiff ?

    # get the chosen survival fitness function
    if survival_fitness_method == "epistasis"
        get_survival_fitnesses = get_survival_fitnesses_epistasis
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        get_survival_fitnesses = get_survival_fitnesses_hetdisadvantage
        short_survFitnessMethod = "Het"
    else
        println("ERROR--no survival fitness method chosen--should be either epistasis or hetdisadvantage")
    end

    # add up and check the total loci
    total_functional_loci = max(maximum(female_mating_trait_loci), maximum(male_mating_trait_loci), maximum(competition_trait_loci), maximum(hybrid_survival_loci))
    functional_loci_range = 1:total_functional_loci
    num_neutral_loci = length(Vector(neutral_loci))
    if total_functional_loci + num_neutral_loci ≠ total_loci
        println("#### WARNING: Please examine your loci numbers and indices, as they don't all match up ####")
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

    # Generate genotype array for population of females:
    # this is a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
    # columns (D2) are loci, and pages (D3) are individuals
    genotypes_F = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    genotypes_M = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    # functional loci are first, followed by neutral loci

    # Generate breeding locations of individuals (vector in same order of individuals as D3 of the genotype array)
    locations_F_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_F_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_F = [locations_F_pop0; locations_F_pop1]
    locations_M_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_M_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_M = [locations_M_pop0; locations_M_pop1]

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])
    ind_locations_if_even_at_K = range(geographic_limits[1], geographic_limits[2], length=K_total)
    function get_density_if_even_at_K(focal_location) # this function calculates local density according to a normal curve
        sum(exp.(-((ind_locations_if_even_at_K .- focal_location) .^ 2) ./ (2 * (sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end
    ideal_densities_at_spaced_locations = map(get_density_if_even_at_K, spaced_locations) # this applies the above function to each geographic location
    # assume both resources have same constant density across range
    ideal_densities_at_spaced_locations_resourceA = ideal_densities_at_spaced_locations ./ 2
    ideal_densities_at_spaced_locations_resourceB = ideal_densities_at_spaced_locations ./ 2

    #plot(spaced_locations, ideal_densities_at_spaced_locations_resourceA)

    extinction = false  # if extinction happens later this will be set true

    # Pick function for choosing potential male mate (these two functions defined above) 
    if sympatry # if sympatry, then choose random male
        pick_potential_mate = choose_random_male
    else # if space matters, then choose closest male
        pick_potential_mate = choose_closest_male
    end

    # loop throught the generations
    for generation in 1:max_generations

        # Prepare for mating and reproduction
        N_F = size(genotypes_F, 3)
        N_M = size(genotypes_M, 3)

        println("generation: ", generation, "; individuals: ", N_F + N_M)

        # calculate mating trait values (T) from genotypes
        female_mating_traits = calc_traits_additive(genotypes_F[:, female_mating_trait_loci, :])
        male_mating_traits = calc_traits_additive(genotypes_M[:, male_mating_trait_loci, :])

        # calculate ecological competition trait values from genotypes
        competition_traits_F = calc_traits_additive(genotypes_F[:, competition_trait_loci, :])
        competition_traits_M = calc_traits_additive(genotypes_M[:, competition_trait_loci, :])

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA_F = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_F) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
        ind_useResourceB_F = competAbility_useResourceB_pop0 .+ (competition_traits_F .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))
        ind_useResourceA_M = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_M) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
        ind_useResourceB_M = competAbility_useResourceB_pop0 .+ (competition_traits_M .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))

        if sympatry
            # sum up the global resource use over all individuals
            total_useResourceA = sum(ind_useResourceA_F) + sum(ind_useResourceA_M)
            total_useResourceB = sum(ind_useResourceB_F) + sum(ind_useResourceB_M)
            # calculate global growth rates due to each resource (according to discrete time logistic growth equation)
            growth_rate_resourceA = (intrinsic_R * K_A) / (K_A + ((total_useResourceA) * (intrinsic_R - 1)))
            growth_rate_resourceB = (intrinsic_R * K_B) / (K_B + ((total_useResourceB) * (intrinsic_R - 1)))

        else  # in spatial model, calculate growth rates based on local resource use
            # determine local resource use for each location across range
            ind_locations_real = [locations_F; locations_M]
            ind_useResourceA_all = [ind_useResourceA_F; ind_useResourceA_M]
            ind_useResourceB_all = [ind_useResourceB_F; ind_useResourceB_M]
            function get_useResourceA_density_real(focal_location) # this function calculates local density according to a normal curve, weighted by individual resource use
                sum(ind_useResourceA_all .* exp.(-((ind_locations_real .- focal_location) .^ 2) ./ (2 * (sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
            end
            real_densities_at_spaced_locations_resourceA = map(get_useResourceA_density_real, spaced_locations) # this applies the above function to each geographic location
            function get_useResourceB_density_real(focal_location) # do the same for resource B
                sum(ind_useResourceB_all .* exp.(-((ind_locations_real .- focal_location) .^ 2) ./ (2 * (sigma_comp^2))))
            end
            real_densities_at_spaced_locations_resourceB = map(get_useResourceB_density_real, spaced_locations) # this applies the above function to each geographic location 
            # calculate local growth rates due to each resource (according to discrete time logistic growth equation)
            local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_spaced_locations_resourceA ./ (ideal_densities_at_spaced_locations_resourceA .+ ((real_densities_at_spaced_locations_resourceA) .* (intrinsic_R - 1)))
            local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_spaced_locations_resourceB ./ (ideal_densities_at_spaced_locations_resourceB .+ ((real_densities_at_spaced_locations_resourceB) .* (intrinsic_R - 1)))
        end

        #plot(spaced_locations, local_growth_rates_resourceA)
        #plot(spaced_locations, local_growth_rates_resourceB)  

        #plot(spaced_locations, local_growth_rates)  

        ### NEED TO CAREFULLY PROOF THE ABOVE, PARTICULARLY IF ECOLDIF > 0


        # Set up structure to record number of matings per male (and female, although almost always 1 for females), 
        # to determine sexual selection due to HI class:
        matings_per_male = zeros(Int8, N_M)
        matings_per_female = zeros(Int8, N_F)

        # create empty data structures for keeping track of numbers offspring of parents
        daughters_per_mother = zeros(Int16, N_F)
        sons_per_mother = zeros(Int16, N_F)
        daughters_per_father = zeros(Int16, N_M)
        sons_per_father = zeros(Int16, N_M)

        # make empty arrays for storing genotypes of daughters and sons
        genotypes_daughters = Array{Int8,3}(undef, 2, total_loci, 0)
        genotypes_sons = Array{Int8,3}(undef, 2, total_loci, 0)

        # make empty arrays for storing locations of daughters and sons
        locations_daughters = Array{Float32,1}(undef, 0)
        locations_sons = Array{Float32,1}(undef, 0)

        # create structures for recording indices of mother and father (for error checking, and potentially for tracking genealogies)
        # first column for mother index (3rd dim of genotypes_F) and second column for father index (3rd dim of genotypes_M) 
        daughter_parent_IDs = Array{Int}(undef, 0, 2)
        son_parent_IDs = Array{Int}(undef, 0, 2)

        # loop through mothers, mating and reproducing
        for mother in 1:N_F
            # initialize tracking variables
            mate = false # becomes true when female is paired with male
            rejects = 0 # will track number of rejected males (in cases there is cost--which there isn't in main HZAM-Sym paper)
            father = [] # will contain the index of the male mate
            # make vector of indices of eligible males
            elig_M = Vector{UInt32}(1:N_M)  # this integer type allows up to more than 4 billion values 
            # determine male mate of female
            while mate == false
                # present female with random male (sympatric case) or closest male (spatial case), and remove him from list:
                focal_male, elig_M = pick_potential_mate(elig_M, locations_M, locations_F[mother])
                # compare male trait with female's trait (preference), and determine
                # whether she accepts; note that match_strength is determined by a
                # Gaussian, with a maximum of 1 and minimum of zero.
                match_strength = (exp(1)^((-(male_mating_traits[focal_male] - female_mating_traits[mother])^2) / (2 * (pref_SD^2))))
                if rand() < match_strength
                    # she accepts male, and mates
                    father = focal_male
                    matings_per_male[focal_male] += 1 # this adds 1 to the matings for that male
                    matings_per_female[mother] += 1
                    mate = true
                else
                    # she rejects male
                    rejects += 1
                    if length(elig_M) == 0
                        break
                    end
                end
            end
            # determine fitness cost due to mate search (number of rejected males)
            search_fitness = (1 - per_reject_cost)^rejects    # (in most of HZAM-Sym paper, per_reject_cost = 0)

            # determine fitness due to female use of available resources
            if sympatry
                growth_rate_of_focal_female = (ind_useResourceA_F[mother] * growth_rate_resourceA) + (ind_useResourceB_F[mother] * growth_rate_resourceB)
            else
                location_index = argmin(abs.(spaced_locations .- locations_F[mother]))
                local_growth_A = local_growth_rates_resourceA[location_index]
                local_growth_B = local_growth_rates_resourceB[location_index]
                growth_rate_of_focal_female = (ind_useResourceA_F[mother] * local_growth_A) + (ind_useResourceB_F[mother] * local_growth_B)
            end

            #combine for total fitness:   
            reproductive_fitness = 2 * growth_rate_of_focal_female * search_fitness  # the 2 is because only females, not males, produce offspring
            # now draw the number of offspring from a poisson distribution with a mean of reproductive_fitness
            if isempty(father)
                offspring = 0  # if no mate (because all males were rejected), then no offspring
            else
                offspring = rand(Poisson(reproductive_fitness))
            end
            # if offspring, generate their genotypes and sexes
            if offspring >= 1
                for kid in 1:offspring
                    kid_info = Array{Int8,2}(undef, 2, total_loci) # place to store genotype of one offspring
                    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
                    for locus in 1:total_loci
                        kid_info[1, locus] = genotypes_F[rand([1 2]), locus, mother]  # for this locus, pick a random allele from the mother
                        kid_info[2, locus] = genotypes_M[rand([1 2]), locus, father]  # and from the father
                    end
                    # determine sex and location of kid
                    if rand() > 0.5 # kid is daughter
                        genotypes_daughters = cat(genotypes_daughters, kid_info, dims=3)
                        daughter_parent_IDs = cat(daughter_parent_IDs, [mother father], dims=1)
                        daughters_per_mother[mother] += 1
                        daughters_per_father[father] += 1
                        locations_daughters = [locations_daughters; disperse_individual(locations_F[mother], sigma_disp, geographic_limits)]
                    else # kid is son
                        genotypes_sons = cat(genotypes_sons, kid_info, dims=3)
                        son_parent_IDs = cat(son_parent_IDs, [mother father], dims=1)
                        sons_per_mother[mother] += 1
                        sons_per_father[father] += 1
                        locations_sons = [locations_sons; disperse_individual(locations_F[mother], sigma_disp, geographic_limits)]
                    end
                end
            end
        end # of loop through mothers

        # check if either no daughters or no sons, and end the simulation if so
        if (size(genotypes_daughters, 3) == 0) || (size(genotypes_sons, 3) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end

        # For someday: add in here the option of tracking fitness?

        # determine survival fitnesses of daughters due to epistasis

        survival_fitness_daughters = get_survival_fitnesses(genotypes_daughters[:, hybrid_survival_loci, :], w_hyb)
        daughters_survive = survival_fitness_daughters .> rand(length(survival_fitness_daughters))
        # same for sons:
        survival_fitness_sons = get_survival_fitnesses(genotypes_sons[:, hybrid_survival_loci, :], w_hyb)
        sons_survive = survival_fitness_sons .> rand(length(survival_fitness_sons))

        # check if either no surviving daughters or no surviving sons, and end the simulation if so
        if (sum(daughters_survive) == 0) || (sum(sons_survive) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end

        if do_plot
            push!(position, plot_density(
                [locations_F; locations_M],
                [calc_traits_additive(genotypes_F); calc_traits_additive(genotypes_M)],
                sigma_comp
            ))
        end

        # assign surviving offspring to new adult population
        genotypes_F = genotypes_daughters[:, :, daughters_survive]
        locations_F = locations_daughters[daughters_survive]
        genotypes_M = genotypes_sons[:, :, sons_survive]
        locations_M = locations_sons[sons_survive]

        # update the plot

    end # of loop through generations
    return position

    return genotypes_F, locations_F, genotypes_M, locations_M, extinction
end

function density_at_point(focal_location, locations, genotypes, population, sigma_comp)
    function in_pop(genotype)
        return genotype == population
    end
    indices = findall(in_pop, genotypes)

    return sum(exp.(-((locations[indices] .- focal_location) .^ 2) ./ (2 * (sigma_comp^2))))
end

function plot_density(locations, genotypes, sigma_comp)
    spaced_locations = -2.5:0.01:2.5
    D_0 = density_at_point.(spaced_locations, Ref(locations), Ref(genotypes), Ref(0), Ref(sigma_comp))
    M = maximum(D_0)
    D_0 = D_0 ./ Ref(M)
    plot(spaced_locations, D_0, xlims=(0, 1), ylims=(0, 1))
    D_1 = density_at_point.(spaced_locations, Ref(locations), Ref(genotypes), Ref(1), Ref(sigma_comp))
    D_1 = D_1 ./ Ref(M)
    plot!(spaced_locations, D_1)

    gui()

    return spaced_locations[argmin(abs.(D_0 .- Ref(0.1)))]
end

position = run_one_HZAM_sim(0, 10, 1, 1.5, K_total=10000, sigma_disp=0.05, max_generations=100)

plot(1:50, position)
gui()
readline()