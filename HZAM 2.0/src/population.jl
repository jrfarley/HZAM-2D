module Population

include("plot_data.jl")

using .Plot_Data

using Statistics

export initialize_population, update_population
export calc_female_growth_rate, choose_closest_male, choose_random_male, get_population, calc_match_strength, generate_offspring_genotype, disperse_individual
export plot_population, update_plot
export calc_distance

export get_genotypes, get_growth_rates, get_locations, get_ideal_densities_taylor # for use in tests only

global sympatry # if sympatry is true then the ranges are one-dimensional
global intrinsic_R # population growth rate
global K_total, K_A, K_B # the carrying capacities of each resource
global genotypes_F, genotypes_M # the genotype of each individual
global locations_F, locations_M # the locations of each individual
global growth_rate_resourceA, growth_rate_resourceB # the growth rates of each female with respect to a resource
global competition_trait_loci, female_mating_trait_loci, male_mating_trait_loci # the loci range for competition trait, female mating trait and male mating trait
global ind_useResourceA_F, ind_useResourceB_F # individual contributions to resource use (female)
global ind_useResourceA_M, ind_useResourceB_M # individual contributions to resource use (male)
global competAbility_useResourceA_pop0, competAbility_useResourceA_pop1, competAbility_useResourceB_pop0, competAbility_useResourceB_pop1 # the ability of each population to use each resource
global geographic_limits # the geographic range of the simulation
global functional_loci_range, hybrid_survival_loci

# generates the starting genotypes/locations
# and calculates the growth rates based on individual resource use
function initialize_population(new_K_total::Int,
    starting_pop_ratio,
    ecolDiff,
    new_sympatry::Bool,
    total_loci::Int,
    new_competition_trait_loci,
    new_female_mating_trait_loci,
    new_male_mating_trait_loci,
    new_hybrid_survival_loci,
    new_intrinsic_R,
    new_sigma_comp;
    new_geographic_limits::Vector{}=[0.0, 1.0],
    starting_range_pop0=[0.0, 0.48],
    starting_range_pop1=[0.52, 1.0])

    # initialize the global variables
    global competition_trait_loci = new_competition_trait_loci
    global female_mating_trait_loci = new_female_mating_trait_loci
    global male_mating_trait_loci = new_male_mating_trait_loci
    global hybrid_survival_loci = new_hybrid_survival_loci
    global sympatry = new_sympatry
    global intrinsic_R = new_intrinsic_R
    global sigma_comp = new_sigma_comp
    global geographic_limits = new_geographic_limits
    global spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])

    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    global competAbility_useResourceA_pop0 = (1 + ecolDiff) / 2    # equals 1 when ecolDiff = 1   
    global competAbility_useResourceB_pop0 = 1 - competAbility_useResourceA_pop0
    global competAbility_useResourceA_pop1 = (1 - ecolDiff) / 2   # equals 0 when ecolDiff = 1
    global competAbility_useResourceB_pop1 = 1 - competAbility_useResourceA_pop1

    # set up carying capacities on each resource, and starting pop sizes of each species
    global K_total = new_K_total
    global K_A = convert(Int64, K_total / 2)  # EVEN NUMBER; carrying capacity (on resource alpha) of entire range (for two sexes combined), regardless of species 
    global K_B = K_A   # EVEN NUMBER; carrying capacity (on resource beta) of entire range (for two sexes combined), regardless of species

    total_functional_loci = max(maximum(female_mating_trait_loci), maximum(male_mating_trait_loci), maximum(competition_trait_loci), maximum(hybrid_survival_loci))
    global functional_loci_range = 1:total_functional_loci

    if sympatry  # if sympatry is true, then set both pop sizes according to full range 
        pop0_starting_N = K_A   # starting N of species 0
        pop0_starting_N_half = Int(pop0_starting_N / 2)  # The "Int" is to ensure no decimal
        pop1_starting_N = Int(round(starting_pop_ratio * K_B))   # starting N of species 1, which can be lower if starting_pop_ratio is below 1)
        pop1_starting_N_half = Int(pop1_starting_N / 2)
    else  # sympatry is false, then set pop sizes assuming no range overlap--this is why the "2" is in the formulae below  
        pop0_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop0) + (K_B * competAbility_useResourceB_pop0)) * (starting_range_pop0[2] - starting_range_pop0[1]) / (geographic_limits[2] - geographic_limits[1]))
        pop0_starting_N_half = Int(pop0_starting_N / 2)  # starting number for each sex
        pop1_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop1) + (K_B * competAbility_useResourceB_pop1)) * (starting_range_pop1[2] - starting_range_pop1[1]) / (geographic_limits[2] - geographic_limits[1]))
        println(pop0_starting_N_half)
        println(pop1_starting_N)
        pop1_starting_N_half = Int(pop1_starting_N / 2)
    end

    # Generate genotype array for population of females:
    # this is a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
    # columns (D2) are loci, and pages (D3) are individuals
    global genotypes_F = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    global genotypes_M = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    # functional loci are first, followed by neutral loci

    # initialize growth rates
    if sympatry
        global growth_rate_resourceA, growth_rate_resourceB = calculate_initial_growth_rates_sympatry(K_A, K_B)
    else
        global locations_F, locations_M = generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)
        global ind_useResourceA_F = [fill(competAbility_useResourceA_pop0, pop0_starting_N_half); fill(competAbility_useResourceA_pop1, pop1_starting_N_half)]
        global ind_useResourceA_M = ind_useResourceA_F
        global ind_useResourceB_F = [fill(competAbility_useResourceB_pop0, pop0_starting_N_half); fill(competAbility_useResourceB_pop1, pop1_starting_N_half)]
        global ind_useResourceB_M = ind_useResourceB_F

        global growth_rate_resourceA, growth_rate_resourceB = calculate_growth_rates_spatial()
    end
end

function generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)
    # Generate breeding locations of individuals (vector in same order of individuals as D3 of the genotype array)
    locations_F_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_F_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_M_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_M_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    return [locations_F_pop0; locations_F_pop1], [locations_M_pop0; locations_M_pop1]
end

# updates the genotypes, locations, and growth rates with the values for the next generation
function update_population(genotypes_daughters::Array{Int8,3},
    genotypes_sons::Array{Int8,3},
    locations_daughters::Array{Float32},
    locations_sons::Array{Float32})

    global genotypes_F = genotypes_daughters
    global genotypes_M = genotypes_sons

    global locations_F = locations_daughters
    global locations_M = locations_sons


    # calculate ecological competition trait values from genotypes
    competition_traits_F = calc_traits_additive(genotypes_F[:, competition_trait_loci, :])
    competition_traits_M = calc_traits_additive(genotypes_M[:, competition_trait_loci, :])

    global growth_rate_resourceA, growth_rate_resourceB = calculate_growth_rates(competition_traits_F, competition_traits_M)
end

# calculates the growth rates at the beginning of the simulation when the populations are in sympatry
function calculate_initial_growth_rates_sympatry(K_A, K_B)
    total_useResourceA = competAbility_useResourceA_pop0 * K_A + competAbility_useResourceA_pop1 * K_B
    total_useResourceB = competAbility_useResourceB_pop0 * K_A + competAbility_useResourceB_pop1 * K_B
    # calculate global growth rates due to each resource (according to discrete time logistic growth equation)
    growth_rate_resourceA = (intrinsic_R * K_A) / (K_A + ((total_useResourceA) * (intrinsic_R - 1)))
    growth_rate_resourceB = (intrinsic_R * K_B) / (K_B + ((total_useResourceB) * (intrinsic_R - 1)))

    return growth_rate_resourceA, growth_rate_resourceB
end

# calculates growth rates
function calculate_growth_rates(competition_traits_F, competition_traits_M)

    # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1

    global ind_useResourceA_F = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_F) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
    global ind_useResourceB_F = competAbility_useResourceB_pop0 .+ (competition_traits_F .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))
    global ind_useResourceA_M = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_M) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
    global ind_useResourceB_M = competAbility_useResourceB_pop0 .+ (competition_traits_M .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))

    if sympatry
        return calculate_growth_rates_sympatry()
    else
        return calculate_growth_rates_spatial()
    end

end

# calculates growth rates when the populations are in sympatry
function calculate_growth_rates_sympatry()
    # calculate growth rate in sympatric model
    # sum up the global resource use over all individuals
    total_useResourceA = sum(ind_useResourceA_F) + sum(ind_useResourceA_M)
    total_useResourceB = sum(ind_useResourceB_F) + sum(ind_useResourceB_M)
    # calculate global growth rates due to each resource (according to discrete time logistic growth equation)
    growth_rate_resourceA = (intrinsic_R * K_A) / (K_A + ((total_useResourceA) * (intrinsic_R - 1)))
    growth_rate_resourceB = (intrinsic_R * K_B) / (K_B + ((total_useResourceB) * (intrinsic_R - 1)))

    return growth_rate_resourceA, growth_rate_resourceB
end

# calculates the ideal density (assuming normal distribution) at the location of each female
function get_ideal_densities_taylor(K_total, sigma_comp, locations_F)
    #approximation of integral from 0 to 1 of K_total*exp(-(x-focal_location)^2/(2*sigma_comp^2)) with respect to x

    maximum_density = (K_total / 24) * (2 * exp(-(1 / (8 * sigma_comp^2))) +
                                        8 * exp(-((0.375^2) / (2 * sigma_comp^2))) +
                                        4 * exp(-(1 / (32 * sigma_comp^2))) +
                                        8 * exp(-(1 / (128 * sigma_comp^2))) + 2)

    first_coefficient_taylor_series = ((K_total * exp(-1 / (8 * sigma_comp^2))) / (2 * sigma_comp^2))

    second_coefficient_taylor_series = K_total * (((3 * exp(-1 / (8 * sigma_comp^2))) / (24 * sigma_comp^4)) -
                                                  (exp(-1 / (8 * sigma_comp^2)) / (96 * sigma_comp^6)))

    function get_ideal_density(focal_location)
        return 0.5 * (maximum_density - first_coefficient_taylor_series * (focal_location - 0.5)^2 + second_coefficient_taylor_series * (focal_location - 0.5)^4)
    end

    return map(get_ideal_density, locations_F)
end

# calculates growth rates in the spatial model
function calculate_growth_rates_spatial()
    # in spatial model, calculate growth rates based on local resource use

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity

    ideal_densities_at_locations_F = get_ideal_densities_taylor(K_total, sigma_comp, locations_F) # this applies the above function to each geographic location

    # assume both resources have same constant density across range

    # determine local resource use for each location across range
    ind_useResourceA_all = [ind_useResourceA_F; ind_useResourceA_M]
    ind_useResourceB_all = [ind_useResourceB_F; ind_useResourceB_M]
    ind_locations_real = [locations_F; locations_M]
    function get_useResourceA_density_real(focal_location) # this function calculates local density according to a normal curve, weighted by individual resource use
        sum(ind_useResourceA_all .* exp.(-((ind_locations_real .- focal_location) .^ 2) ./ (2 * (sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end

    real_densities_at_locations_F_resourceA = map(get_useResourceA_density_real, locations_F) # this applies the above function to each geographic location


    function get_useResourceB_density_real(focal_location) # do the same for resource B
        sum(ind_useResourceB_all .* exp.(-((ind_locations_real .- focal_location) .^ 2) ./ (2 * (sigma_comp^2))))
    end
    real_densities_at_locations_F_resourceB = map(get_useResourceB_density_real, locations_F) # this applies the above function to each geographic location 

    # calculate local growth rates due to each resource (according to discrete time logistic growth equation)

    local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_locations_F ./ (ideal_densities_at_locations_F .+ ((real_densities_at_locations_F_resourceA) .* (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_locations_F ./ (ideal_densities_at_locations_F .+ ((real_densities_at_locations_F_resourceB) .* (intrinsic_R - 1)))
    #println(calc_local_growth_rate(ind_locations_real[29]))

    #println(local_growth_rates_resourceA[29])

    return local_growth_rates_resourceA, local_growth_rates_resourceB
end

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
    N = size(genotypes, 3)
    #traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64

    function mean(array)
        return sum(array) / length(array)
    end

    traits = map(x -> mean(genotypes[:, :, x]), 1:N)
    return traits
end

function calc_mating_trait_female(n)
    return mean(genotypes_F[:, female_mating_trait_loci, n])
end

function calc_mating_trait_male(n)
    return mean(genotypes_M[:, female_mating_trait_loci, n])
end

# These next two functions define the way potential male mates are chosen.
# elig_M is a vector of indices of possible male mates.
# When in sympatric model, choose random male (note the second and third arguments not used but allows function to be called in same way as in spatial model):
function choose_random_male(elig_M::Vector{UInt32}, female_index::Real)
    focal_male = splice!(elig_M, rand(eachindex(elig_M))) # this gets the index of a random male, and removes that male from the list in elig_M
    return focal_male, elig_M
end

# When in spatial model, choose closest male:
function choose_closest_male(elig_M::Vector{UInt32}, female_index::Real)
    focal_male = splice!(elig_M, argmin(abs.(locations_M[elig_M] .- locations_F[female_index]))) # this gets the index of a closest male, and removes that male from the list in elig_M
    return focal_male, elig_M
end

# calculates the distance between the mother and the father
function calc_distance(mother, father)
    return abs(locations_F[mother] - locations_M[father])
end

# compare male trait with female's trait (preference), and determine
# whether she accepts; note that match_strength is determined by a
# Gaussian, with a maximum of 1 and minimum of zero
function calc_match_strength(focal_female_index, focal_male_index, pref_SD)
    mating_trait_dif = calc_mating_trait_male(focal_male_index) - calc_mating_trait_female(focal_female_index)
    return exp((-(mating_trait_dif^2)) / (2 * (pref_SD^2)))
end

# determine fitness due to female use of available resources
function calc_female_growth_rate(focal_female_index)
    if sympatry
        return (ind_useResourceA_F[focal_female_index] * growth_rate_resourceA) + (ind_useResourceB_F[focal_female_index] * growth_rate_resourceB)
    else
        return (ind_useResourceA_F[focal_female_index] * growth_rate_resourceA[focal_female_index]) + (ind_useResourceB_F[focal_female_index] * growth_rate_resourceB[focal_female_index])
    end
end

# returns the total number of females and the total number of males 
function get_population()
    return size(genotypes_F, 3), size(genotypes_M, 3)
end

# generates the genotype for an offspring based on its parents' genotypes
function generate_offspring_genotype(female_index, male_index, total_loci)
    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
    kid_info = Array{Int8,2}(undef, 2, total_loci)
    for locus in 1:total_loci
        kid_info[1, locus] = genotypes_F[rand([1 2]), locus, female_index]  # for this locus, pick a random allele from the mother
        kid_info[2, locus] = genotypes_M[rand([1 2]), locus, male_index]  # and from the father
    end
    return kid_info
end

# This function determines breeding location of one individual,
# based on a normal distribution with width sigma_disp,
# centred on birth location. Constrained to be within range. 
# geographic_limits should be a vector with two numbers.
function disperse_individual(female_index::Int, sigma_disp::Real, geographic_limits::Vector{})::Float32
    while true
        new_location = locations_F[female_index] + (sigma_disp * randn())
        if (new_location >= geographic_limits[1]) & (new_location <= geographic_limits[2]) # checks if location is in range
            return new_location
        end
    end
end

# returns the hybrid index of every individual in the population
function calc_hybrid_indices()
    return [calc_traits_additive(genotypes_F[:, functional_loci_range, :]); calc_traits_additive(genotypes_M[:, functional_loci_range, :])]
end


function plot_population()
    create_new_plot(calc_hybrid_indices(), [locations_F; locations_M])
end

function update_plot(generation)
    update_population_plot(calc_hybrid_indices(), [locations_F; locations_M], generation)
end

function get_genotypes()
    return genotypes_F, genotypes_M
end

function get_growth_rates()
    return growth_rate_resourceA, growth_rate_resourceB
end

function get_locations()
    return locations_F, locations_M
end

end