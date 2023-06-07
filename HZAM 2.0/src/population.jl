module Population

include("plot_data.jl")

using .Plot_Data

using Statistics

export initialize_population, update_population
export calc_female_growth_rate, choose_closest_male, choose_random_male, get_population, calc_match_strength, generate_offspring_genotype, disperse_individual
export plot_population, update_plot
export active_F, active_M, left_boundary, right_boundary
export calc_distance, get_num_hybridization_events
export get_results
using Test

global sympatry # if sympatry is true then the ranges are one-dimensional
global intrinsic_R # population growth rate
global K_total, K_A, K_B # the carrying capacities of each resource
global genotypes_F, genotypes_M # the genotype of each individual
global locations_F, locations_M # the locations of each individual
global growth_rate_resourceA, growth_rate_resourceB # the growth rates of each female with respect to a resource
global competition_trait_loci, female_mating_trait_loci, male_mating_trait_loci, neutral_loci # the loci ranges for competition trait, female mating trait, male mating trait, and the neutral loci
global ind_useResourceA_F, ind_useResourceB_F # individual contributions to resource use (female)
global ind_useResourceA_M, ind_useResourceB_M # individual contributions to resource use (male)
global competAbility_useResourceA_pop0, competAbility_useResourceA_pop1, competAbility_useResourceB_pop0, competAbility_useResourceB_pop1 # the ability of each population to use each resource
global geographic_limits # the geographic range of the simulation
global functional_loci_range, hybrid_survival_loci
global num_hybridization_events = 0
global active_region_limits = [0.4, 0.6]
global active_F, active_M = [], []
global inactive_F, inactive_M = [], []
global buffer_M = [[], []]
global individuals_per_zone_F = [[] for i in 1:10]
global individuals_per_zone_M = [[] for i in 1:10]
global left_boundary = 5
global right_boundary = 6
global fixed_locations_F, fixed_locations_M, fixed_genotypes_F, fixed_genotypes_M
global mitochondria_F, mitochondria_M, fixed_mitochondria_F, fixed_mitochondria_M


function get_num_hybridization_events()
    return num_hybridization_events
end

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
    global neutral_loci = setdiff(collect(1:total_loci), [competition_trait_loci; female_mating_trait_loci; male_mating_trait_loci; hybrid_survival_loci])
    global sympatry = new_sympatry
    global intrinsic_R = new_intrinsic_R
    global sigma_comp = new_sigma_comp
    global geographic_limits = new_geographic_limits
    #global spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])

    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    global competAbility_useResourceA_pop0 = (1 + ecolDiff) / 2    # equals 1 when ecolDiff = 1   
    global competAbility_useResourceB_pop0 = 1 - competAbility_useResourceA_pop0
    global competAbility_useResourceA_pop1 = competAbility_useResourceB_pop0   # equals 0 when ecolDiff = 1
    global competAbility_useResourceB_pop1 = competAbility_useResourceA_pop0

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
        pop1_starting_N_half = Int(pop1_starting_N / 2)
    end

    # Generate genotype array for population of females:
    # this is a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
    # columns (D2) are loci, and pages (D3) are individuals
    global genotypes_F = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    global genotypes_M = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)

    global mitochondria_F = [fill(0.25, pop0_starting_N_half); fill(0.75, pop1_starting_N_half)]
    global mitochondria_M = [fill(0.25, pop0_starting_N_half); fill(0.75, pop1_starting_N_half)]
    # functional loci are first, followed by neutral loci

    # initialize growth rates
    if sympatry
        global growth_rate_resourceA, growth_rate_resourceB = calculate_initial_growth_rates_sympatry(K_A, K_B)
        global active_F = 1:length(genotypes_F)
        global active_M = 1:length(genotypes_M)
    else
        global locations_F, locations_M = generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)
        global individuals_per_zone_F, individuals_per_zone_M = assign_zones(locations_F, locations_M)
        global active_F = vcat(individuals_per_zone_F[5], individuals_per_zone_F[6])
        global inactive_F = vcat(individuals_per_zone_F[1],
            individuals_per_zone_F[2],
            individuals_per_zone_F[3],
            individuals_per_zone_F[4],
            individuals_per_zone_F[7],
            individuals_per_zone_F[8],
            individuals_per_zone_F[9],
            individuals_per_zone_F[10])
        global buffer_M = [individuals_per_zone_M[4], individuals_per_zone_M[7]]
        global active_M = vcat(individuals_per_zone_M[5], individuals_per_zone_M[6], buffer_M[1], buffer_M[2])
        global inactive_M = vcat(individuals_per_zone_M[1],
            individuals_per_zone_M[2],
            individuals_per_zone_M[3],
            individuals_per_zone_M[8],
            individuals_per_zone_M[9],
            individuals_per_zone_M[10])

        global ind_useResourceA_F = [fill(competAbility_useResourceA_pop0, pop0_starting_N_half); fill(competAbility_useResourceA_pop1, pop1_starting_N_half)]
        global ind_useResourceA_M = ind_useResourceA_F
        global ind_useResourceB_F = [fill(competAbility_useResourceB_pop0, pop0_starting_N_half); fill(competAbility_useResourceB_pop1, pop1_starting_N_half)]
        global ind_useResourceB_M = ind_useResourceB_F

        global growth_rate_resourceA, growth_rate_resourceB = calculate_growth_rates_spatial()
        global fixed_locations_F, fixed_locations_M = (copy(locations_F), copy(locations_M))
    end

    global fixed_genotypes_F, fixed_genotypes_M = (copy(genotypes_F), copy(genotypes_M))
    global fixed_mitochondria_F, fixed_mitochondria_M = (copy(mitochondria_F), copy(mitochondria_M))
end


# Generate breeding locations of individuals (vector in same order of individuals as D3 of the genotype array)
function generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)
    locations_F_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_F_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_M_pop0 = Array{Float32,1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1])
    locations_M_pop1 = Array{Float32,1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    return [locations_F_pop0; locations_F_pop1], [locations_M_pop0; locations_M_pop1]
end

# arranges the indices of each individual by location
# zone 1 is [0,1), zone 2 is [1, 2), etc.
function assign_zones(locations_F, locations_M)
    new_individuals_per_zone_F = [[] for i in 1:10]
    new_individuals_per_zone_M = [[] for i in 1:10]

    # assigns females
    for indv in eachindex(locations_F)
        zone = trunc(Int, 10 * locations_F[indv]) + 1
        push!(new_individuals_per_zone_F[zone], indv)
    end

    #assigns males
    for indv in eachindex(locations_M)
        zone = trunc(Int, 10 * locations_M[indv]) + 1
        push!(new_individuals_per_zone_M[zone], indv)
    end

    return new_individuals_per_zone_F, new_individuals_per_zone_M
end

# updates the genotypes, locations, and growth rates with the values for the next generation
function update_population(genotypes_daughters,
    genotypes_sons,
    mitochondria_daughters,
    mitochondria_sons,
    locations_daughters::Array{Float32},
    locations_sons::Array{Float32},
    expand_left,
    expand_right,
    generation)

    new_active_F = [] # keeps track of inactive females that will be active in the next generation
    new_active_M = [] # keeps track of inactive males that will be active in the next generation

    if generation % 10 == 0
        # every 10 generations the active and inactive zones are reset to reflect current distribution of genotypes

        # compile list of genotypes/locations of offspring and fixed individuals
        global genotypes_F = vcat(genotypes_daughters, fixed_genotypes_F[inactive_F])
        global genotypes_M = vcat(genotypes_sons, fixed_genotypes_M[buffer_M[1]], fixed_genotypes_M[buffer_M[2]], fixed_genotypes_M[inactive_M])

        global mitochondria_F = vcat(mitochondria_daughters, fixed_mitochondria_F[inactive_F])
        global mitochondria_M = vcat(mitochondria_sons, fixed_mitochondria_M[buffer_M[1]], fixed_mitochondria_M[buffer_M[2]], fixed_mitochondria_M[inactive_M])

        global locations_F = vcat(locations_daughters, fixed_locations_F[inactive_F])
        global locations_M = vcat(locations_sons, fixed_locations_M[buffer_M[1]], fixed_locations_M[buffer_M[2]], fixed_locations_M[inactive_M])

        # update fixed individuals to current population
        global fixed_locations_F, fixed_locations_M = copy(locations_F), copy(locations_M)
        global fixed_genotypes_F, fixed_genotypes_M = copy(genotypes_F), copy(genotypes_M)
        global fixed_mitochondria_F, fixed_mitochondria_M = copy(mitochondria_F), copy(mitochondria_M)

        # group population into zones
        global individuals_per_zone_F, individuals_per_zone_M = assign_zones(locations_F, locations_M)

        dead_zones = find_inactive_zones(individuals_per_zone_F, individuals_per_zone_M)
        active_zones = setdiff!(collect(1:10), dead_zones)

        # compile list of new inactive individuals
        if dead_zones != []
            global inactive_F = reduce(vcat, individuals_per_zone_F[dead_zones])
            global inactive_M = reduce(vcat, individuals_per_zone_M[dead_zones])
            global buffer_M = [[], []]
        else
            global inactive_F = []
            global inactive_M = []
        end

        # compile list of new active individuals
        if active_zones != []
            global active_F = reduce(vcat, individuals_per_zone_F[active_zones])
            global active_M = reduce(vcat, individuals_per_zone_M[active_zones])
            global left_boundary = minimum(vcat(10, active_zones))
            global right_boundary = maximum(vcat(1, active_zones))
        else
            global active_F = []
            global active_M = []
        end
    else
        # expands active zone and updates buffer if necessary
        if expand_left && left_boundary > 1 # expands the left boundary if not already at edge of range
            global left_boundary -= 1
            new_active_F = individuals_per_zone_F[left_boundary]
            new_active_M = individuals_per_zone_M[left_boundary]
            if left_boundary > 1
                global buffer_M[1] = individuals_per_zone_M[left_boundary-1]
            else
                global buffer_M[1] = []
            end
        end

        if expand_right && right_boundary < 10 # expands the right boundary if not already at edge of range
            global right_boundary += 1
            new_active_F = vcat(new_active_F, individuals_per_zone_F[right_boundary])
            new_active_M = vcat(new_active_M, individuals_per_zone_M[right_boundary])
            if right_boundary < 10
                global buffer_M[2] = individuals_per_zone_M[right_boundary+1]
            else
                global buffer_M[2] = []
            end
        end

        new_active_M = vcat(new_active_M, buffer_M[1], buffer_M[2]) # add males in the buffer zone to the list of active males

        num_fixed_F = reduce(+, map(length, individuals_per_zone_F))
        num_fixed_M = reduce(+, map(length, individuals_per_zone_M))

        @test num_fixed_F == length(fixed_locations_F)
        @test num_fixed_M == length(fixed_locations_M)
        @test maximum(vcat(0, new_active_F)) <= length(fixed_locations_F)
        @test maximum(vcat(0, new_active_M)) <= length(fixed_locations_M)
        @test maximum(vcat(0, inactive_F)) <= length(fixed_locations_F)
        @test maximum(vcat(0, inactive_M)) <= length(fixed_locations_M)

        # updates list of inactive individuals
        global inactive_F = setdiff(inactive_F, new_active_F)
        global inactive_M = setdiff(inactive_M, new_active_M)

        # creates new lists of genotypes and locations so that the active individuals are first and the inactive individuals last 
        global genotypes_F = vcat(genotypes_daughters, fixed_genotypes_F[new_active_F], fixed_genotypes_F[inactive_F])
        global genotypes_M = vcat(genotypes_sons, fixed_genotypes_M[new_active_M], fixed_genotypes_M[inactive_M])

        global mitochondria_F = vcat(mitochondria_daughters, fixed_mitochondria_F[new_active_F], fixed_mitochondria_F[inactive_F])
        global mitochondria_M = vcat(mitochondria_sons, fixed_mitochondria_M[new_active_M], fixed_mitochondria_M[inactive_M])

        global locations_F = vcat(locations_daughters, fixed_locations_F[new_active_F], fixed_locations_F[inactive_F])
        global locations_M = vcat(locations_sons, fixed_locations_M[new_active_M], fixed_locations_M[inactive_M])


        # updates list of active individuals
        num_active_F = length(locations_daughters) + length(new_active_F)
        num_active_M = length(locations_sons) + length(new_active_M)

        global active_F = collect(1:num_active_F)
        global active_M = collect(1:num_active_M)

    end

    # calculate new competition traits
    competition_traits_F = calc_traits_additive(genotypes_F, competition_trait_loci)
    competition_traits_M = calc_traits_additive(genotypes_M, competition_trait_loci)

    # calculate new growth rates
    global growth_rate_resourceA, growth_rate_resourceB = calculate_growth_rates(competition_traits_F, competition_traits_M)
end

# finds which parts of the range contain only individuals of one genotype
function find_inactive_zones(indv_per_zone_F, indv_per_zone_M)
    dead_zones = []
    output = []
    zone_classification = []

    # check in each zone if all individuals have the same genotypes
    for zone in 1:10
        if length(indv_per_zone_M[zone]) > 0 && length(indv_per_zone_F[zone]) > 0
            all_equal = true
            genotype1 = genotypes_M[indv_per_zone_M[zone][1]] # genotype of the first male in the zone
            mitochondria1 = mitochondria_M[indv_per_zone_M[zone][1]]

            push!(zone_classification, genotype1[1, 1])

            # check each female
            for indv in indv_per_zone_F[zone]
                if genotypes_F[indv] != genotype1 || mitochondria_F[indv] != mitochondria1
                    all_equal = false
                    break
                end
            end
            if all_equal
                # check each male
                for indv in indv_per_zone_M[zone]
                    if genotypes_M[indv] != genotype1 || mitochondria_M[indv] != mitochondria1
                        all_equal = false
                        break
                    end
                end
                if all_equal
                    push!(dead_zones, zone)
                end
            end
        else
            push!(dead_zones, zone)
            push!(zone_classification, 0.5)
        end
    end

    # goes through each zone where all individuals have the same genotype and keeps the zones where 
    # all the individuals in the neighbouring zones also have the same genotype
    for zone in dead_zones
        left_zone = max(1, zone - 1)
        right_zone = min(10, zone + 1)
        if left_zone in dead_zones && right_zone in dead_zones
            if abs(zone_classification[left_zone] - zone_classification[zone]) < 1 &&
               abs(zone_classification[right_zone] - zone_classification[zone]) < 1
                push!(output, zone)
            end
        end
    end

    return output
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

# calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
function calculate_ind_useResource(competition_traits_F, competition_traits_M)
    ecolDiff = competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1

    ind_useResourceA_F = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_F) .* ecolDiff)
    ind_useResourceB_F = competAbility_useResourceB_pop0 .+ (competition_traits_F .* ecolDiff)
    ind_useResourceA_M = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_M) .* ecolDiff)
    ind_useResourceB_M = competAbility_useResourceB_pop0 .+ (competition_traits_M .* ecolDiff)

    return ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M
end

# calculates growth rates
function calculate_growth_rates(competition_traits_F, competition_traits_M)

    # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1

    global ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = calculate_ind_useResource(competition_traits_F, competition_traits_M)

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
function get_ideal_densities(K_total, sigma_comp, locations_F)
    function erf(x) # approximation of the error function
        return tanh((sqrt(pi) * log(2)) * x)
    end

    function calc_ideal_density(location) # integral from 0 to 1 of K_total*exp(-(x-focal_location)^2/(2*sigma_comp^2)) with respect to x
        return K_total * sqrt(pi / 2) * sigma_comp * (erf((1 - location) / (sqrt(2) * sigma_comp)) + erf(location / (sqrt(2) * sigma_comp)))
    end

    return map(calc_ideal_density, locations_F)
end

# calculates growth rates in the spatial model
function calculate_growth_rates_spatial()
    # in spatial model, calculate growth rates based on local resource use

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity

    ideal_densities_at_locations_F = get_ideal_densities(K_total, sigma_comp, locations_F) # this applies the above function to each geographic location

    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

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

    local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_locations_F_resourceA ./ (ideal_densities_at_locations_F_resourceA .+ ((real_densities_at_locations_F_resourceA) .* (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_locations_F_resourceB ./ (ideal_densities_at_locations_F_resourceB .+ ((real_densities_at_locations_F_resourceB) .* (intrinsic_R - 1)))

    return local_growth_rates_resourceA, local_growth_rates_resourceB
end

# This function sets up the genotypes of the starting population
# in a 3D array, where:
# rows (D1) are alleles (row 1 from mother, row 2 from father),
# columns (D2) are loci, 
# pages (D3) are individuals.  
function generate_genotype_array(N_pop0::Integer, N_pop1::Integer, loci::Integer)
    total_N = N_pop0 + N_pop1
    genotypes = vcat([fill(0, 2, loci) for i in 1:N_pop0], [fill(1, 2, loci) for j in (N_pop0+1):total_N])

    return genotypes
end

# This function calculates each mean values of the genotypes passed to it (for each individual).
# Used to determine trait values in an additive way.
# Only those loci that are additive trait loci should be passed to this function.
function calc_traits_additive(genotypes, loci)::Vector{Float32} #=::Array{Int8,3}=#
    N = length(genotypes)
    #traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64

    function mean(genotype, loci)
        return sum(genotype[:, loci]) / (2 * length(loci))
    end

    traits = map(x -> mean(genotypes[x], loci), 1:N)
    return traits
end

function calc_mating_trait_female(n)
    return mean(genotypes_F[n][:, female_mating_trait_loci])
end

function calc_mating_trait_male(n)
    return mean(genotypes_M[n][:, female_mating_trait_loci])
end

# These next two functions define the way potential male mates are chosen.
# elig_M is a vector of indices of possible male mates.
# When in sympatric model, choose random male (note the second and third arguments not used but allows function to be called in same way as in spatial model):
function choose_random_male(elig_M, female_index::Real)
    focal_male = splice!(elig_M, rand(eachindex(elig_M))) # this gets the index of a random male, and removes that male from the list in elig_M
    return focal_male, elig_M
end

# When in spatial model, choose closest male:
function choose_closest_male(elig_M, female_index::Real)
    focal_male = splice!(elig_M, argmin(abs.(locations_M[elig_M] .- locations_F[female_index]))) # this gets the index of a closest male, and removes that male from the list in elig_M
    if (calc_distance(female_index, focal_male) >= 0.1)
        return focal_male, []
    end
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
        growth_rateA = ind_useResourceA_F[focal_female_index] * growth_rate_resourceA
        growth_rateB = ind_useResourceB_F[focal_female_index] * growth_rate_resourceB
        return growth_rateA + growth_rateB
    else
        growth_rateA = ind_useResourceA_F[focal_female_index] * growth_rate_resourceA[focal_female_index]
        growth_rateB = ind_useResourceB_F[focal_female_index] * growth_rate_resourceB[focal_female_index]
        return growth_rateA + growth_rateB
    end
end

# returns the total number of females and the total number of males 
function get_population()
    return length(genotypes_F), length(genotypes_M)
end

# generates the genotype for an offspring based on its parents' genotypes
function generate_offspring_genotype(female_index, male_index, total_loci)
    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
    kid_genotype = Array{Int8,2}(undef, 2, total_loci)
    kid_mitochondria = mitochondria_F[female_index]
    for locus in 1:total_loci
        kid_genotype[1, locus] = genotypes_F[female_index][rand([1 2]), locus]  # for this locus, pick a random allele from the mother
        kid_genotype[2, locus] = genotypes_M[male_index][rand([1 2]), locus]  # and from the father
    end

    if (kid_genotype != genotypes_F[female_index][:, :])
        global num_hybridization_events += 1
    end

    return kid_genotype, kid_mitochondria
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
function calc_hybrid_indices(genotypes)
    return calc_traits_additive(genotypes, functional_loci_range)
end

# creates a plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function plot_population()
    create_new_plot(calc_hybrid_indices([genotypes_F[active_F]; genotypes_M[active_M]]),
        [mitochondria_F[active_F]; mitochondria_M[active_M]],
        [locations_F[active_F]; locations_M[active_M]],
        calc_hybrid_indices([fixed_genotypes_F[inactive_F]; fixed_genotypes_M[inactive_M]]),
        [fixed_mitochondria_F[inactive_F]; fixed_mitochondria_M[inactive_M]],
        [fixed_locations_F[inactive_F]; fixed_locations_M[inactive_M]])
end

# updates the plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function update_plot(generation)
    genotypes_active = [genotypes_F[active_F]; genotypes_M[active_M]]
    locations_active = [locations_F[active_F]; locations_M[active_M]]
    genotypes_inactive = [fixed_genotypes_F[inactive_F]; fixed_genotypes_M[inactive_M]]
    locations_inactive = [fixed_locations_F[inactive_F]; fixed_locations_M[inactive_M]]
    mitochondria_active = [mitochondria_F[active_F]; mitochondria_M[active_M]]
    mitochondria_inactive = [fixed_mitochondria_F[inactive_F]; fixed_mitochondria_M[inactive_M]]
    update_population_plot(calc_hybrid_indices(genotypes_active),
        mitochondria_active,
        locations_active,
        calc_hybrid_indices(genotypes_inactive),
        mitochondria_inactive,
        locations_inactive,
        generation)
end

# updates the male and female locations to the given values (FOR TESTING ONLY)
function set_locations(new_locations_F, new_locations_M)
    global locations_F = new_locations_F
    global locations_M = new_locations_M
end

# updates the male and female genotypes to the given values (FOR TESTING ONLY)
function set_genotypes(new_genotypes_F, new_genotypes_M)
    global genotypes_F = new_genotypes_F
    global genotypes_M = new_genotypes_M
end

function get_results(extinction)
    functional_HI_all_inds = [calc_traits_additive(genotypes_F, functional_loci_range); calc_traits_additive(genotypes_M, functional_loci_range)]
    HI_NL_all_inds = [calc_traits_additive(genotypes_F, neutral_loci); calc_traits_additive(genotypes_M, neutral_loci)]
    # calculate proportion of all individuals who are species0 or species1 (defined as low and high 10% of HI distribution, respectively)

    return extinction, functional_HI_all_inds, HI_NL_all_inds
end

end