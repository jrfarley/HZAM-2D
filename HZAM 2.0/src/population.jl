module Population

include("plot_data.jl")

using .Plot_Data

using Test
using SpecialFunctions
using QuadGK

export PopulationData, Location
export initialize_population, update_population
export choose_closest_male, calc_match_strength, generate_offspring_genotype
export plot_population, update_plot
export calc_traits_additive
export assign_zone
export NUM_DEMES

NUM_DEMES = 10

# stores a location with an x coordinate and a y coordinate
struct Location
    x::Float32
    y::Float32

    # creates a location the with given x and y coordinates
    function Location(x::Float32, y::Float32)
        new(x, y)
    end

    # generates a random location within the range given
    function Location(starting_range::Vector{Location})
        x = rand() * (starting_range[2].x - starting_range[1].x) + starting_range[1].x
        y = rand() * (starting_range[2].y - starting_range[1].y) + starting_range[1].y
        new(x, y)
    end

    # generates a new location
    # based on a normal distribution with width sigma_disp,
    # centred on the given location location. Constrained to be within range. 
    # geographic_limits should be a vector with two locations.
    function Location(starting_location::Location, sigma_disp::Real, geographic_limits::Vector{Location})
        x = -1
        y = -1
        while (x < geographic_limits[1].x) || (x > geographic_limits[2].x) || (y < geographic_limits[1].y) || (y > geographic_limits[2].y)
            dist = sigma_disp * randn()
            dir = rand() * 2 * pi
            x = starting_location.x + dist * cos(dir)
            y = starting_location.y + dist * sin(dir)
        end
        new(x, y)
    end
end

# stores genotype, location, mitochondria, and resource use data for all individuals in part of the range
struct Deme
    genotypes_F::Vector{Matrix{Int8}} # the female genotypes. rows are alleles (row 1 from mother, row 2 from father) and columns are loci 
    genotypes_M::Vector{Matrix{Int8}} # the male genotypes
    locations_F::Vector{Location} # female locations
    locations_M::Vector{Location} # male locations
    mitochondria_F::Vector{Int8} # female mitochondria types (0 for species 0 and 1 for species 1)
    mitochondria_M::Vector{Int8} # male mitochondria types

    ind_useResourceA_F::Vector{Float64}
    ind_useResourceA_M::Vector{Float64}
    ind_useResourceB_F::Vector{Float64}
    ind_useResourceB_M::Vector{Float64}

    # initialize population in deme
    function Deme(starting_N,
        total_loci,
        location,
        size,
        population, ecolDiff)

        N_half = trunc(Int, starting_N / 2)

        deme_range = [location, Location(min(location.x + size, 0.999f0), min(location.y + size, 0.999f0))]

        genotypes = fill(fill(population, 2, total_loci), N_half)

        locations = [[Location(deme_range) for i in 1:N_half] for i in 1:2]

        mitochondria = fill(population, N_half)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA = calculate_ind_useResourceA.(mitochondria, ecolDiff) # at the beginning the mitochondria value is equal to the competition trait
        ind_useResourceB = calculate_ind_useResourceB.(mitochondria, ecolDiff)

        new(genotypes,
            genotypes,
            locations[1],
            locations[1],
            mitochondria,
            mitochondria,
            ind_useResourceA,
            ind_useResourceA,
            ind_useResourceB,
            ind_useResourceB)
    end

    # updates population in deme with the offspring genotypes, locations, and mitochondria
    function Deme(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M, competition_trait_loci, ecolDiff)
        # calculate new competition traits
        competition_traits_F = calc_traits_additive(genotypes_F, competition_trait_loci)
        competition_traits_M = calc_traits_additive(genotypes_M, competition_trait_loci)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA_F = calculate_ind_useResourceA.(competition_traits_F, ecolDiff)
        ind_useResourceA_M = calculate_ind_useResourceA.(competition_traits_M, ecolDiff)
        ind_useResourceB_F = calculate_ind_useResourceB.(competition_traits_F, ecolDiff)
        ind_useResourceB_M = calculate_ind_useResourceB.(competition_traits_M, ecolDiff)

        new(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M,
            ind_useResourceA_F,
            ind_useResourceA_M,
            ind_useResourceB_F,
            ind_useResourceB_M)
    end
end


# stores all the population data (genotypes, locations, etc.) that get updated each generation
struct PopulationData
    population::Matrix{Deme}

    growth_rates_F::Matrix{Vector{Float64}} # table of the female growth rates associated with each active female index

    # initializes the genotypes, locations, mitochondria, and growth rates of the simulation
    function PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp)

        num_individuals_per_deme = K_total / (NUM_DEMES^2) # number of individuals per deme (innitially constant throughout range)

        intervals = collect(0.0f0:Float32(1 / NUM_DEMES):0.99f0) # locations of the demes along an axis

        deme_locations = Location.(intervals, intervals') # location of each deme (lower left corner)

        deme_populations = map(l -> l.x < 0.5 ? 0 : 1, deme_locations) # assigns 0 to each deme occupied by population A and 1 to each deme occupied by population B (dividing line down the middle at the beginning)

        # initializes the genotypes, locations, and mitochondria for each deme
        demes = Deme.(Ref(num_individuals_per_deme),
            Ref(total_loci),
            deme_locations,
            Ref(Float32(1 / NUM_DEMES)),
            deme_populations,
            Ref(ecolDiff))

        # empty matrix for storing a list of growth rates of every female per deme
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_DEMES, NUM_DEMES)

        # initializes growth rates for every female
        for deme_index in CartesianIndices((1:NUM_DEMES, 1:NUM_DEMES))
            growth_rates_F[deme_index] = calculate_growth_rates(demes,
                deme_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(demes,
            growth_rates_F)
    end

    # creates new PopulationData with the genotypes, locations, and mitochondria from the offspring
    function PopulationData(genotypes_daughters,
        genotypes_sons,
        mitochondria_daughters,
        mitochondria_sons,
        locations_daughters,
        locations_sons,
        competition_trait_loci,
        K_total,
        sigma_comp,
        intrinsic_R,
        ecolDiff)

        # set up empty matrices for storing the deme data (genotypes, locations, and mitochondria for every individual in the deme) and the female growth rates
        demes = Matrix{Deme}(undef, NUM_DEMES, NUM_DEMES)
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_DEMES, NUM_DEMES)

        # initialize each deme with the offspring data
        for deme_index in CartesianIndices((1:NUM_DEMES, 1:NUM_DEMES))
            demes[deme_index] = Deme(genotypes_daughters[deme_index],
                genotypes_sons[deme_index],
                locations_daughters[deme_index],
                locations_sons[deme_index],
                mitochondria_daughters[deme_index],
                mitochondria_sons[deme_index],
                competition_trait_loci,
                ecolDiff)
        end

        # calculate growth rates for every female in each deme
        for deme_index in CartesianIndices((1:NUM_DEMES, 1:NUM_DEMES))
            growth_rates_F[deme_index] = calculate_growth_rates(demes,
                deme_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(demes,
            growth_rates_F)

    end
end

# determine which deme a location falls in
function assign_zone(location::Location)
    zone_x = convert(Integer, trunc(NUM_DEMES * location.x) + 1)
    zone_y = convert(Integer, trunc(NUM_DEMES * location.y) + 1)
    return CartesianIndex(zone_x, zone_y)
end

# calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
function calculate_ind_useResourceA(competition_traits, ecolDiff)
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceA = competAbility .+ ((1 .- competition_traits) .* ecolDiff)

    return ind_useResourceA
end

function calculate_ind_useResourceB(competition_traits, ecolDiff)
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceB = competAbility .+ (competition_traits .* ecolDiff)

    return ind_useResourceB
end

# calculate the distance from a point along an angle to the limit of the range
function max_radius_squared(x, y, t)
    if y > 0.97 && t < pi
        y_dist = ((1 - y) / sin(t))^2
    elseif y < 0.03 && t > pi
        y_dist = (y / sin(t))^2
    else
        y_dist = 0.0009
    end

    if x > 0.97 && (t < pi / 2 || t > 3 * pi / 2)
        x_dist = ((1 - x) / cos(t))^2
    elseif x < 0.03 && (pi / 2 < t < 3 * pi / 2)
        x_dist = (x / cos(pi - t))^2
    else
        x_dist = 0.0009
    end

    return min(x_dist, y_dist)
end


# calculates the ideal density (assuming normal distribution) at the location of each female
function get_ideal_densities(K_total, sigma_comp, locations_F)
    function calc_ideal_density(location) # integral of the gaussian distribution centred on the female's location with a standard deviation of sigma_comp and cut off at the range boundaries
        if 0.03 < location.x < 0.97 && 0.03 < location.y <= 0.97 # the density is treated as constant further than 3 standard deviations from the edge of the range 
            return 1 + K_total * 2 * pi * (1 - exp(-0.00045 / (sigma_comp^2))) * (sigma_comp^2)
        else
            return 1 + K_total * (sigma_comp^2) * (2 * pi - quadgk(t -> exp(-(max_radius_squared(location.x, location.y, t) / (2 * sigma_comp^2))), 0, 2 * pi, rtol=0.04)[1])
        end
    end

    return map(calc_ideal_density, locations_F)
end

# calculates the distance between two points
function distance(location1, location2)
    return sqrt((location1.x - location2.x)^2 + (location1.y - location2.y)^2)
end

# calculates the squared distances from a set of points to a single point
function get_squared_distances(location_list, focal_location)
    dif_x = [l.x for l in location_list] .- focal_location.x
    dif_y = [l.y for l in location_list] .- focal_location.y

    return dif_x .^ 2 .+ dif_y .^ 2
end

# calculates female growth rates
function calculate_growth_rates(population,
    deme_index,
    K_total,
    sigma_comp,
    intrinsic_R)
    # in spatial model, calculate growth rates based on local resource use

    locations_F = population[deme_index].locations_F

    # determines which demes are needed to calculate growth rates
    # demes that are further than 0.03 units away from all females in the current deme are unnecessary since the individuals there would have a negligible effect
    lower_left = max(deme_index - CartesianIndex(1, 1), CartesianIndex(1, 1))
    upper_right = min(deme_index + CartesianIndex(1, 1), CartesianIndex(NUM_DEMES, NUM_DEMES))
    neighbourhood = lower_left:upper_right


    #=ind_useResourceA_all = vcat([[d.ind_useResourceA_F; d.ind_useResourceA_M] for d in neighbourhood]...)
    ind_useResourceB_all = vcat([[d.ind_useResourceB_F; d.ind_useResourceB_M] for d in neighbourhood]...)
    locations_all = vcat([[d.locations_F; d.locations_M] for d in neighbourhood]...)=#
    locations_F = population[deme_index].locations_F

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    ideal_densities_at_locations_F = get_ideal_densities(K_total, sigma_comp, locations_F) # this applies the above function to each geographic location

    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

    # assume both resources have same constant density across range

    # calculate the resource use density for each resource at a point
    function get_useResource_densities(focal_location)
        useResource_densities = [0.0, 0.0]
        for deme_index in neighbourhood
            squared_distances = get_squared_distances(population[deme_index].locations_F, focal_location)
            for i in eachindex(squared_distances)
                if squared_distances[i] <= 0.03^2
                    useResource_densities[1] += population[deme_index].ind_useResourceA_F[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
                    useResource_densities[2] += population[deme_index].ind_useResourceB_F[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
                end
            end
            squared_distances = get_squared_distances(population[deme_index].locations_M, focal_location)
            for i in eachindex(squared_distances)
                if squared_distances[i] <= 0.03^2
                    useResource_densities[1] += population[deme_index].ind_useResourceA_M[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
                    useResource_densities[2] += population[deme_index].ind_useResourceB_M[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
                end
            end
        end
        return useResource_densities
    end

    # calculates the densities for each resource at the location of each female
    real_densities = get_useResource_densities.(locations_F)

    # splits the densities up by resource
    real_densities_at_locations_F_resourceA = [d[1] for d in real_densities]
    real_densities_at_locations_F_resourceB = [d[2] for d in real_densities]

    # calculate local growth rates due to each resource (according to discrete time logistic growth equation)

    local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_locations_F_resourceA ./ (ideal_densities_at_locations_F_resourceA .+ ((real_densities_at_locations_F_resourceA) .* (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_locations_F_resourceB ./ (ideal_densities_at_locations_F_resourceB .+ ((real_densities_at_locations_F_resourceB) .* (intrinsic_R - 1)))

    growth_rateA = population[deme_index].ind_useResourceA_F .* local_growth_rates_resourceA
    growth_rateB = population[deme_index].ind_useResourceB_F .* local_growth_rates_resourceB

    growth_rates = growth_rateA .+ growth_rateB

    return growth_rates
end


# calculates the mean value of the genotype across the list of loci given
function mean(genotype, loci)
    return sum(genotype[:, loci]) / (2 * length(loci))
end

# This function calculates each mean values of the genotypes passed to it (for each individual).
# Used to determine trait values in an additive way.
# Only those loci that are additive trait loci should be passed to this function.
function calc_traits_additive(genotypes, loci)::Vector{Float32} #=::Array{Int8,3}=#
    N = length(genotypes)
    #traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64

    traits = map(x -> mean(genotypes[x], loci), 1:N)
    return traits
end

# Finds the closest male given a list of eligible males
function choose_closest_male(elig_M::Vector{Int64}, locations_M::Vector{Location}, location_mother::Location)
    focal_male = splice!(elig_M, argmin(get_squared_distances(locations_M[elig_M], location_mother))) # this gets the index of a closest male, and removes that male from the list in elig_M
    if (distance(location_mother, locations_M[focal_male]) >= 0.1)
        return focal_male, []
    end
    return focal_male, elig_M
end

# finds the closest male and updates list of eligible males
function choose_closest_male(demes::Matrix{Deme}, deme_indices::Vector{CartesianIndex{2}}, elig_M::Dict{CartesianIndex,Vector{Int64}}, location_mother::Location, neighbourhood_size::Float32)
    elig_M_per_deme = Dict{CartesianIndex,Vector{Int64}}() # initializes empty dict to store a list of remaining eligible males in each deme
    males = Dict{CartesianIndex,Int64}() # initializes empty dict to store the closest male in each deme
    for deme_index in deme_indices # loop through the demes that are nearby
        if length(elig_M[deme_index]) > 0 # checks if there are any remaining males that have not been passed over already
            deme_focal_male, deme_elig_M = choose_closest_male(elig_M[deme_index], demes[deme_index].locations_M, location_mother) # finds the closest male in that deme to the female
            if distance(demes[deme_index].locations_M[deme_focal_male], location_mother) < neighbourhood_size
                # keeps track of the index and remaining eligible males if the male is within the cutoff distance
                elig_M_per_deme[deme_index] = deme_elig_M
                males[deme_index] = deme_focal_male
            else
                # if male is beyond the cutoff distance
                males[deme_index] = -1
                elig_M_per_deme[deme_index] = []
            end
        end
    end

    # filter closest males to only keep the ones that are within the cutoff distance
    males = filter(((k, v),) -> v ≠ -1, males)

    if length(males) > 0 # checks that there are males within the cutoff distance
        # get the deme with the closest male
        index = reduce((x, y) -> distance(demes[x].locations_M[males[x]], location_mother) ≤ distance(demes[y].locations_M[males[y]], location_mother) ? x : y, keys(males))

        # update the list of eligible males from that deme to remove the closest since it has already been chosen
        output_elig_M = copy(elig_M)
        output_elig_M[index] = elig_M_per_deme[index]

        return males[index], output_elig_M, index # returns the index of the closest male, the remaining eligible males, and the index of the deme that the male was found in
    else # if there are no males within the cutoff distance returns -1
        return -1, elig_M, -1
    end

end

# compare male trait with female's trait (preference), and determine
# whether she accepts; note that match_strength is determined by a
# Gaussian, with a maximum of 1 and minimum of zero
function calc_match_strength(female_genotype, male_genotype, pref_SD, female_mating_trait_loci, male_mating_trait_loci)
    mating_trait_male = mean(male_genotype, male_mating_trait_loci)
    mating_trait_female = mean(female_genotype, female_mating_trait_loci)

    mating_trait_dif = mating_trait_male - mating_trait_female
    return exp((-(mating_trait_dif^2)) / (2 * (pref_SD^2)))
end

# generates the genotype for an offspring based on its parents' genotypes
function generate_offspring_genotype(mother_genotype, father_genotype)
    total_loci = size(mother_genotype, 2)
    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
    kid_genotype = Array{Int8,2}(undef, 2, total_loci)
    for locus in 1:total_loci
        kid_genotype[1, locus] = mother_genotype[rand([1 2]), locus]  # for this locus, pick a random allele from the mother
        kid_genotype[2, locus] = father_genotype[rand([1 2]), locus]  # and from the father
    end

    return kid_genotype
end


# creates a plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function plot_population(pd, functional_loci_range)
    genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
    locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]

    create_new_plot(calc_traits_additive(genotypes, functional_loci_range),
        [],
        locations)
end

# updates the plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function update_plot(pd, generation, functional_loci_range)
    genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
    locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]
    update_population_plot(calc_traits_additive(genotypes, functional_loci_range),
        [],
        locations,
        generation)
end
end