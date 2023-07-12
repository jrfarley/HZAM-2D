module Population

using Test
using SpecialFunctions
using QuadGK

export PopulationData, Location, Deme
export initialize_population, update_population
export plot_population, update_plot
export calc_traits_additive
export assign_zone
export NUM_DEMES
export mean

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
function max_radius_squared(x, y, t, max_dist) #x,y are the coordinates, t is the angle, and max_dist is the maximum length of the line
    if y > (1 - max_dist) && t < pi
        y_dist = ((1 - y) / sin(t))^2
    elseif y < max_dist && t > pi
        y_dist = (y / sin(t))^2
    else
        y_dist = max_dist^2
    end

    if x > (1 - max_dist) && (t < pi / 2 || t > 3 * pi / 2)
        x_dist = ((1 - x) / cos(t))^2
    elseif x < max_dist && (pi / 2 < t < 3 * pi / 2)
        x_dist = (x / cos(pi - t))^2
    else
        x_dist = max_dist^2
    end

    return min(x_dist, y_dist)
end


# calculates the ideal density (assuming normal distribution) at the location of each female
function get_ideal_densities(K_total, sigma_comp, locations_F, max_dist)
    function calc_ideal_density(location) # integral of the gaussian distribution centred on the female's location with a standard deviation of sigma_comp and cut off at the range boundaries
        if max_dist < location.x < (1-max_dist) && max_dist < location.y <= (1-max_dist) # the density is treated as constant further than 3 standard deviations from the edge of the range 
            return 1 + K_total * 2 * pi * (1 - exp(-((max_dist^2) / 2) / (sigma_comp^2))) * (sigma_comp^2)
        else
            return 1 + K_total * (sigma_comp^2) * (2 * pi - quadgk(t -> exp(-(max_radius_squared(location.x, location.y, t, max_dist) / (2 * sigma_comp^2))), 0, 2 * pi)[1])
        end
    end
    #return 1 + K_total * 2 * pi * (1 - exp(-0.00045 / (sigma_comp^2))) * (sigma_comp^2)

    return map(calc_ideal_density, locations_F)
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

    max_dist = 3 * sigma_comp #  maximum distance that the density calculation takes into account (3 standard deviations of the normal distribution)

    # in spatial model, calculate growth rates based on local resource use

    locations_F = population[deme_index].locations_F

    # determines which demes are needed to calculate growth rates
    # demes that are further than 0.03 units away from all females in the current deme are unnecessary since the individuals there would have a negligible effect
    lower_left = max(deme_index - CartesianIndex(1, 1), CartesianIndex(1, 1))
    upper_right = min(deme_index + CartesianIndex(1, 1), CartesianIndex(NUM_DEMES, NUM_DEMES))
    neighbourhood = population[lower_left:upper_right]


    ind_useResourceA_all = vcat([[d.ind_useResourceA_F; d.ind_useResourceA_M] for d in neighbourhood]...)
    ind_useResourceB_all = vcat([[d.ind_useResourceB_F; d.ind_useResourceB_M] for d in neighbourhood]...)
    locations_all = vcat([[d.locations_F; d.locations_M] for d in neighbourhood]...)
    locations_F = population[deme_index].locations_F

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    ideal_densities_at_locations_F = get_ideal_densities(K_total, sigma_comp, locations_F, max_dist) # this applies the above function to each geographic location

    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

    # assume both resources have same constant density across range

    # calculate the resource use density for each resource at a point
    function get_useResource_densities(focal_location)
        squared_distances = get_squared_distances(locations_all, focal_location)
        useResource_densities = [0.0, 0.0]
        for i in eachindex(squared_distances)
            if squared_distances[i] <= max_dist^2
                useResource_densities[1] += ind_useResourceA_all[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
                useResource_densities[2] += ind_useResourceB_all[i] * exp(-squared_distances[i] / (2 * (sigma_comp^2)))
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
end