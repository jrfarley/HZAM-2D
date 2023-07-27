"Functions and datatypes for managing population data (locations, genotypes, 
and growth rates)"
module Population
import ..DataAnalysis.calc_traits_additive
import QuadGK.quadgk

export PopulationData, Location, Zone

"The width of the (square) grid of zones that divides the population into more manageable 
chunks. Default is 10x10."
NUM_ZONES = 10

"A location with an x coordinate and a y coordinate."
struct Location
    "The location's x coordinate."
    x::Float32
    "The location's y coordinate."
    y::Float32

    """
        Location(x::Float32, y::Float32)

    Create a new location with the given x and y coordinates.
    """
    function Location(x::Float32, y::Float32)
        new(x, y)
    end

    """
        Location(starting_range::Vector{Location})

    Create a random location within the given range.

    # Arguments
    - `starting_range::Vector{Location}`: a vector of two locations specifying the lower 
    left and the upper right corners of the range.
    """
    function Location(starting_range::Vector{Location})
        x = rand() * (starting_range[2].x - starting_range[1].x) + starting_range[1].x
        y = rand() * (starting_range[2].y - starting_range[1].y) + starting_range[1].y
        new(x, y)
    end

    """
        Location(starting_location::Location, sigma_disp::Real)

    Generate an offspring's location based on the mother's location and the dispersal 
    distance. 
    
    The offspring's location is based on a normal distribution with the dispersal 
    distance being the standard deviation and the center being the given location.

    # Arguments
    - `starting_location::Location`: the center of the normal distribution.
    - `sigma_disp::Real`: the width of the normal distribution.
    """
    function Location(starting_location::Location, sigma_disp::Real)
        x = -1
        y = -1
        # Keep generating new locations until there's one that's within the range.
        while ~(0 < x < 0.999 && 0 < y < 0.999) # Check if the location is within the range.
            dist = sigma_disp * randn()
            dir = rand() * 2 * pi
            x = starting_location.x + dist * cos(dir)
            y = starting_location.y + dist * sin(dir)
        end
        new(x, y)
    end
end

"The population data (genotypes, mitochondria, locations, and resource use) for a subset of 
the total population."
struct Zone
    "The female genotypes. rows are alleles (row 1 from mother, row 2 from father) and 
    columns are loci."
    genotypes_F::Vector{Matrix{Int8}}
    "The male genotypes."
    genotypes_M::Vector{Matrix{Int8}}
    "The female locations."
    locations_F::Vector{Location}
    "The male locations."
    locations_M::Vector{Location}
    "The female mitochondria (a single locus genotype)."
    mitochondria_F::Vector{Int8}
    "The male mitochondria."
    mitochondria_M::Vector{Int8}

    "Each female's contribution to the use of resource A."
    ind_useResourceA_F::Vector{Float64}
    "Each male's contribution to the use of resource A."
    ind_useResourceA_M::Vector{Float64}
    "Each female's contribution to the use of resource B."
    ind_useResourceB_F::Vector{Float64}
    "Each male's contribution to the use of resource B."
    ind_useResourceB_M::Vector{Float64}

    """
        Zone(
            starting_N::Integer,
            total_loci::Integer,
            location::Location,
            size::Float32,
            population::Real,
            ecolDiff::Real
        )

    Set up the initial locations, genotypes, mitochondria, and resource use of every 
    individual in the zone.

    # Arguments
    - `starting_N::Integer`: the starting number of individuals in the zone.
    - `total_loci::Integer`: the number of loci in the genotypes.
    - `location::Location`: the location of the lower left corner of the zone.
    - `size::Float32`: the width of the zone.
    - `species::Real`: the 'species' identifier for the initial occupants of the zone.
    - `ecolDiff::Real`: the ecological difference between the two species.
    """
    function Zone(
        starting_N::Integer,
        total_loci::Integer,
        location::Location,
        size::Float32,
        species::Real,
        ecolDiff::Real
    )
        species = Int8(species)

        N_half = trunc(Int, starting_N / 2) # the number of individuals in each sex

        zone_range = [
            location,
            Location(min(location.x + size, 0.999f0), min(location.y + size, 0.999f0))
        ] # the geographic limits of the zone

        genotypes = fill(fill(species, 2, total_loci), N_half)

        locations = [[Location(zone_range) for i in 1:N_half] for i in 1:2]

        mitochondria = fill(species, N_half)


        #=
        Calculate individual contributions to resource use, according to linear gradient 
        between use of species 0 and species 1.

        At the beginning the mitochondria value can be used as a stand in for the 
        competition trait.
        =#
        ind_useResourceA = calc_ind_useResourceA(mitochondria, ecolDiff)
        ind_useResourceB = calc_ind_useResourceB(mitochondria, ecolDiff)

        new(
            genotypes,
            genotypes,
            locations[1],
            locations[1],
            mitochondria,
            mitochondria,
            ind_useResourceA,
            ind_useResourceA,
            ind_useResourceB,
            ind_useResourceB
        )
    end

    """
        Zone(
            genotypes_F::Vector{Matrix{Int8}},
            genotypes_M::Vector{Matrix{Int8}},
            locations_F::Vector{Location},
            locations_M::Vector{Location},
            mitochondria_F::Vector{Int8},
            mitochondria_M::Vector{Int8},
            competition_trait_loci,
            ecolDiff
        )

    Create a new Zone using the genotypes, locations, and mitochondria of the offspring.

    # Arguments
    - `genotypes_F::Vector{Matrix{Int8}}`: the female genotypes.
    - `genotypes_M::Vector{Matrix{Int8}}`: the male genotypes.
    - `locations_F::Vector{Location}`: the female locations.
    - `locations_M::Vector{Location}`: the male locations.
    - `mitochondria_F::Vector{Int8}`: the female mitochondria.
    - `mitochondria_M::Vector{Int8}`: the male mitochondria.
    - `competition_trait_loci`: the loci specifying the ecological trait (used in 
    fitness related to resource use).
    - `ecolDiff::Real`: the ecological difference between the two species.
    """
    function Zone(
        genotypes_F::Vector{Matrix{Int8}},
        genotypes_M::Vector{Matrix{Int8}},
        locations_F::Vector{Location},
        locations_M::Vector{Location},
        mitochondria_F::Vector{Int8},
        mitochondria_M::Vector{Int8},
        competition_trait_loci,
        ecolDiff::Real
    )
        # calculate new competition traits
        competition_traits_F = calc_traits_additive(genotypes_F, competition_trait_loci)
        competition_traits_M = calc_traits_additive(genotypes_M, competition_trait_loci)

        #=
        calculate individual contributions to resource use, according to linear gradient 
        between use of species 0 and species 1
        =#
        ind_useResourceA_F = calc_ind_useResourceA(competition_traits_F, ecolDiff)
        ind_useResourceA_M = calc_ind_useResourceA(competition_traits_M, ecolDiff)
        ind_useResourceB_F = calc_ind_useResourceB(competition_traits_F, ecolDiff)
        ind_useResourceB_M = calc_ind_useResourceB(competition_traits_M, ecolDiff)

        new(
            genotypes_F,
            genotypes_M,
            locations_F,
            locations_M,
            mitochondria_F,
            mitochondria_M,
            ind_useResourceA_F,
            ind_useResourceA_M,
            ind_useResourceB_F,
            ind_useResourceB_M
        )
    end
end


"The population data (genotypes, locations, mitochondria, resource use, and growth rates) 
for all individuals in the simulation. The data is subdivided by 'zone' into a matrix for 
calculation efficiency."
struct PopulationData
    "All of the zones (containing the genotypes, locations, mitochondria, and resource use) 
    in the simulation."
    population::Matrix{Zone}

    "A list of the growth rates of every female for each zone."
    growth_rates_F::Matrix{Vector{Float64}}

    """
        calc_zone_population(K_total::Integer, ecolDiff::Real)

    Compute the initial number of individuals in each zone. 
    
    When ecolDiff=0 this will just 
    be the carrying capacity (K_total) divided by the number of zones. But when there is an 
    ecological difference between the two species the effective carrying capacity is reduced 
    since not all of the available resources will be used.
    """
    function calc_zone_population(K_total::Integer, ecolDiff::Real)
        effective_K = K_total / (1 + ecolDiff)

        floor(Int, effective_K / (NUM_ZONES^2))
    end

    # initializes the genotypes, locations, mitochondria, and growth rates of the simulation
    """
        PopulationData(
            K_total::Integer,
            ecolDiff::Real,
            total_loci::Integer,
            intrinsic_R::Real,
            sigma_comp::Real
        )

    Set up the initial population for the simulation.

    # Arguments
    - `K_total::Integer`: the carrying capacity.
    - `ecolDiff::Real`: the ecological difference between the two species.
    - `total_loci::Integer`: the total number of loci for the genotypes.
    - `intrinsic_R::Real`: the intrinsic growth rate.
    - `sigma_comp::Real`: the standard deviation for the normal curve used in calculating 
    local density.
    """
    function PopulationData(
        K_total::Integer,
        ecolDiff::Real,
        total_loci::Integer,
        intrinsic_R::Real,
        sigma_comp::Real
    )

        # the number of individuals per zone (innitially constant throughout range)
        num_individuals_per_zone = calc_zone_population(K_total, ecolDiff)

        # the locations of the zones along an axis
        spaced_locations = collect(0.0f0:Float32(1 / NUM_ZONES):0.99f0)

        # the location of each zone (lower left corner)
        zone_locations = Location.(spaced_locations, spaced_locations')

        #= 
        Assign each zone the species identifier of the species that will initially
        occupy it. The dividing line between species will initially be x=0.5.
        =#
        zone_species = map(l -> l.x < 0.5 ? 0 : 1, zone_locations)

        # initialize the genotypes, locations, and mitochondria for each zone
        zones = Zone.(
            Ref(num_individuals_per_zone),
            Ref(total_loci),
            zone_locations,
            Ref(Float32(1 / NUM_ZONES)),
            zone_species,
            Ref(ecolDiff)
        )

        # empty matrix for storing a list of growth rates of every female per zone
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_ZONES, NUM_ZONES)

        # initialize the growth rates for every female
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            growth_rates_F[zone_index] = calc_growth_rates(zones,
                zone_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(zones, growth_rates_F)
    end

    """
        PopulationData(
            genotypes_daughters::Matrix{Vector{Matrix{Int8}}},
            genotypes_sons::Matrix{Vector{Matrix{Int8}}},
            mitochondria_daughters::Matrix{Vector{Int8}},
            mitochondria_sons::Matrix{Vector{Int8}},
            locations_daughters::Matrix{Vector{Location}},
            locations_sons::Matrix{Vector{Location}},
            competition_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
            K_total::Integer,
            sigma_comp::Real,
            intrinsic_R::Real,
            ecolDiff::Real
        )

    Create a new PopulationData using the genotypes, locations, and mitochondria of the 
    offspring.

    # Arguments
    - `genotypes_daughters::Matrix{Vector{Matrix{Int8}}}`: the female genotypes.
    - `genotypes_sons::Matrix{Vector{Matrix{Int8}}}`: the male genotypes.
    - `mitochondria_daughters::Matrix{Vector{Int8}}`: the female mitochondria.
    - `mitochondria_sons::Matrix{Vector{Int8}}`: the male mitochondria.
    - `locations_daughters::Matrix{Vector{Location}}`: the female locations.
    - `locations_sons::Matrix{Vector{Location}}`: the male locations.
    - `competition_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: list of the 
    loci specifying the ecological trait (used in fitness related to resource use).
    - `K_total::Integer`: the carrying capacity of the environment
    - `sigma_comp::Real`: the standard deviation for the normal curve used in calculating 
    local density.
    - `intrinsic_R::Real`: the intrinsic growth rate.
    - `ecolDiff::Real`: the ecological difference between the two species.
    """
    function PopulationData(
        genotypes_daughters::Matrix{Vector{Matrix{Int8}}},
        genotypes_sons::Matrix{Vector{Matrix{Int8}}},
        mitochondria_daughters::Matrix{Vector{Int8}},
        mitochondria_sons::Matrix{Vector{Int8}},
        locations_daughters::Matrix{Vector{Location}},
        locations_sons::Matrix{Vector{Location}},
        competition_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
        K_total::Integer,
        sigma_comp::Real,
        intrinsic_R::Real,
        ecolDiff::Real
    )

        #= 
        Set up empty matrices for storing the zone data (genotypes, locations, and 
        mitochondria for every individual in the zone) and the female growth rates.
        =#
        zones = Matrix{Zone}(undef, NUM_ZONES, NUM_ZONES)
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_ZONES, NUM_ZONES)

        # initialize each zone with the offspring data
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            zones[zone_index] = Zone(genotypes_daughters[zone_index],
                genotypes_sons[zone_index],
                locations_daughters[zone_index],
                locations_sons[zone_index],
                mitochondria_daughters[zone_index],
                mitochondria_sons[zone_index],
                competition_trait_loci,
                ecolDiff)
        end

        # calculate growth rates for every female in each zone
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            growth_rates_F[zone_index] = calc_growth_rates(zones,
                zone_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(zones, growth_rates_F)

    end
end

"""
    assign_zone(location::Location)

Determine which zone a location falls in.
"""
function assign_zone(location::Location)
    zone_x = convert(Integer, trunc(NUM_ZONES * location.x) + 1) # x index of the zone
    zone_y = convert(Integer, trunc(NUM_ZONES * location.y) + 1) # y index of the zone
    return CartesianIndex(zone_x, zone_y)
end

"""
    calc_ind_useResourceA(competition_traits::Vector{Real}, ecolDiff::Real)

Compute the individual contributions to the use of resource A according to a linear gradient 
between the resource use of species 0 and species 1.

# Arguments
- `competition_traits::Vector{Real}`: the ecological competition trait values
- `ecolDiff::Real`: the ecological difference between the species. 
"""
function calc_ind_useResourceA(competition_traits::Vector, ecolDiff::Real)
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceA = competAbility .+ ((1 .- competition_traits) .* ecolDiff)

    return ind_useResourceA
end


"""
    calc_ind_useResourceB(competition_traits::Vector{Real}, ecolDiff::Real)

Compute the individual contributions to the use of resource B according to a linear gradient 
between the resource use of species 0 and species 1.

# Arguments
- `competition_traits::Vector{Real}`: the ecological competition trait values
- `ecolDiff::Real`: the ecological difference between the species. 
"""
function calc_ind_useResourceB(competition_traits::Vector, ecolDiff::Real)
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceB = competAbility .+ (competition_traits .* ecolDiff)

    return ind_useResourceB
end

"""
    max_radius_squared(x::Real, y::Real, t::Real, max_dist::Real)

Compute the distance from a point along an angle to the limit of the range 
(defined by x∈[0,1), y∈[0,1)). 

The distance gets cut off at max_dist. Used in calculating 
the ideal densities.

# Arguments
- `location::Location`: the starting location.
- `t::Real`: the angle.
- `max_dist`: the maximum distance.
"""
function max_radius_squared(location::Location, t::Real, max_dist::Real)
    x = location.x
    y = location.y

    if y > (1 - max_dist) && t < pi
        y_dist = (1 - y) / sin(t)
    elseif y < max_dist && t > pi
        y_dist = y / sin(t)
    else
        y_dist = max_dist
    end

    if x > (1 - max_dist) && (t < pi / 2 || t > 3 * pi / 2)
        x_dist = (1 - x) / cos(t)
    elseif x < max_dist && (pi / 2 < t < 3 * pi / 2)
        x_dist = x / cos(pi - t)
    else
        x_dist = max_dist
    end

    return (min(x_dist, y_dist))^2
end

"""
    calc_ideal_densities(
        K_total::Integer,
        sigma_comp::Real,
        locations_F::Vector{Location},
        max_dist::Real
    )

Compute expected local densities, based on geographically even distribution of individuals 
at carrying capacity.

# Arguments
- `K_total::Integer`: the carrying capacity.
- `sigma_comp::Real`: the standard deviation for the normal curve used in calculating 
local density.
- `locations_F::Vector{Location}`: a list of the locations where the ideal density is to be 
computed.
- `max_dist::Real`: the cutoff for the furthest away an individual can be and affect the 
density calculation. Should be 3x sigma_comp.
"""
function calc_ideal_densities(
    K_total::Integer,
    sigma_comp::Real,
    locations_F::Vector{Location},
    max_dist::Real
)

    """
    Compute densities using the equation:

    ``K_{total}\\int\\limits_0^{2\\pi}\\int\\limits_0^{0.03} 
    \\exp(-\\frac{r^2}{2σ_{comp}^2})rdrdθ``

    where the density function is equal to 0 outside of the range limits.
    """
    function calc_ideal_density(location)
        if max_dist < location.x < (1 - max_dist) &&
           max_dist < location.y <= (1 - max_dist)
            # the ideal density is constant further than 3 sigma_comps from the range limits 
            return begin
                1 +
                K_total * 2 * pi *
                (1 - exp(-((max_dist^2) / 2) / (sigma_comp^2))) * (sigma_comp^2)
            end
        else
            # density function
            f(t) = exp(-(max_radius_squared(location, t, max_dist) / (2 * sigma_comp^2)))
            return begin
                1 +
                K_total * (sigma_comp^2) *
                (2 * pi -
                 quadgk(f, 0, 2 * pi)[1])
            end
        end
    end

    return map(calc_ideal_density, locations_F)
end

"""
    calc_squared_distances(location_list::Vector{Location}, focal_location::Location)

Compute the squared distances from a set of locations to a single location.

# Example
```jldoctest
julia> calc_squared_distances([Location(0.5f0, 0.5f0), Location(0.4f0, 0.4f0), 
Location(0f0,0f0)], Location(0.5f0, 0.5f0))
3-element Vector{Float32}:
 0.0
 0.019999998
 0.5
 ```
"""
function calc_squared_distances(location_list::Vector{Location}, focal_location::Location)
    dif_x = [l.x for l in location_list] .- focal_location.x
    dif_y = [l.y for l in location_list] .- focal_location.y

    return dif_x .^ 2 .+ dif_y .^ 2
end

"""
    calc_growth_rates(
        population::Matrix{Zone},
        zone_index::CartesianIndex,
        K_total::Integer,
        sigma_comp::Real,
        intrinsic_R::Real
    )

Compute the female growth rates.

# Arguments
- `population::Matrix{Zone}`: the matrix of zones storing all the locations, and genotypes 
in the simulation.
- `zone_index::CartesianIndex`: the index of the focal zone.
- `K_total::Integer`: the carrying capacity.
- `sigma_comp::Real`: the standard deviation for the normal curve used in calculating 
local density.
- `intrinsic_R::Real`: the intrinsic growth rate.
"""
function calc_growth_rates(
    population::Matrix{Zone},
    zone_index::CartesianIndex,
    K_total::Integer,
    sigma_comp::Real,
    intrinsic_R::Real
)
    # maximum distance that the density calculation takes into account
    max_dist = 3 * sigma_comp

    locations_F = population[zone_index].locations_F

    #= 
    Determine which zones are needed to calculate growth rates. Zones that are further than 
    0.03 units away from all females in the current zone are not used in the growth rate 
    calculations.
    =#
    lower_left = max(zone_index - CartesianIndex(1, 1), CartesianIndex(1, 1))
    upper_right = min(
        zone_index + CartesianIndex(1, 1),
        CartesianIndex(NUM_ZONES, NUM_ZONES)
    )
    neighbourhood = population[lower_left:upper_right]


    # calculate the individual contributions to the use of both resources.
    ind_useResourceA_all = vcat(
        [[z.ind_useResourceA_F; z.ind_useResourceA_M] for z in neighbourhood]...
    )
    ind_useResourceB_all = vcat(
        [[z.ind_useResourceB_F; z.ind_useResourceB_M] for z in neighbourhood]...
    )
    locations_all = vcat([[z.locations_F; z.locations_M] for z in neighbourhood]...)
    locations_F = population[zone_index].locations_F

    #= 
    Set up the expected local densities, based on a geographically even distribution of 
    individuals at the carrying capacity.
    =#
    ideal_densities_at_locations_F = calc_ideal_densities(
        K_total,
        sigma_comp,
        locations_F,
        max_dist
    )

    # assume both resources have same constant density across range
    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

    """
    Calculate the resource use density for both resources at the given location.
    """
    function calc_useResource_densities(focal_location)
        squared_distances = calc_squared_distances(locations_all, focal_location)
        useResource_densities = [0.0, 0.0]
        for i in eachindex(squared_distances)
            if squared_distances[i] <= max_dist^2
                useResource_densities[1] += ind_useResourceA_all[i] *
                                            exp(-squared_distances[i] /
                                                (2 * (sigma_comp^2)))
                useResource_densities[2] += ind_useResourceB_all[i] *
                                            exp(-squared_distances[i] /
                                                (2 * (sigma_comp^2)))
            end
        end
        return useResource_densities
    end

    # calculate the densities for both resources at the location of each female
    real_densities = calc_useResource_densities.(locations_F)

    # split the densities up by resource
    real_densities_at_locations_F_resourceA = [d[1] for d in real_densities]
    real_densities_at_locations_F_resourceB = [d[2] for d in real_densities]

    #= 
    Calculate the local growth rates due to each resource according to discrete time 
    logistic growth equation.
    =#
    local_growth_rates_resourceA = intrinsic_R .*
                                   ideal_densities_at_locations_F_resourceA ./
                                   (ideal_densities_at_locations_F_resourceA .+
                                    ((real_densities_at_locations_F_resourceA) .*
                                     (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .*
                                   ideal_densities_at_locations_F_resourceB ./
                                   (ideal_densities_at_locations_F_resourceB .+
                                    ((real_densities_at_locations_F_resourceB) .*
                                     (intrinsic_R - 1)))

    growth_rateA = population[zone_index].ind_useResourceA_F .* local_growth_rates_resourceA
    growth_rateB = population[zone_index].ind_useResourceB_F .* local_growth_rates_resourceB

    growth_rates = growth_rateA .+ growth_rateB

    return growth_rates
end

"""
    genotype_mean(genotype::Matrix{<:Integer}, loci::Union{UnitRange{<:Integer},Vector{<:Integer}})

Compute the mean value of the genotype across the list of loci given.

# Example 
```jldoctest
julia> genotype_mean([0 1 0; 0 1 0], [1,3])
0.0
```
"""
function genotype_mean(
    genotype::Matrix{<:Integer},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    return sum(genotype[:, loci]) / (2 * length(loci))
end

"""
    calc_survival_fitness_hetdisadvantage(genotype::Matrix{<:Integer}, w_hyb::Real)

Compute the survival fitness of an individual according to heterozygosity.

# Arguments
- `genotype:Matrix{<:Integer}`: the individual's genotype.
- `w_hyb::Real`: the simulation's hybrid fitness value.
"""
function calc_survival_fitness_hetdisadvantage(
    genotype::Matrix{<:Integer},
    w_hyb::Real
)::Vector{Float32}
    num_loci = size(genotype, 2)
    s_per_locus = 1 - w_hyb^(1 / num_loci)  # loss in fitness due to each heterozygous locus 
    num_hetloci = sum(genotype[1, :] .≠ genotype[2, :])

    hetdisadvantage_fitness = (1 - s_per_locus)^num_hetloci
    return hetdisadvantage_fitness
end

"""
    calc_survival_fitness_epistasis(
        genotype::Matrix{<:Integer}, 
        hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}, 
        w_hyb::Real, 
        beta::Real=1
    )

Compute the survival fitness of an individual according to epistasis, with the beta 
parameter set to one as a default.

# Arguments
- `genotype:Matrix{<:Integer}`: the individual's genotype.
- `hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci 
responsible for heterozygote disadvantage in computing the survival fitness.
- `w_hyb::Real`: the simulation's hybrid fitness value.
- `beta::Real`
"""
function calc_survival_fitness_epistasis(genotype::Matrix{<:Integer},
    hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
    w_hyb::Real,
    beta::Real=1
)::Float32
    survival_HI = genotype_mean(genotype, hybrid_survival_loci)
    epistasis_fitnesses = 1 - (1 - w_hyb) * (4 * survival_HI * (1 - survival_HI))^beta
    return epistasis_fitnesses
end
end