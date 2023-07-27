"Functions involved in mate selection and generating the offspring genotype."
module Mating
using ..Population

"""
    distance(location1::Location, location2::Location)

Compute the distance between `location1` and `location2`.
"""
function distance(location1::Location, location2::Location)
    return sqrt((location1.x - location2.x)^2 + (location1.y - location2.y)^2)
end

"""
    choose_closest_male(
        zones::Matrix{Zone}, 
        zone_indices::Vector{CartesianIndex{2}}, 
        elig_M::Dict{CartesianIndex,<:Vector{<:Integer}}, 
        location_mother::Location, 
        neighbourhood_size::Float32
    )

Find the index of the closest male within an area covering multiple zones from a list of 
eligible males. 

Return both the male index and its zone index. Return -1 for both values if 
no eligible male is found.

# Arguments
- `zones::Matrix{Zone}`: the matrix storing the population data for each zone.
- `zone_indices::Vector{CartesianIndex{2}}`: the indices of the zones to be checked for the 
closest male.
- `elig_M::Dict{CartesianIndex,<:Vector{<:Integer}}`: the indices of all the eligible males in 
each zone.
- `location_mother::Location`: the location from which the males' distances are computed.
- `neighbourhood_size::Float32`: the distance cutoff for the search.
"""
function choose_closest_male(
    zones::Matrix{Zone},
    zone_indices::Vector{CartesianIndex{2}},
    elig_M::Dict{CartesianIndex,<:Vector{<:Integer}},
    location_mother::Location,
    neighbourhood_size::Float32
)
    shortest_distance = neighbourhood_size # males further than this distance are ignored
    output_male = -1
    output_zone = -1
    for zone_index in zone_indices # loop through the zones that are nearby

        # check if there are any remaining males that have not been passed over already
        if length(elig_M[zone_index]) == 0
            continue
        end

        # find the closest male in that zone to the female
        male_index = choose_closest_male_from_zone(
            elig_M[zone_index],
            zones[zone_index].locations_M,
            location_mother
        )

        # calculate distance from the mother's location to the closest male
        male_distance = distance(zones[zone_index].locations_M[male_index], location_mother)

        # check if the distance is smaller than the previous closest distance
        if male_distance < shortest_distance
            output_male = male_index
            output_zone = zone_index
            shortest_distance = male_distance
        end
    end

    return output_male, output_zone
end

"""
    choose_closest_male_from_zone(
        elig_M::Vector{<:Integer}, 
        locations_M::Vector{Location}, 
        location_mother::Location
    )

Find the index of the closest male within a single zone from a list of eligible males.

# Arguments
- `elig_M::Vector{<:Integer}`: the indices of all the eligible males.
- `locations_M::Vector{Location}`: the locations of all the males in the zone.
- `location_mother::Location`: the location from which the males' distances are computed.
"""
function choose_closest_male_from_zone(
    elig_M::Vector{<:Integer},
    locations_M::Vector{Location},
    location_mother::Location)
    return elig_M[argmin(
        Population.calc_squared_distances(
            locations_M[elig_M],
            location_mother
        )
    )]
end

"""
    calc_match_strength(
        female_genotype::Matrix{<:Integer}, 
        male_genotype::Matrix{<:Integer}, 
        pref_SD::Float64, 
        female_mating_trait_loci, 
        male_mating_trait_loci
    )

Compare the male's mating trait with the female mating trait to determine the match strength 
along a Gaussian acceptance curve. 

The acceptance curve is bounded between 0 and 1.

# Arguments
- `female_genotype::Matrix{<:Integer}`: the genotype of the focal female
- `male_genotype::Matrix{<:Integer}`: the genotype of the focal male
- `pref_SD::Float64`: width of the Gaussian acceptance curve
- `female_mating_trait_loci`: the loci contributing to the female's mate preference
- `male_mating_trait_loci`: the loci contributing to the male's trait
"""
function calc_match_strength(
    female_genotype::Matrix{<:Integer},
    male_genotype::Matrix{<:Integer},
    pref_SD::Float64,
    female_mating_trait_loci,
    male_mating_trait_loci
)
    # calculate the mean value of the male and female mating traits
    mating_trait_M = Population.genotype_mean(male_genotype, male_mating_trait_loci)
    mating_trait_F = Population.genotype_mean(female_genotype, female_mating_trait_loci)

    mating_trait_dif = mating_trait_M - mating_trait_F
    return exp((-(mating_trait_dif^2)) / (2 * (pref_SD^2)))
end

"""
    generate_offspring_genotype(
        mother_genotype::Matrix{<:Integer}, 
        father_genotype::Matrix{<:Integer}
    )::Matrix{<:Integer}

Generate the offspring genotype from the parent genotypes. 

For each locus (column), first row for allele from mother, second row for allele from 
father. The number of loci in each genotype must be equal.
    
# Example
```jldoctest
julia> generate_offspring_genotype([0 0 0; 0 0 0], [1 1 1; 1 1 1])
2×3 Matrix:
 0  0  0
 1  1  1
```
"""
function generate_offspring_genotype(
    mother_genotype::Matrix{<:Integer},
    father_genotype::Matrix{<:Integer}
)::Matrix{<:Integer}
    total_loci = size(mother_genotype, 2)
    kid_genotype = Matrix{Integer}(undef, 2, total_loci)
    for locus in 1:total_loci
        # pick a random allele from the mother
        kid_genotype[1, locus] = mother_genotype[rand([1 2]), locus]
        kid_genotype[2, locus] = father_genotype[rand([1 2]), locus]  # and from the father
    end

    return kid_genotype
end
end