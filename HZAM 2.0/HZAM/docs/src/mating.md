# Mating
Functions involved in mate selection and generating the offspring genotype.

## Functions
```@docs
HZAM.Mating.distance(location1::HZAM.Population.Location, location2::HZAM.Population.Location)
HZAM.Mating.choose_closest_male(
        zones::Matrix{HZAM.Population.Zone}, 
        zone_indices::Vector{CartesianIndex{2}}, 
        elig_M::Dict{CartesianIndex,<:Vector{<:Integer}}, 
        location_mother::HZAM.Population.Location, 
        neighbourhood_size::Float32
    )
HZAM.Mating.choose_closest_male_from_zone(
        elig_M::Vector{<:Integer}, 
        locations_M::Vector{HZAM.Population.Location}, 
        location_mother::HZAM.Population.Location
    )
HZAM.Mating.calc_match_strength(
        female_genotype::Matrix{<:Integer}, 
        male_genotype::Matrix{<:Integer}, 
        pref_SD::Float64, 
        female_mating_trait_loci, 
        male_mating_trait_loci
    )
HZAM.Mating.generate_offspring_genotype(
        mother_genotype::Matrix{<:Integer}, 
        father_genotype::Matrix{<:Integer}
    )
```