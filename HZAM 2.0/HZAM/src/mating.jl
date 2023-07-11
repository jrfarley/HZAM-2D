module Mating
using ..Population
export choose_closest_male, calc_match_strength, generate_offspring_genotype,a

# calculates the distance between two points
function distance(location1, location2)
    return sqrt((location1.x - location2.x)^2 + (location1.y - location2.y)^2)
end

# finds the closest male and updates list of eligible males
function choose_closest_male(demes::Matrix{Deme}, deme_indices::Vector{CartesianIndex{2}}, elig_M::Dict{CartesianIndex,Vector{Int64}}, location_mother::Location, neighbourhood_size::Float32)

    shortest_distance = neighbourhood_size # males further than this distance are ignored
    output_male = -1
    output_deme = -1
    for deme_index in deme_indices # loop through the demes that are nearby
        if length(elig_M[deme_index]) == 0 # checks if there are any remaining males that have not been passed over already
            continue
        end

        male_index = choose_closest_male_from_deme(elig_M[deme_index], demes[deme_index].locations_M, location_mother) # finds the closest male in that deme to the female
        male_distance = distance(demes[deme_index].locations_M[male_index], location_mother) # calculates distance from the mother's location to the closest male
        if male_distance < shortest_distance # checks if the distance is smaller than the previous closest distance
            output_male = male_index
            output_deme = deme_index
            shortest_distance = male_distance
        end
        # keeps track of the index and remaining eligible males if the male is within the cutoff distance
    end

    return output_male, output_deme
end

# Finds the closest male given a list of eligible males
function choose_closest_male_from_deme(elig_M::Vector{Int64}, locations_M::Vector{Location}, location_mother::Location)
    return elig_M[argmin(Population.get_squared_distances(locations_M[elig_M], location_mother))] # this gets the index of a closest male, and removes that male from the list in elig_M
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
end