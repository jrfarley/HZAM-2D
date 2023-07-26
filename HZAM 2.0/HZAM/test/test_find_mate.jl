using Random

@testset "choose_closest_male_no_zones" begin
    elig_M = collect(1:1001)
    location_mother = Location(0.5f0, 0.5f0)

    locations_M = Location.(Float32.(rand(1000)), Float32.(rand(1000)))
    male_loc = Location(0.5f0, 0.5f0)
    locations_M = shuffle(push!(locations_M, male_loc))

    closest_male = Mating.choose_closest_male_from_zone(elig_M, locations_M, location_mother)

    @test locations_M[closest_male] == male_loc

    for i in 1:100
        male_loc = Location(0.55f0, 0.45f0)
        locations_M = Location.(Float32.(rand(1000) .* Ref(0.4)), Float32.(rand(1000) .* Ref(0.4) .+ Ref(0.55)))
        locations_M = shuffle(push!(locations_M, male_loc))

        closest_male = Mating.choose_closest_male_from_zone(elig_M, locations_M, location_mother)

        @test locations_M[closest_male] == male_loc
    end
end

@testset "choose_closest_male_one_zone" begin
    intervals = collect(0.0f0:0.1f0:0.99f0) # locations of the zones along an axis

    zone_locations = Location.(intervals, intervals') # location of each zone (lower left corner)

    zone_populations = map(l -> l.x < 0.5 ? 1 : 1, zone_locations) # assigns 0 to each zone occupied by population A and 1 to each zone occupied by population B (dividing line down the middle at the beginning)

    # initializes the genotypes, locations, and mitochondria for each zone
    zones = Population.Zone.(Ref(1000),
        Ref(10),
        zone_locations,
        Ref(0.1f0),
        zone_populations,
        Ref(0))

    location_mother = zones[5,5].locations_F[1]

    elig_M = Dict{CartesianIndex,Vector{Int64}}()

    for i in eachindex(IndexCartesian(), zones)
        elig_M[i] = 1:length(zones[i].locations_M)
    end

    closest_male, zone_index = choose_closest_male(zones, [CartesianIndex(5, 5)], elig_M, location_mother, 0.1f0)

    @test zones[zone_index].locations_M[closest_male] == location_mother
    @test zone_index == CartesianIndex(5, 5)
end

@testset "choose_closest_male_multiple_zones" begin
    intervals = collect(0.0f0:0.1f0:0.99f0) # locations of the zones along an axis

    zone_locations = Location.(intervals, intervals') # location of each zone (lower left corner)

    zone_populations = map(l -> l.x < 0.5 ? 1 : 1, zone_locations) # assigns 0 to each zone occupied by population A and 1 to each zone occupied by population B (dividing line down the middle at the beginning)

    for i in 1:10
        # initializes the genotypes, locations, and mitochondria for each zone
        zones = Population.Zone.(Ref(10),
            Ref(10),
            zone_locations,
            Ref(0.1f0),
            zone_populations,
            Ref(0))

        location_mother = Location(0.45f0, 0.45f0)

        locations_M = vcat(zones[i, i].locations_M, location_mother)

        genotypes_F = zones[i, i].genotypes_F
        genotypes_M = zones[i, i].genotypes_M
        push!(genotypes_M, genotypes_M[1])
        locations_F = zones[i, i].locations_F
        mitochondria_F = zones[i, i].mitochondria_F
        mitochondria_M = vcat(zones[i, i].mitochondria_M, Int8(1))

        zones[i, i] = Population.Zone(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M, collect(1:2), 0)

        elig_M = Dict{CartesianIndex,Vector{Int64}}()

        for i in eachindex(IndexCartesian(), zones)
            elig_M[i] = 1:length(zones[i].locations_M)
        end
        neighbourhood = filter(e -> length(elig_M[e]) > 0, CartesianIndex(1, 1):CartesianIndex(10, 10))

        closest_male, zone_index = choose_closest_male(zones, neighbourhood, elig_M, location_mother, 1.0f0)

        @test zones[zone_index].locations_M[closest_male] == location_mother
        @test zone_index == CartesianIndex(i, i)
    end
end



@testset "generate_offspring_genotype" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    @test generate_offspring_genotype(genotypes[1], genotypes[2]) == [0 0 0; 1 1 1]
    @test generate_offspring_genotype(genotypes[1], genotypes[1]) == [0 0 0; 0 0 0]
end


@testset "calc_match_strength" begin
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 3:3
    pref_SD = sqrt(-1 / (2 * log(1 / (1 + 10^(-15)))))

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    match_strength1 = calc_match_strength(genotypes[1], genotypes[1], pref_SD, female_mating_trait_loci, male_mating_trait_loci)
    match_strength2 = calc_match_strength(genotypes[2], genotypes[2], pref_SD, female_mating_trait_loci, male_mating_trait_loci)
    match_strength3 = calc_match_strength(genotypes[1], genotypes[2], pref_SD, female_mating_trait_loci, male_mating_trait_loci)

    @test match_strength1 == match_strength2
    @test match_strength1 > match_strength3
end