using Random

@testset "choose_closest_male_no_zones" begin
    elig_M = collect(1:1001)
    location_mother = Location(0.5f0, 0.5f0)

    locations_M = Location.(Float32.(rand(1000)), Float32.(rand(1000)))
    male_loc = Location(0.5f0, 0.5f0)
    locations_M = shuffle(push!(locations_M, male_loc))

    closest_male =
        Mating.choose_closest_male_from_zone(
            elig_M,
            locations_M,
            location_mother
        )

    @test locations_M[closest_male] == male_loc

    for i in 1:100
        male_loc = Location(0.55f0, 0.45f0)
        locations_M =
            Location.(
                Float32.(rand(1000) .* Ref(0.4)),
                Float32.(rand(1000) .* Ref(0.4) .+ Ref(0.55))
            )

        locations_M = shuffle(push!(locations_M, male_loc))

        closest_male =
            Mating.choose_closest_male_from_zone(
                elig_M,
                locations_M,
                location_mother
            )

        @test locations_M[closest_male] == male_loc
    end
end

@testset "choose_closest_male_one_zone" begin
    # locations of the zones along an axis
    intervals = collect(0.0f0:0.1f0:0.99f0)

    # location of each zone (lower left corner)
    zone_locations = Location.(intervals, intervals')

    # assign the population number to each zone
    zone_populations = map(l -> l.x < 0.5 ? 1 : 1, zone_locations)

    # initialize the genotypes and locations and for each zone
    zones = Zone.(Ref(1000),
        Ref(10),
        zone_locations,
        Ref(0.1f0),
        zone_populations,
        Ref(0))

    location_mother = zones[5, 5].locations_F[1]

    elig_M = Dict{CartesianIndex,Vector{Int64}}()

    for i in eachindex(IndexCartesian(), zones)
        elig_M[i] = 1:length(zones[i].locations_M)
    end

    closest_male, zone_index =
        Mating.choose_closest_male(
            zones, [CartesianIndex(5, 5)],
            elig_M, location_mother, 0.1f0
        )

    @test zones[zone_index].locations_M[closest_male] == location_mother
    @test zone_index == CartesianIndex(5, 5)
end

@testset "choose_closest_male_multiple_zones" begin
    # locations of the zones along an axis
    intervals = collect(0.0f0:0.1f0:0.99f0)

    # location of each zone (lower left corner)
    zone_locations = Location.(intervals, intervals')

    # assign the population number to each zone
    zone_populations = map(l -> l.x < 0.5 ? 1 : 1, zone_locations)

    for i in 1:10
        # initializes the genotypes and locations for each zone
        zones = Zone.(Ref(10),
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

        zones[i, i] = Zone(
            genotypes_F,
            genotypes_M,
            locations_F,
            locations_M,
            collect(1:2),
            0
        )

        elig_M = Dict{CartesianIndex,Vector{Int64}}()

        for i in eachindex(IndexCartesian(), zones)
            elig_M[i] = 1:length(zones[i].locations_M)
        end
        neighbourhood =
            filter(e -> length(elig_M[e]) > 0, CartesianIndex(1, 1):CartesianIndex(10, 10))

        closest_male, zone_index =
            Mating.choose_closest_male(
                zones,
                neighbourhood,
                elig_M,
                location_mother,
                1.0f0
            )

        @test zones[zone_index].locations_M[closest_male] == location_mother
        @test zone_index == CartesianIndex(i, i)
    end
end

@testset "choose_closest_male_empty_zone" begin
    zone = Zone(0, 0, Location(0.0f0, 0.0f0), 0.1f0, 0, 0)
    zones = [zone zone; zone zone]
    elig_M = Dict{CartesianIndex,Vector{Int64}}()
    elig_M[CartesianIndex(1, 1)] = []

    closest_male = Mating.choose_closest_male(
        zones,
        [CartesianIndex(1, 1)],
        elig_M,
        Location(0.5f0, 0.6f0),
        0.5f0
    )
end



@testset "generate_offspring_genotype" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 0 1; 1 0 1]]

    @test Mating.generate_offspring_genotype(genotypes[1], genotypes[2]) == [0 0 0; 1 1 1]
    @test Mating.generate_offspring_genotype(genotypes[1], genotypes[1]) == [0 0 0; 0 0 0]
    @test Mating.generate_offspring_genotype(genotypes[3], genotypes[1]) == [1 0 1; 0 0 0]
end


@testset "calc_match_strength_S_AM1" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2
    pref_SD = sqrt(-1 / (2 * log(1 / (1 + 10^(-15)))))

    pref_SD = Inf

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 1
    @test match_strength12 ≈ 1
    @test match_strength13 ≈ 1
end

@testset "calc_match_strength_S_AM_inf" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2

    pref_SD = 10^(-15)

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 0
    @test match_strength12 ≈ 0
    @test match_strength13 ≈ 0
end

@testset "calc_match_strength_S_AM_10" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2

    pref_SD = sqrt(-1 / (2 * log(1 / 10)))

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 0.5623413251903492
    @test match_strength12 ≈ 0.1
    @test match_strength13 ≈ 0.5623413251903492
end

@testset "distance" begin
    loc1 = Location(0.5f0, 0.5f0)
    loc2 = Location(0.0f0, 0.0f0)
    loc3 = Location(1.0f0, 1.0f0)

    @test Mating.distance(loc1, loc1) == 0
    @test Mating.distance(loc1, loc2) ≈ sqrt(0.5)
    @test Mating.distance(loc1, loc3) ≈ sqrt(0.5)
    @test Mating.distance(loc2, loc3) ≈ sqrt(2)
end