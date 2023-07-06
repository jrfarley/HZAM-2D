using Test
using BenchmarkTools
using Random
include("../src/population.jl")

using .Population

@testset "choose_closest_male_no_demes" begin
    elig_M = collect(1:1001)
    location_mother = Location(0.5f0, 0.5f0)

    locations_M = Location.(Float32.(rand(1000)), Float32.(rand(1000)))
    male_loc = Location(0.5f0, 0.5f0)
    locations_M = shuffle(push!(locations_M, male_loc))

    closest_male = choose_closest_male(elig_M, locations_M, location_mother)

    @test locations_M[closest_male] == male_loc

    for i in 1:100
        male_loc = Location(0.55f0, 0.45f0)
        locations_M = Location.(Float32.(rand(1000) .* Ref(0.4)), Float32.(rand(1000) .* Ref(0.4) .+ Ref(0.55)))
        locations_M = shuffle(push!(locations_M, male_loc))

        closest_male = choose_closest_male(elig_M, locations_M, location_mother)

        @test locations_M[closest_male] == male_loc
    end
end

@testset "choose_closest_male_one_deme" begin
    intervals = collect(0.0f0:0.1f0:0.99f0) # locations of the demes along an axis

    deme_locations = Location.(intervals, intervals') # location of each deme (lower left corner)

    deme_populations = map(l -> l.x < 0.5 ? 1 : 1, deme_locations) # assigns 0 to each deme occupied by population A and 1 to each deme occupied by population B (dividing line down the middle at the beginning)

    # initializes the genotypes, locations, and mitochondria for each deme
    demes = Population.Deme.(Ref(1000),
        Ref(10),
        deme_locations,
        Ref(0.1f0),
        deme_populations,
        Ref(0))

    location_mother = demes[5,5].locations_F[1]

    elig_M = Dict{CartesianIndex,Vector{Int64}}()

    for i in eachindex(IndexCartesian(), demes)
        elig_M[i] = 1:length(demes[i].locations_M)
    end

    closest_male, deme_index = choose_closest_male(demes, [CartesianIndex(5, 5)], elig_M, location_mother, 0.1f0)

    @test demes[deme_index].locations_M[closest_male] == location_mother
    @test deme_index == CartesianIndex(5, 5)
end

@testset "choose_closest_male_multiple_demes" begin
    intervals = collect(0.0f0:0.1f0:0.99f0) # locations of the demes along an axis

    deme_locations = Location.(intervals, intervals') # location of each deme (lower left corner)

    deme_populations = map(l -> l.x < 0.5 ? 1 : 1, deme_locations) # assigns 0 to each deme occupied by population A and 1 to each deme occupied by population B (dividing line down the middle at the beginning)

    for i in 1:10
        # initializes the genotypes, locations, and mitochondria for each deme
        demes = Population.Deme.(Ref(10),
            Ref(10),
            deme_locations,
            Ref(0.1f0),
            deme_populations,
            Ref(0))

        location_mother = Location(0.45f0, 0.45f0)

        locations_M = vcat(demes[i, i].locations_M, location_mother)

        genotypes_F = demes[i, i].genotypes_F
        genotypes_M = demes[i, i].genotypes_M
        push!(genotypes_M, genotypes_M[1])
        locations_F = demes[i, i].locations_F
        mitochondria_F = demes[i, i].mitochondria_F
        mitochondria_M = vcat(demes[i, i].mitochondria_M, 1)

        demes[i, i] = Population.Deme(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M, collect(1:2), 0)

        elig_M = Dict{CartesianIndex,Vector{Int64}}()

        for i in eachindex(IndexCartesian(), demes)
            elig_M[i] = 1:length(demes[i].locations_M)
        end
        neighbourhood = filter(e -> length(elig_M[e]) > 0, CartesianIndex(1, 1):CartesianIndex(10, 10))

        closest_male, deme_index = choose_closest_male(demes, neighbourhood, elig_M, location_mother, 1.0f0)

        @test demes[deme_index].locations_M[closest_male] == location_mother
        @test deme_index == CartesianIndex(i, i)
    end
end