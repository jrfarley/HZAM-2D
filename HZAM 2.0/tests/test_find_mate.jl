using Test
using BenchmarkTools
include("../src/population.jl")

using .Population

locations = Matrix{Vector{Location}}(undef, 10, 10)

for i in eachindex(IndexCartesian(), locations)
    locations[i] = Location.(fill([Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)], 10000))
end

function foo(locations)
    locations_M = vcat(locations[CartesianIndex(4, 4):CartesianIndex(6, 6)]...)
    elig_M = collect(1:length(locations_M))
    for i in 1:10
        mate, elig_M = choose_closest_male(elig_M,
            locations_M,
            Location(0.5f0, 0.5f0))
    end
end
a = [[] for i in eachindex(IndexCartesian(), locations)]
demes = Population.Deme.(a, a, a, locations, a, a, Ref(0:0), Ref(0))

function bar(locations, demes)
    count = 1
    neighbourhood_size = 0
    unused_deme_indices = collect(CartesianIndex(1, 1):CartesianIndex(10, 10))
    elig_M = Dict{CartesianIndex,Vector{Int64}}()
    while count <= 2
        lower_left = max(CartesianIndex(6, 6) - CartesianIndex(neighbourhood_size, neighbourhood_size), CartesianIndex(1, 1))
        upper_right = min(CartesianIndex(6, 6) + CartesianIndex(neighbourhood_size, neighbourhood_size), CartesianIndex(10, 10))
        neighbourhood = filter(e -> e ∈ unused_deme_indices, lower_left:upper_right) # remove used demes

        for i in neighbourhood
            elig_M[i] = 1:length(locations[i])
        end

        if sum(map(length, values(elig_M))) == 0
            continue
        end

        unused_deme_indices = filter(e -> e ∉ neighbourhood, unused_deme_indices) # remove used demes from list
        for i in 1:5
            focal_male, elig_M, male_deme = choose_closest_male(demes, neighbourhood, elig_M, Location(0.5f0, 0.5f0))
        end
        neighbourhood_size += 1
        count += 1
    end
end


@btime foo(locations)
@btime bar(locations, demes)