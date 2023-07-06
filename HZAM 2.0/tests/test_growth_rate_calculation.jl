using Test
include("../src/population.jl")

using .Population
using QuadGK


@testset "get_ideal_densities" begin
    x_locations_F = [0.0f0, 0.3f0, 0.7f0, 0.8f0, 1.0f0]

    y_locations_F = fill(0.5f0, 5)

    locations_F = Location.(x_locations_F, y_locations_F)

    K_total = 40233

    sigma_comp = 0.01

    ideal_densities = Population.get_ideal_densities(K_total, sigma_comp, locations_F) .- Ref(1)

    max_density = 1000 * quadgk(x -> exp(-(x - 0.5)^2 / (2 * sigma_comp^2)), 0.47, 0.53)[1]

    @test abs(ideal_densities[1] - max_density / 2) < 0.001

    @test abs(ideal_densities[2] - max_density) < 0.001

    @test abs(ideal_densities[3] - max_density) < 0.001

    @test abs(ideal_densities[4] - max_density) < 0.001

    @test abs(ideal_densities[5] - max_density / 2) < 0.001
end

@testset "calculate_ind_useResource_no_ecolDiff" begin
    ecolDiff = 0.0
    competition_traits_F = [0, 1]

    ind_useResourceA_F = Population.calculate_ind_useResourceA(competition_traits_F, ecolDiff)
    ind_useResourceB_F = Population.calculate_ind_useResourceB(competition_traits_F, ecolDiff)

    @test ind_useResourceA_F == [0.5, 0.5]
    @test ind_useResourceB_F == [0.5, 0.5]
end

@testset "calculate_ind_useResource_with_ecolDiff" begin
    ecolDiff = 1.0
    competition_traits_F = [0, 1]

    ind_useResourceA_F = Population.calculate_ind_useResourceA(competition_traits_F, ecolDiff)
    ind_useResourceB_F = Population.calculate_ind_useResourceB(competition_traits_F, ecolDiff)


    @test ind_useResourceA_F == [1, 0]
    @test ind_useResourceB_F == [0, 1]
end


@testset "calc_growth_rates" begin
    K_total = Integer(trunc(1000 / sqrt(2 * pi * 0.01^2)))
    ecolDiff = 0.0
    competition_trait_loci = 1:1
    intrinsic_R = 1.1
    sigma_comp = 0.01

    deme = Population.Deme([fill(1, 2, 3), fill(1, 2, 3)], [], [Location(0.5f0, 0.5f0), Location(0.0f0, 0.0f0)], [], [1], [], [1], 0)
    empty_deme = Population.Deme([], [], [], [], [1], [], [1:1], 0)

    growth_rates_F = Population.calculate_growth_rates([empty_deme empty_deme empty_deme; empty_deme deme empty_deme; empty_deme empty_deme empty_deme],
        CartesianIndex(2, 2),
        K_total,
        sigma_comp,
        intrinsic_R)

    @test abs(growth_rates_F[1] - 1.1) < 0.005

    @test growth_rates_F[2] < 1.1
end
