using QuadGK

@testset "compare_ideal_densities_to_1D" begin
    x_locations_F = [0.0f0, 0.3f0, 0.7f0, 0.8f0, 1.0f0]

    y_locations_F = fill(0.5f0, 5)

    locations_F = Location.(x_locations_F, y_locations_F)

    K_total = 40233

    sigma_comp = 0.01

    ideal_densities =
        Population.calc_ideal_densities(
            K_total,
            sigma_comp,
            locations_F,
            0.03
        )

    ideal_densities .-= 1

    max_density = 1000 * quadgk(x -> exp(-(x - 0.5)^2 / (2 * sigma_comp^2)), 0.47, 0.53)[1]

    @test abs(ideal_densities[1] - max_density / 2) < 0.001

    @test abs(ideal_densities[2] - max_density) < 0.001

    @test abs(ideal_densities[3] - max_density) < 0.001

    @test abs(ideal_densities[4] - max_density) < 0.001

    @test abs(ideal_densities[5] - max_density / 2) < 0.001
end


@testset "ideal_density_approximation" begin
    x_locations_F = [0.0f0, 0.3f0, 0.7f0, 0.8f0, 1.0f0, 1.0f0]

    y_locations_F = [0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.0f0]

    locations_F = Location.(x_locations_F, y_locations_F)

    K_total = 10000

    sigma_comp = 0.01

    ideal_densities =
        Population.calc_ideal_densities(
            K_total,
            sigma_comp,
            locations_F,
            0.03
        )

    ideal_densities .-= 1

    max_density =
        K_total *
        quadgk(
            t -> quadgk(r -> r * exp(-(r^2) / (2 * sigma_comp^2)), 0, 0.03)[1],
            0,
            2 * pi
        )[1]

    @test ideal_densities[1] ≈ max_density / 2

    @test ideal_densities[2] ≈ max_density

    @test ideal_densities[3] ≈ max_density

    @test ideal_densities[4] ≈ max_density

    @test ideal_densities[5] ≈ max_density / 2

    @test ideal_densities[6] ≈ max_density / 4
end

@testset "calc_ind_useResource_no_ecolDiff" begin
    ecolDiff = 0.0
    competition_traits_F = [0, 1, 0.5, 0.25]

    ind_useResourceA_F = Population.calc_ind_useResourceA(competition_traits_F, ecolDiff)
    ind_useResourceB_F = Population.calc_ind_useResourceB(competition_traits_F, ecolDiff)

    @test ind_useResourceA_F == [0.5, 0.5, 0.5, 0.5]
    @test ind_useResourceB_F == [0.5, 0.5, 0.5, 0.5]
end

@testset "calc_ind_useResource_ecolDiff1" begin
    ecolDiff = 1.0
    competition_traits_F = [0, 1, 0.5, 0.25]

    ind_useResourceA_F = Population.calc_ind_useResourceA(competition_traits_F, ecolDiff)
    ind_useResourceB_F = Population.calc_ind_useResourceB(competition_traits_F, ecolDiff)


    @test ind_useResourceA_F == [1, 0, 0.5, 0.75]
    @test ind_useResourceB_F == [0, 1, 0.5, 0.25]
end


@testset "calc_growth_rates" begin
    K_total = Integer(trunc(1000 / sqrt(2 * pi * 0.01^2)))
    ecolDiff = 0.0
    competition_trait_loci = 1:1
    intrinsic_R = 1.1
    sigma_comp = 0.01

    zone = Population.Zone(
        [fill(Int8(1), 2, 3), fill(Int8(1), 2, 3)],
        Matrix{Int8}[],
        [Location(0.5f0, 0.5f0),
            Location(0.0f0, 0.0f0)],
        Location[],
        [1],
        0
    )
    empty_zone = Population.Zone(
        Matrix{Int8}[],
        Matrix{Int8}[],
        Location[],
        Location[],
        1:1,
        0
    )

    growth_rates_F = Population.calc_growth_rates(
        [empty_zone empty_zone empty_zone; 
        empty_zone zone empty_zone; 
        empty_zone empty_zone empty_zone],
        CartesianIndex(2, 2),
        K_total,
        sigma_comp,
        intrinsic_R)

    @test abs(growth_rates_F[1] - 1.1) < 0.005

    @test growth_rates_F[2] < 1.1
end

@testset "calc_real_densities_ecolDiff0" begin
    K_total = Integer(trunc(1000 / sqrt(2 * pi * 0.01^2)))
    ecolDiff = 0.0
    competition_trait_loci = 1:1
    intrinsic_R = 1.1
    sigma_comp = 0.01

    zone = Population.Zone(
        [fill(Int8(1), 2, 3), fill(Int8(1), 2, 3)],
        Matrix{Int8}[],
        [Location(0.5f0, 0.5f0), Location(0.0f0, 0.0f0)],
        Location[], [Int8(1)],
        0
    )
    empty_zone = Population.Zone(
        Matrix{Int8}[],
        Matrix{Int8}[],
        Location[],
        Location[],
        1:1,
        0
    )

    neighbourhood = [
        empty_zone empty_zone empty_zone
        empty_zone zone empty_zone
        empty_zone empty_zone empty_zone
    ]
    locations_F = [Location(1.0f0, 0.0f0), Location(0.5f0, 0.5f0), Location(0.49f0, 0.5f0)]

    densities_A, densities_B = Population.calc_real_densities(
        neighbourhood,
        locations_F,
        0.03,
        0.01
    )

    @test densities_A[1] == 0
    @test densities_B[1] == 0
    @test densities_A[2] == 0.5
    @test densities_B[2] == 0.5
    @test abs(densities_A[3] - (1 / 2) * exp(-1 / 2)) < 0.00001
    @test abs(densities_B[3] - (1 / 2) * exp(-1 / 2)) < 0.00001
end

@testset "calc_real_densities_ecolDiff1" begin
    K_total = Integer(trunc(1000 / sqrt(2 * pi * 0.01^2)))
    ecolDiff = 1
    competition_trait_loci = 1:1
    intrinsic_R = 1.1
    sigma_comp = 0.01

    zone = Population.Zone(
        [fill(Int8(1), 2, 3), fill(Int8(1), 2, 3)],
        Matrix{Int8}[],
        [Location(0.5f0, 0.5f0), Location(0.0f0, 0.0f0)],
        Location[],
        [1],
        1
    )
    empty_zone = Population.Zone(
        Matrix{Int8}[],
        Matrix{Int8}[],
        Location[],
        Location[],
        1:1,
        1
    )

    neighbourhood = [
        empty_zone empty_zone empty_zone
        empty_zone zone empty_zone
        empty_zone empty_zone empty_zone
    ]
    locations_F = [Location(1.0f0, 0.0f0), Location(0.5f0, 0.5f0), Location(0.49f0, 0.5f0)]

    densities_A, densities_B =
        Population.calc_real_densities(neighbourhood, locations_F, 0.03, 0.01)

    @test densities_A[1] == 0
    @test densities_B[1] == 0
    @test densities_A[2] == 0
    @test densities_B[2] == 1
    @test densities_A[3] == 0
    @test abs(densities_B[3] - exp(-1 / 2)) < 0.00001
end

@testset "max_radius_squared" begin
    locations = Location.(
        [0.98f0, 0.98f0, 0.01f0, 0.01f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.01f0, 0.02f0,
            0.98f0, 0.99f0],
        [0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0, 0.98f0, 0.98f0, 0.01f0, 0.01f0, 0.01f0, 0.021f0,
            0.99f0, 0.98f0]
    )

    angles = [0, 3 * pi / 2, 2 * pi / 3, pi / 4, 2, pi / 2, 0, 5 * pi / 4, pi + 0.01,
        3 * pi / 4, 3 * pi / 4, pi / 2, pi / 2]

    max_radii_squared = Population.max_radius_squared.(locations, angles, Ref(0.03))

    expected = [0.02^2, 0.03^2, 0.02^2, 0.03^2, 0.03^2, 0.02^2, 0.03^2, 0.0002, 0.03^2,
        0.0002, 0.0008, 0.0001, 0.0004]
    for i in eachindex(max_radii_squared)
        @test abs(max_radii_squared[i] - expected[i]) < 0.0000001
    end
end
