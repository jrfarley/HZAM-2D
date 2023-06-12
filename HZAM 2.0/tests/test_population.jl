using Test
include("../src/population.jl")

using .Population

@testset "initialize_population_1" begin
    K_total = 8
    ecolDiff = 0
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]
    starting_range_pop0 = [0.0, 0.5]
    starting_range_pop1 = [0.5, 1.0]

    pd = PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp,
        geographic_limits,
        starting_range_pop0,
        starting_range_pop1,
        false)

    @test pd.genotypes_F == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 1 1 1]]
    @test pd.genotypes_M == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    @test pd.mitochondria_F == [0, 0, 1, 1]
    @test pd.mitochondria_M == [0, 0, 1, 1]

    @test length(pd.locations_F) == 4
    @test length(pd.locations_M) == 4

    @test length(pd.growth_rates_F) == 4
end

@testset "initialize_population_2" begin
    K_total = 8
    ecolDiff = 0
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]
    starting_range_pop0 = [0.0, 0.75]
    starting_range_pop1 = [0.75, 1.0]

    pd = PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp,
        geographic_limits,
        starting_range_pop0,
        starting_range_pop1,
        false)

    @test pd.genotypes_F == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1]]
    @test pd.genotypes_M == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    @test pd.mitochondria_F == [0, 0, 0, 1]
    @test pd.mitochondria_M == [0, 0, 0, 1]

    @test length(pd.locations_F) == 4
    @test length(pd.locations_M) == 4

    @test length(pd.growth_rates_F) == 4
end

@testset "initialize_population_3" begin
    K_total = 8
    ecolDiff = 0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]
    starting_range_pop0 = [0.0, 0.5]
    starting_range_pop1 = [0.5, 1.0]

    pd = PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp,
        geographic_limits,
        starting_range_pop0,
        starting_range_pop1,
        true)

    @test pd.genotypes_F == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 1 1 1]]
    @test pd.genotypes_M == [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    @test pd.mitochondria_F == [0, 0, 1, 1]
    @test pd.mitochondria_M == [0, 0, 1, 1]

    @test length(pd.locations_F) == 4
    @test length(pd.locations_M) == 4

    @test length(pd.growth_rates_F) == length(pd.active_F)

    @test pd.left_boundary == 5
    @test pd.right_boundary == 6

    @test pd.inactive_locations_F == pd.locations_F
    @test pd.inactive_genotypes_F == pd.genotypes_F
    @test pd.inactive_mitochondria_F == pd.mitochondria_F

    @test length(pd.active_F) == 0 || (maximum(pd.locations_F[pd.active_F]) <= 0.6 && minimum(pd.locations_F[pd.active_F]) >= 0.4)
    @test length(pd.active_M) == 0 || (maximum(pd.locations_M[pd.active_M]) <= 0.7 && minimum(pd.locations_M[pd.active_M]) >= 0.3)

    @test length(pd.buffer_M[1]) == 0 || (maximum(pd.locations_M[pd.buffer_M[1]]) <= 0.4 && minimum(pd.locations_M[pd.buffer_M[1]]) >= 0.3)
    @test length(pd.buffer_M[2]) == 0 || (maximum(pd.locations_M[pd.buffer_M[2]]) <= 0.7 && minimum(pd.locations_M[pd.buffer_M[2]]) >= 0.6)

end

@testset "calc_growth_rates" begin
    K_total = 100
    ecolDiff = 0.0
    competition_trait_loci = 1:1
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = Population.generate_genotype_array(8, 8, 3)


    locations_M = collect(Float32, 0.36:0.02:0.66)
    locations_F = collect(Float32, 0.35:0.02:0.65)

    competition_traits_F = calc_traits_additive(genotypes, competition_trait_loci)
    competition_traits_M = calc_traits_additive(genotypes, competition_trait_loci)

    ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = Population.calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)

    growth_rates_F = Population.calculate_growth_rates(ind_useResourceA_F,
        ind_useResourceB_F,
        ind_useResourceA_M,
        ind_useResourceB_M,
        locations_F,
        locations_M,
        K_total,
        sigma_comp,
        intrinsic_R,
        collect(1:16))

    @test growth_rates_F[1] > growth_rates_F[8]

    @test growth_rates_F[16] > growth_rates_F[8]

    @test growth_rates_F[1] > 1 && growth_rates_F[1] < 1.1

    @test growth_rates_F[8] > 0.995 && growth_rates_F[8] < 1.005
end

@testset "generate_genotype_array" begin
    genotypes = Population.generate_genotype_array(1, 1, 3)

    @test genotypes[1] == [0 0 0; 0 0 0]
    @test genotypes[2] == [1 1 1; 1 1 1]

    genotypes = Population.generate_genotype_array(1, 2, 2)

    @test genotypes[1] == [0 0; 0 0]
    @test genotypes[2] == [1 1; 1 1]
    @test genotypes[3] == [1 1; 1 1]
end

@testset "calc_traits_additive" begin
    genotypes = Population.generate_genotype_array(1, 1, 3)

    @test calc_traits_additive(genotypes, 1:3) == [0, 1]

    genotypes = [[1 0 0; 0 0 0], [0 1 0; 1 0 0], [1 1 1; 1 1 1]]

    @test calc_traits_additive(genotypes, 1:2) == [0.25, 0.5, 1]
end

@testset "disperse_individual" begin
    geographic_limits = [0.0, 1.0]
    sigma_disp = 0.01

    num_fail = 0

    for i in 1:1000
        new_location = Population.disperse_individual(0.1, sigma_disp, geographic_limits)

        if new_location < 0.1 - 2 * sigma_disp || new_location > 0.1 + 2 * sigma_disp
            num_fail += 1
        end
    end

    @test num_fail < 70
end

@testset "generate_offspring_genotype" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    @test Population.generate_offspring_genotype(genotypes[1], genotypes[2]) == [0 0 0; 1 1 1]
    @test Population.generate_offspring_genotype(genotypes[1], genotypes[1]) == [0 0 0; 0 0 0]
end

@testset "assign_zones" begin
    locations_F, locations_M = ([0.1, 0.3, 0.5, 0.7], [0.2, 0.4, 0.6, 0.8])

    females_per_zone, males_per_zone = Population.assign_zones(locations_F, locations_M)

    @test females_per_zone[1] == []
    @test females_per_zone[4] == [2]
    @test males_per_zone[5] == [2]
end

@testset "choose_closest_male" begin
    locations_M = [0.15, 0.4, 0.6, 0.8]

    closest_male, elig_M = choose_closest_male(collect(1:4), locations_M, 0.1)

    @test closest_male == 1
    @test elig_M == [2, 3, 4]

    closest_male, elig_M = choose_closest_male(elig_M, locations_M, 0.1)
    @test closest_male == 2
    @test elig_M == []
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

@testset "get_ideal_densities" begin
    locations_F = [0.0, 0.3, 0.7, 0.8, 1.0]

    K_total = 1000
    sigma_comp = 0.01

    ideal_densities = Population.get_ideal_densities(K_total, sigma_comp, locations_F)

    @test ideal_densities[1] > 12.5 && ideal_densities[1] < 12.55

    @test ideal_densities[2] > 25.066 && ideal_densities[2] < 25.067

    @test ideal_densities[3] > 25.066 && ideal_densities[3] < 25.067

    @test ideal_densities[4] > 25.066 && ideal_densities[4] < 25.067

    @test ideal_densities[5] > 12.5 && ideal_densities[5] < 12.55
end

@testset "calculate_ind_useResource_no_ecolDiff" begin
    ecolDiff = 0.0
    competition_traits_F, competition_traits_M = ([0, 1], [0, 1])

    ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = Population.calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)

    @test ind_useResourceA_F == [0.5, 0.5]
    @test ind_useResourceB_F == [0.5, 0.5]
    @test ind_useResourceA_M == [0.5, 0.5]
    @test ind_useResourceB_M == [0.5, 0.5]
end

@testset "calculate_ind_useResource_with_ecolDiff" begin
    ecolDiff = 1.0
    competition_traits_F, competition_traits_M = ([0, 1], [0, 1])

    ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = Population.calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)

    @test ind_useResourceA_F == [1, 0]
    @test ind_useResourceB_F == [0, 1]
    @test ind_useResourceA_M == [1, 0]
    @test ind_useResourceB_M == [0, 1]
end

@testset "find_inactive_zones1" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0, 100)
    mitochondria_pop1 = fill(1, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]


    genotypes = [genotypes_pop0; genotypes_pop1]

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == [1, 2, 3, 4, 7, 8, 9, 10]
end

@testset "find_inactive_zones2" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 99)
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0, 99)
    mitochondria_pop1 = fill(1, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    genotypes = [genotypes_pop0; genotypes_pop1]


    locations_F = collect(Float32, 0.0:0.005:0.99)

    locations_M = collect(Float32, 0.001:0.005:0.991)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == [1, 2, 3, 7, 8, 9, 10]
end

@testset "find_inactive_zones3" begin
    genotypes_pop0 = vcat(fill([0 0 0; 0 0 0], 45), [fill(1, 2, 3)], fill([0 0 0; 0 0 0], 54))
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0, 100)
    mitochondria_pop1 = fill(1, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == [1, 7, 8, 9, 10]
end

@testset "find_inactive_zones4" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = vcat(fill(0, 50), 1, fill(0, 49))
    mitochondria_pop1 = fill(1, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == [1, 7, 8, 9, 10]
end

@testset "find_inactive_zones5" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = vcat(fill([1 1 1; 1 1 1], 45), [fill(0, 2, 3)], fill([1 1 1; 1 1 1], 54))
    mitochondria_pop0 = vcat(fill(0, 50), 1, fill(0, 49))
    mitochondria_pop1 = fill(1, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == [1, 10]
end

@testset "find_inactive_zones6" begin
    genotypes_pop0 = [[0 0 0; 0 0 0]]
    genotypes_pop1 = [[0 0 0; 0 0 0]]
    mitochondria_pop0 = [0]
    mitochondria_pop1 = [0]

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    locations_F = [0.01, 0.99]

    locations_M = [0.02, 0.98]

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F,
        indv_per_zone_M,
        genotypes,
        genotypes,
        mitochondria,
        mitochondria)

    @test dead_zones == []
end

@testset "mean" begin
    genotypes = [[1 0 1; 0 1 0], [0 0 1; 1 0 0], [1 1 1; 1 1 1]]

    @test Population.mean(genotypes[1], collect(1:2)) == 0.5

    @test Population.mean(genotypes[2], [2]) == 0
    @test Population.mean(genotypes[3], collect(1:3)) == 1
end