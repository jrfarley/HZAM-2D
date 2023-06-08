using Test
include("../src/population.jl")

using .Population

#import .Population : get_genotypes, get_growth_rates, get_locations, get_ideal_densities, calculate_growth_rates_spatial, set_locations, generate_genotype_array, calc_traits_additive


@testset "initialize_population_sympatry" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = true
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp)

    @test get_population() == (2, 2)

    @test (Population.genotypes_F, Population.genotypes_M) == (genotypes, genotypes)

    @test (Population.growth_rate_resourceA, Population.growth_rate_resourceB) == (1.0, 1.0)
end


@testset "initialize_population_sympatry_dif_pop_ratio" begin
    K_total = 8
    starting_pop_ratio = 0.5
    ecolDiff = 0.0
    sympatry = true
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 0.48],
        starting_range_pop1=[5.2, 1.0])

    @test get_population() == (3, 3)

    @test (Population.genotypes_F, Population.genotypes_M) == (genotypes, genotypes)

    @test (Population.growth_rate_resourceA, Population.growth_rate_resourceB) == (1.0, 1.0)
end


@testset "initialize_population_sympatry_ecolDiff_0.5" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0.5
    sympatry = true
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp)

    @test get_population() == (2, 2)

    @test (Population.genotypes_F, Population.genotypes_M) == (genotypes, genotypes)

    @test (Population.growth_rate_resourceA, Population.growth_rate_resourceB) == (1.0, 1.0)
end

@testset "initialize_population_allopatry" begin
    K_total = 8
    starting_pop_ratio = 1.0
    ecolDiff = 1.0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 0.48],
        starting_range_pop1=[0.52, 1.0])

    @test get_population() == (2, 2)

    @test (Population.genotypes_F, Population.genotypes_M) == (genotypes, genotypes)
end

@testset "calc_growth_rates" begin
    K_total = 100
    starting_pop_ratio = 1.0
    ecolDiff = 0.0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 0.48],
        starting_range_pop1=[0.52, 1.0])



    spaced_locations_M = collect(Float32, 0.02:0.02:0.96)
    spaced_locations_F = collect(Float32, 0.01:0.02:0.95)

    Population.set_locations(spaced_locations_F, spaced_locations_M)

    growth_rates_A, growth_rates_B = Population.calculate_growth_rates_spatial()


    @test growth_rates_A[1] > 1.01 && growth_rates_A[1] < 1.03

    @test growth_rates_A[6] > 0.999 && growth_rates_A[6] < 1.0

end

@testset "generate_genotype_array" begin
    genotypes = Population.generate_genotype_array(1, 1, 3)

    @test genotypes[1] == [0 0 0; 0 0 0]
    @test genotypes[2] == [1 1 1; 1 1 1]
end

@testset "calc_traits_additive" begin
    genotypes = Population.generate_genotype_array(1, 1, 3)

    @test Population.calc_traits_additive(genotypes, 1:3) == [0, 1]
end

@testset "disperse_individual" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 1.0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    geographic_limits = [0.0, 1.0]
    sigma_disp = 0.01

    Population.set_locations([0.1, 0.7], [0.2, 0.4])

    num_fail = 0

    for i in 1:1000
        new_location = Population.disperse_individual(1, sigma_disp, geographic_limits)

        if new_location < 0.1 - 2 * sigma_disp || new_location > 0.1 + 2 * sigma_disp
            num_fail += 1
        end
    end

    @test num_fail < 70
end


@testset "calc_hybrid_indices" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0.0
    sympatry = true
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 0.48],
        starting_range_pop1=[5.2, 1.0])

    println(Population.calc_hybrid_indices([Population.genotypes_F; Population.genotypes_M]))

    @test Population.calc_hybrid_indices([Population.genotypes_F; Population.genotypes_M]) == [0.0, 1.0, 0.0, 1.0]

end

@testset "generate_offspring_genotype" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = true
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp)

    @test (Population.genotypes_F, Population.genotypes_M) == (genotypes, genotypes)

    @test Population.generate_offspring_genotype(1, 2, 3) == ([0 0 0; 1 1 1], 0.25)
end

@testset "assign_zones" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    locations_F, locations_M = ([0.1, 0.3, 0.5, 0.7], [0.2, 0.4, 0.6, 0.8])

    females_per_zone, males_per_zone = Population.assign_zones(locations_F, locations_M)

    @test females_per_zone[1] == []
    @test females_per_zone[4] == [2]
    @test males_per_zone[5] == [2]
end

@testset "choose_closest_male" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    locations_F, locations_M = ([0.1, 0.3, 0.5, 0.7], [0.15, 0.4, 0.6, 0.8])

    Population.set_locations(locations_F, locations_M)

    closest_male, elig_M = choose_closest_male(collect(1:4), 1)

    @test closest_male == 1
    @test elig_M == [2, 3, 4]

    closest_male, elig_M = choose_closest_male(elig_M, 1)
    @test closest_male == 2
    @test elig_M == []
end

@testset "calc_match_strength" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01
    pref_SD = sqrt(-1 / (2 * log(1 / (1 + 10^(-15)))))

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    locations_F, locations_M = ([0.1, 0.3, 0.5, 0.7], [0.15, 0.4, 0.6, 0.8])

    Population.set_locations(locations_F, locations_M)

    match_strength1 = calc_match_strength(1, 1, pref_SD)

    match_strength2 = calc_match_strength(2, 2, pref_SD)

    match_strength3 = calc_match_strength(1, 2, pref_SD)

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
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    competition_traits_F, competition_traits_M = ([0, 1], [0, 1])

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 0.50],
        starting_range_pop1=[0.50, 1.0])

    ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = Population.calculate_ind_useResource(competition_traits_F, competition_traits_M)

    @test ind_useResourceA_F == [0.5, 0.5]
    @test ind_useResourceB_F == [0.5, 0.5]
    @test ind_useResourceA_M == [0.5, 0.5]
    @test ind_useResourceB_M == [0.5, 0.5]
end

@testset "calculate_ind_useResource_with_ecolDiff" begin
    K_total = 4
    starting_pop_ratio = 1.0
    ecolDiff = 1.0
    sympatry = false
    total_loci = 3
    competition_trait_loci = 1:1
    female_mating_trait_loci = 2:2
    male_mating_trait_loci = 2:2
    hybrid_survival_loci = 3:3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    competition_traits_F, competition_traits_M = ([0, 1], [0, 1])

    initialize_population(K_total,
        starting_pop_ratio,
        ecolDiff,
        sympatry,
        total_loci,
        competition_trait_loci,
        female_mating_trait_loci,
        male_mating_trait_loci,
        hybrid_survival_loci,
        intrinsic_R,
        sigma_comp;
        new_geographic_limits=[0.0, 1.0],
        starting_range_pop0=[0, 1.0],
        starting_range_pop1=[0.0, 1.0])

    ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = Population.calculate_ind_useResource(competition_traits_F, competition_traits_M)

    @test ind_useResourceA_F == [1, 0]
    @test ind_useResourceB_F == [0, 1]
    @test ind_useResourceA_M == [1, 0]
    @test ind_useResourceB_M == [0, 1]
end

@testset "find_inactive_zones1" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0.25, 100)
    mitochondria_pop1 = fill(0.75, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == [1, 2, 3, 4, 7, 8, 9, 10]
end

@testset "find_inactive_zones2" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 99)
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0.25, 99)
    mitochondria_pop1 = fill(0.75, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = collect(Float32, 0.0:0.005:0.99)

    locations_M = collect(Float32, 0.001:0.005:0.991)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == [1, 2, 3, 7, 8, 9, 10]
end

@testset "find_inactive_zones3" begin
    genotypes_pop0 = vcat(fill([0 0 0; 0 0 0], 45), [fill(1, 2, 3)], fill([0 0 0; 0 0 0], 54))
    genotypes_pop1 = fill([1 1 1; 1 1 1], 100)
    mitochondria_pop0 = fill(0.25, 100)
    mitochondria_pop1 = fill(0.25, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == [1, 7, 8, 9, 10]
end

@testset "find_inactive_zones4" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = genotypes_pop0
    mitochondria_pop0 = vcat(fill(0.25, 50), 0.75, fill(0.25, 49))
    mitochondria_pop1 = fill(0.25, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == [1, 5, 6, 7, 8, 9, 10]
end

@testset "find_inactive_zones5" begin
    genotypes_pop0 = fill([0 0 0; 0 0 0], 100)
    genotypes_pop1 = vcat(fill([0 0 0; 0 0 0], 45), [fill(1, 2, 3)], fill([0 0 0; 0 0 0], 54))
    mitochondria_pop0 = vcat(fill(0.25, 50), 0.75, fill(0.25, 49))
    mitochondria_pop1 = fill(0.25, 100)

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = collect(Float32, 0.0:0.005:0.995)

    locations_M = collect(Float32, 0.001:0.005:0.996)

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == [1, 5, 6, 10]
end

@testset "find_inactive_zones6" begin
    genotypes_pop0 = [[0 0 0; 0 0 0]]
    genotypes_pop1 = [[0 0 0; 0 0 0]]
    mitochondria_pop0 = [0.25]
    mitochondria_pop1 = [0.25]

    mitochondria = [mitochondria_pop0; mitochondria_pop1]

    Population.set_mitochondria(mitochondria, mitochondria)

    println(length(genotypes_pop0))

    genotypes = [genotypes_pop0; genotypes_pop1]

    Population.set_genotypes(genotypes, genotypes)

    locations_F = [0.01, 0.99]

    locations_M = [0.02, 0.98]

    indv_per_zone_F, indv_per_zone_M = Population.assign_zones(locations_F, locations_M)

    dead_zones = Population.find_inactive_zones(indv_per_zone_F, indv_per_zone_M)

    @test dead_zones == []
end