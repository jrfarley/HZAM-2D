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
    geographic_limits = [Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)]
    starting_range_pop0 = [Location(0.0f0, 0.0f0), Location(0.5f0, 1.0f0)]
    starting_range_pop1 = [Location(0.5f0, 0.0f0), Location(1.0f0, 1.0f0)]

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
    geographic_limits = [Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)]
    starting_range_pop0 = [Location(0.0f0, 0.0f0), Location(0.75f0, 1.0f0)]
    starting_range_pop1 = [Location(0.75f0, 0.0f0), Location(1.0f0, 1.0f0)]

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
    geographic_limits = [Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)]
    starting_range_pop0 = [Location(0.0f0, 0.0f0), Location(0.5f0, 1.0f0)]
    starting_range_pop1 = [Location(0.5f0, 0.0f0), Location(1.0f0, 1.0f0)]

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

    x_locations_F = [l.x for l in pd.locations_F[pd.active_F]]
    x_locations_M = [l.x for l in pd.locations_M[pd.active_M]]
    x_locations_M_buffer1 = [l.x for l in pd.locations_M[pd.buffer_M[1]]]
    x_locations_M_buffer2 = [l.x for l in pd.locations_M[pd.buffer_M[2]]]

    @test length(pd.active_F) == 0 || (maximum(x_locations_F) <= 0.6 && minimum(x_locations_F) >= 0.4)
    @test length(pd.active_M) == 0 || (maximum(x_locations_M) <= 0.7 && minimum(x_locations_M) >= 0.3)

    @test length(pd.buffer_M[1]) == 0 || (maximum(x_locations_M_buffer1) <= 0.4 && minimum(x_locations_M_buffer1) >= 0.3)
    @test length(pd.buffer_M[2]) == 0 || (maximum(x_locations_M_buffer2) <= 0.7 && minimum(x_locations_M_buffer2) >= 0.6)

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
    locations = Location[]
    geographic_limits = [Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)]
    sigma_disp = 0.01

    num_fail = 0

    for i in 1:1000
        new_location = Location(Location(0.5f0, 0.5f0), sigma_disp, geographic_limits)
        push!(locations, new_location)
        if new_location.x < 0.5 - 2 * sigma_disp || new_location.x > 0.5 + 2 * sigma_disp || new_location.y < 0.5 - 2 * sigma_disp || new_location.y > 0.5 + 2 * sigma_disp
            num_fail += 1
        end
    end

    average_x = sum([l.x for l in locations]) / 1000

    average_y = sum([l.y for l in locations]) / 1000

    @test num_fail < 70

    @test average_x > 0.495 && average_x < 0.505
    @test average_y > 0.495 && average_y < 0.505

end

@testset "generate_offspring_genotype" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    @test Population.generate_offspring_genotype(genotypes[1], genotypes[2]) == [0 0 0; 1 1 1]
    @test Population.generate_offspring_genotype(genotypes[1], genotypes[1]) == [0 0 0; 0 0 0]
end

@testset "assign_zones" begin
    x_locations_F, x_locations_M = ([0.1f0, 0.3f0, 0.5f0, 0.7f0], [0.2f0, 0.4f0, 0.75f0, 0.72f0])

    y_locations_F, y_locations_M = ([0.1f0, 0.4f0, 0.2f0, 0.6f0], [0.8f0, 0.3f0, 0.64f0, 0.69f0])

    locations_F = Location.(x_locations_F, y_locations_F)
    locations_M = Location.(x_locations_M, y_locations_M)

    females_per_zone, males_per_zone = Population.assign_zones(locations_F, locations_M)

    @test females_per_zone[1,1] == []
    @test females_per_zone[4,5] == [2]
    @test males_per_zone[5,4] == [2]

    @test males_per_zone[8,7] == [3,4]
end

@testset "choose_closest_male" begin
    x_locations_M = [0.15f0, 0.15f0, 0.05f0, 0.8f0]

    y_locations_M = [0.1f0, 0.15f0, 0.2f0, 0.9f0]
    
    locations_M = Location.(x_locations_M, y_locations_M)

    location_F = Location(0.1f0, 0.1f0)

    closest_male, elig_M = choose_closest_male(collect(1:4), locations_M, location_F)

    @test closest_male == 1
    @test elig_M == [2, 3, 4]

    closest_male, elig_M = choose_closest_male(elig_M, locations_M, location_F)
    @test closest_male == 2
    @test elig_M == [3,4]

    
    closest_male, elig_M = choose_closest_male(elig_M, locations_M, location_F)
    @test closest_male == 3
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



#=
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
end=#

@testset "mean" begin
    genotypes = [[1 0 1; 0 1 0], [0 0 1; 1 0 0], [1 1 1; 1 1 1]]

    @test Population.mean(genotypes[1], collect(1:2)) == 0.5

    @test Population.mean(genotypes[2], [2]) == 0
    @test Population.mean(genotypes[3], collect(1:3)) == 1
end