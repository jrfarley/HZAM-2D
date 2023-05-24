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

    genotypes = [0 0 0; 0 0 0;;; 1 1 1; 1 1 1]

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

    @test Population.get_genotypes() == (genotypes, genotypes)

    @test Population.get_growth_rates() == (1.0, 1.0)
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

    genotypes = [0 0 0; 0 0 0;;; 0 0 0; 0 0 0;;; 1 1 1; 1 1 1]

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

    @test Population.get_genotypes() == (genotypes, genotypes)

    @test Population.get_growth_rates() == (1.0, 1.0)
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

    genotypes = [0 0 0; 0 0 0;;; 1 1 1; 1 1 1]

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

    @test Population.get_genotypes() == (genotypes, genotypes)

    @test Population.get_growth_rates() == (1.0, 1.0)
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

    genotypes = [0 0 0; 0 0 0;;; 1 1 1; 1 1 1]

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

    @test Population.get_genotypes() == (genotypes, genotypes)
end

@testset "get_ideal_densities" begin
    geographic_limits = [0.0, 1.0]
    spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])
    ideal_densities = Population.get_ideal_densities(1000, 0.01, spaced_locations)

    @test ideal_densities[1] > 12.5 && ideal_densities[1] < 12.55

    @test ideal_densities[500] > 25.066 && ideal_densities[500] < 25.067
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

    genotypes = [0 0 0; 0 0 0;;; 1 1 1; 1 1 1]

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


    @test growth_rates_A[1]>1.01 && growth_rates_A[1]<1.03

    @test growth_rates_A[6]>0.999 && growth_rates_A[6]<1.0

end

@testset "generate_genotype_array" begin
    genotypes = Population.generate_genotype_array(1,1,3)

    @test genotypes[:,:,1] == [0 0 0; 0 0 0]
    @test genotypes[:,:,2] == [1 1 1; 1 1 1]
end

@testset "calc_traits_additive" begin
    genotypes = Population.generate_genotype_array(1,1,3)

    @test Population.calc_traits_additive(genotypes) == [0,1]
end