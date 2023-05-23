using Test
include("../src/population.jl")

using .Population

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

    @test get_genotypes() == (genotypes, genotypes)

    @test get_growth_rates() == (1.0, 1.0)
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

    @test get_genotypes() == (genotypes, genotypes)

    @test get_growth_rates() == (1.0, 1.0)
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

    @test get_genotypes() == (genotypes, genotypes)

    @test get_growth_rates() == (1.0, 1.0)
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

    @test get_genotypes() == (genotypes, genotypes)

    #@test get_growth_rates() == (1.0, 1.0)

    #println(get_locations())

    spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])

    ind_locations_if_even_at_K = range(geographic_limits[1], geographic_limits[2], length=K_total)
    function get_density_if_even_at_K(focal_location) # this function calculates local density according to a normal curve
        sum(exp.(-((ind_locations_if_even_at_K .- focal_location).^2)./(2*(sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end
    ideal_densities_at_spaced_locations = map(get_density_if_even_at_K, spaced_locations)

    println(ideal_densities_at_spaced_locations)

    #@test get_ideal_densities_taylor(K_total, sigma_comp, spaced_locations) == ideal_densities_at_spaced_locations

    
end