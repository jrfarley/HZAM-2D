@testset "calc_all_trait_correlations" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 1], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)

    for tc in trait_correlations
        @test tc[3] == 1
    end

    genotypes = [[1 1 1; 1 1 1], [0 1 0; 0 1 0], [0 1 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)
    expected = [1, 0.9333, 0.93333, 0.9333, 1, 1, 0.9333, 1, 1]

    trait_correlations = [tc[3] for tc in trait_correlations]

    @test all(abs.(trait_correlations .- expected) .< 0.001)
end

@testset "calc_trait_correlation" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    @test abs(0.6708 - DataAnalysis.calc_trait_correlation(genotypes, [1, 3], [2])) < 0.0001
    @test abs(1 - DataAnalysis.calc_trait_correlation(genotypes, 1:3, 1:3)) < 0.0001
    @test abs(0.4472 - DataAnalysis.calc_trait_correlation(genotypes, [1], [3])) < 0.0001
    @test abs(0.8407 - DataAnalysis.calc_trait_correlation(genotypes, 1:2, 2:3)) < 0.0001

    genotypes = [[1 0 1; 0 1 1], [0 1 1; 1 0 1], [0 0 1; 1 1 1]]

    correlation = DataAnalysis.calc_trait_correlation(genotypes, [1, 2, 3], [3])

    @test correlation == 0
end

@testset "calc_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 1)

    @test linkage_diseq == 1

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 2)

    @test linkage_diseq ≈ 1 / 9

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 3)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 3, 1)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 2, 3)

    @test linkage_diseq ≈ 1 / 2
end

@testset "calc_average_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]
    loci = (overall=1:3, l1=1:2, l2=[1, 3])

    average_linkage_diseq = DataAnalysis.calc_average_linkage_diseq(genotypes, loci)

    expected = [2576 / 5184, 615 / 1296, 103 / 216, 205 / 432, 5 / 9, 125 / 288, 103 / 216, 125 / 288, 9 / 16]

    averages = vcat([a[3] for a in average_linkage_diseq]...)

    @test averages ≈ expected
end

@testset "calc_length" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [[fill(0, i); fill(1, 100 - i)] for i in 40:5:85]
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    length = DataAnalysis.calc_length(sigmoid_curves)

    @test abs(length - 1.106) < 0.001

    hybrid_indices = fill([fill(0, 50); fill(1, 50)], 10)
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    length = DataAnalysis.calc_length(sigmoid_curves)

    @test length ≈ 1
end

@testset "calc_position" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [[fill(0, i); fill(1, 100 - i)] for i in 40:5:85]
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    position = DataAnalysis.calc_position(sigmoid_curves)

    @test position ≈ 0.621
end

@testset "calc_variance" begin
    positions = collect(0:0.01:0.99)

    @test DataAnalysis.calc_variance(positions) ≈ 0.01

    positions = zeros(500)

    @test DataAnalysis.calc_variance(positions) == 0
end

@testset "count_phenotypes_at_loci" begin
    genotypes = [fill([0 1; 0 1], 100); fill([1 0; 0 0], 100)]

    c1 = DataAnalysis.count_phenotypes_at_loci(genotypes, [1])
    c2 = DataAnalysis.count_phenotypes_at_loci(genotypes, [1, 2])

    @test c1[0.0f0] == 100
    @test c1[0.5f0] == 100
    @test c1[1.0f0] == 0
    @test c2[0.0f0] == 0
    @test c2[0.25f0] == 100
    @test c2[0.5f0] == 100
    @test c2[0.75f0] == 0
    @test c2[1.0f0] == 0
end

@testset "find_fixed_alleles" begin
    genotypes = [[1 1 1; 1 1 1], [1 1 1; 1 1 1], [1 0 0; 1 0 0], [1 0 1; 1 1 0]]

    @test DataAnalysis.find_fixed_alleles(genotypes) == [1]

    genotypes = fill([1 1 1; 1 1 1], 50)

    @test DataAnalysis.find_fixed_alleles(genotypes) == [1, 2, 3]

    genotypes = [[1 1; 1 1], [0 0; 0 0]]

    @test DataAnalysis.find_fixed_alleles(genotypes) == []
end

@testset "calc_chi_squared" begin
    genotypes = [[1 1 1; 1 1 1], [1 0 1; 1 1 1], [1 1 0; 1 1 0], [0 0 0; 0 0 0],
        [0 1 0; 1 0 1]]

    x = DataAnalysis.calc_chi_squared(genotypes, 1)

    @test abs(x - 1.3719) < 0.001

    x = DataAnalysis.calc_chi_squared(genotypes, 2)

    @test abs(x - 0.1389) < 0.001

    x = DataAnalysis.calc_chi_squared(genotypes, 3)

    @test abs(x - 1.8) < 0.001
end


@testset "calc_all_trait_correlations" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 1], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)

    for tc in trait_correlations
        @test tc[3] == 1
    end

    genotypes = [[1 1 1; 1 1 1], [0 1 0; 0 1 0], [0 1 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)
    expected = [1, 0.9333, 0.93333, 0.9333, 1, 1, 0.9333, 1, 1]

    trait_correlations = [tc[3] for tc in trait_correlations]

    @test all(abs.(trait_correlations .- expected) .< 0.001)
end

@testset "calc_trait_correlation" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    @test abs(0.6708 - DataAnalysis.calc_trait_correlation(genotypes, [1, 3], [2])) < 0.0001
    @test abs(1 - DataAnalysis.calc_trait_correlation(genotypes, 1:3, 1:3)) < 0.0001
    @test abs(0.4472 - DataAnalysis.calc_trait_correlation(genotypes, [1], [3])) < 0.0001
    @test abs(0.8407 - DataAnalysis.calc_trait_correlation(genotypes, 1:2, 2:3)) < 0.0001

    genotypes = [[1 0 1; 0 1 1], [0 1 1; 1 0 1], [0 0 1; 1 1 1]]

    correlation = DataAnalysis.calc_trait_correlation(genotypes, [1, 2, 3], [3])

    @test correlation == 0
end

@testset "calc_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 1)

    @test linkage_diseq == 1

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 2)

    @test linkage_diseq ≈ 1 / 9

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 3)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 3, 1)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 2, 3)

    @test linkage_diseq ≈ 1 / 2
end

@testset "calc_average_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]
    loci = (overall=1:3, l1=1:2, l2=[1, 3])

    average_linkage_diseq = DataAnalysis.calc_average_linkage_diseq(genotypes, loci)

    expected = [2576 / 5184, 615 / 1296, 103 / 216, 205 / 432, 5 / 9, 125 / 288, 103 / 216, 125 / 288, 9 / 16]

    averages = vcat([a[3] for a in average_linkage_diseq]...)

    @test averages ≈ expected
end









@testset "calc_width" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [fill(0, 50); fill(1, 50)]
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test abs(width) <= 0.005

    hybrid_indices = vcat(fill([0, 1], 50)...)
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test width ≈ 1


    hybrid_indices = collect(0:0.01:0.99)
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test abs(width - 0.9) <= 0.01
end