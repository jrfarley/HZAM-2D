@testset "calculate_spatial_data" begin
    locations_x = [0.499f0, 0.495f0, 0.48f0, 0.52f0, 0.505f0, 0.501f0]
    locations_y = [0.0f0, 0.3f0, 0.6f0, 0.7f0, 0.8f0, 0.9f0]
    locations = Location.(locations_x, locations_y)
    genotypes = [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]
    loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3,], [1, 2], [3], 1:2, 1:2, 1:2, 1:2))

    output_data = DataAnalysis.SpatialData(
        locations,
        0.05,
        genotypes,
        loci,
        true
    )


    gene_flows = collect(values(output_data.gene_flows))

    [@test gene_flows[i] == 0.0 for i in eachindex(gene_flows)]
    readline()

    cline_widths = collect(values(output_data.cline_widths))

    @test length(cline_widths) == 7
    @test all(cline_widths[1] .== cline_widths)

    cline_positions = collect(values(output_data.cline_positions))

    @test length(cline_positions) == 7
    @test all(cline_positions[1] .== cline_positions)

    @test output_data.bimodality == 1
end

@testset "average_output_data" begin
    locations_x = [0.499f0, 0.495f0, 0.48f0, 0.52f0, 0.505f0, 0.501f0]
    locations_y = [0.0f0, 0.3f0, 0.6f0, 0.7f0, 0.8f0, 0.9f0]
    locations = Location.(locations_x, locations_y)
    genotypes = [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [1 1 1; 1 1 1],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]
    loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3,], [1, 2], [3], 1:2, 1:2, 1:2, 1:2))

    output_data = DataAnalysis.SpatialData(
        locations,
        0.05,
        genotypes,
        loci,
        true
    )

    outputs = fill(output_data, 20)
    average_output = DataAnalysis.SpatialData(outputs)

    @test average_output.bimodality ≈ output_data.bimodality
    for n in keys(loci)
        @test getfield(average_output.gene_flows, n) ≈ getfield(output_data.gene_flows, n)
        @test getfield(average_output.cline_widths, n) ≈ getfield(output_data.cline_widths, n)
        @test getfield(average_output.cline_positions, n) ≈ getfield(output_data.cline_positions, n)
    end
    @test average_output.variance ≈ 0
end

@testset "PopulationTrackingData" begin
    genotypes = [[0 0; 1 1], [0 1; 0 0], [0 0; 0 0], [1 1; 1 1]]
    locations_x = [0.1f0, 0.6f0, 0.2f0, 0.7f0]
    locations_y = fill(0.51f0, 4)
    locations = Location.(locations_x, locations_y)

    loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3,], [1, 2], [3], 1:2, 1:2, 1:2, 1:2))

    tracking_data = DataAnalysis.PopulationTrackingData(genotypes, locations, loci)

    @test tracking_data.population_size == 4

    @test tracking_data.hybridness ≈ 1.5 / 4

    @test tracking_data.overlap ≈ 0

    @test tracking_data.width > 0.95
end

@testset "average_per_phenotype" begin
    genotypes = [[0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 1; 0 0 0], [1 0 0; 0 0 1],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    num_offspring = [1, 2, 0, 1, 1, 1]

    different_hybrid_indices = (0:6) .* (1 / 6)
    different_hybrid_indices = round.(different_hybrid_indices, digits=4)
    different_hybrid_indices = Float32.(different_hybrid_indices)
    avg_offspring_per_phenotype = [1.5, 0, 1, NaN, NaN, NaN, 1]
    expected_fitnesses = Dict(different_hybrid_indices .=> avg_offspring_per_phenotype)

    fitnesses = DataAnalysis.average_data_per_phenotype(num_offspring, genotypes, 1:3)
    for k in collect(keys(fitnesses))
        @test isequal(fitnesses[k], expected_fitnesses[k])
    end
end

@testset "calc_traits_additive" begin
    genotypes = [[0 1 1; 0 1 0], [0 0 0; 0 0 0], [0 0 1; 0 0 0], [1 0 0; 0 0 1],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    t1 = DataAnalysis.calc_traits_additive(genotypes, [1])

    @test t1 ≈ [0, 0, 0, 0.5, 1, 1]

    t2 = DataAnalysis.calc_traits_additive(genotypes, [1, 3])

    @test t2 ≈ [0.25, 0, 0.25, 0.5, 1, 1]

    t3 = DataAnalysis.calc_traits_additive(genotypes, Int8[])

    @test all(isnan.(t3))
end

@testset "average_gene_data" begin
    loci = NamedTuple{(:l1, :l2, :l3, :l4),NTuple{4,Any}}((1:2, 3:3, 1:3, [1, 3]))
    function generate_gene_data(loci_range, counter)
        length(loci_range) + counter
    end

    gene_data = NamedTuple[]

    for i in 1:5
        push!(gene_data, map(x -> generate_gene_data(x, i), loci))
    end

    avg_gene_data = DataAnalysis.average_gene_data(gene_data)

    @test avg_gene_data == (l1=5, l2=4, l3=6, l4=5)
end

@testset "sort_locations" begin
    A = [0, 0, 0.5, 1, 0.76, 0.23]
    bin_size = 0.5

    sorted_indices = DataAnalysis.sort_locations(A, bin_size)

    @test sorted_indices == [[1, 2, 6], [3, 5]]

    bin_size = 0.1

    sorted_indices = DataAnalysis.sort_locations(A, bin_size)

    @test sorted_indices == [[1, 2], [], [6], [], [], [3], [], [5], [], []]
end

@testset "sort_y" begin
    @test DataAnalysis.sort_y([0.01, 0.5, 0.24, 0.9]) ==
          [[1], [], [3], [], [], [2], [], [], [], [4]]
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

@testset "calc_bimodality_in_range" begin
    locations_x = collect(0:0.01:0.99)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    bimodality = DataAnalysis.calc_bimodality_in_range(
        sigmoid_curve,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality == 1

    hybrid_indices = [fill(0, 49); [0.5]; fill(1, 50)]

    bimodality = DataAnalysis.calc_bimodality_in_range(
        sigmoid_curve,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality ≈ 4 / 5
end

@testset "calc_bimodality_overall" begin
    locations_x = collect(0:0.01:0.99)
    hybrid_indices = [fill(0, 50); fill(1, 50)]
    sorted_indices = [collect(((i-1)*100+1):(((i-1)*100+1)+99)) for i in 1:10]
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)
    sigmoid_curves = fill(sigmoid_curve, 10)

    locations_x = vcat(fill(locations_x, 10)...)
    hybrid_indices = vcat(fill(hybrid_indices, 10)...)

    bimodality = DataAnalysis.calc_bimodality_overall(
        sigmoid_curves,
        sorted_indices,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality == 1

    h0 = [fill(0, 50); fill(1, 50)]
    h1 = [fill(0, 49); [0.95]; fill(1, 50)]

    hybrid_indices = vcat(h0, h1, h0, h1, h0, h1, h0, h1, h0, h1)

    bimodality = DataAnalysis.calc_bimodality_overall(
        sigmoid_curves,
        sorted_indices,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality ≈ 0.9
end

@testset "calc_overlap_in_range" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    overlap = DataAnalysis.calc_overlap_in_range(
        locations_x,
        hybrid_indices
    )

    @test overlap == 0

    hybrid_indices = fill(0.99, 100)

    overlap = DataAnalysis.calc_overlap_in_range(
        locations_x,
        hybrid_indices
    )

    @test overlap == 0

    hybrid_indices = vcat(fill([0.0f0, 1.0f0], 50)...)

    overlap = DataAnalysis.calc_overlap_in_range(
        locations_x,
        hybrid_indices
    )

    @test overlap ≈ 0.1
end

@testset "calc_overlap_overall" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [fill(0, 50); fill(1, 50)]
    sorted_indices = [collect(((i-1)*100+1):(((i-1)*100+1)+99)) for i in 1:10]

    locations_x = vcat(fill(locations_x, 10)...)
    hybrid_indices = vcat(fill(hybrid_indices, 10)...)

    overlap = DataAnalysis.calc_overlap_overall(
        locations_x,
        hybrid_indices,
        sorted_indices
    )

    @test overlap ≈ 0

    h0 = [fill(0, 50); fill(1, 50)]
    h1 = vcat(fill([0.0f0, 1.0f0], 50)...)

    hybrid_indices = vcat(h0, h1, h0, h1, h0, h1, h0, h1, h0, h1)

    overlap = DataAnalysis.calc_overlap_overall(
        locations_x,
        hybrid_indices,
        sorted_indices
    )

    @test overlap ≈ 0.5
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

@testset "average_width" begin
    locations_x = collect(0.001:0.01:0.999)
    h0 = [fill(0, 50); fill(1, 50)]
    h1 = vcat(fill([0, 1], 50)...)
    s0 = DataAnalysis.calc_sigmoid_curve(locations_x, h0)
    s1 = DataAnalysis.calc_sigmoid_curve(locations_x, h1)

    avg = DataAnalysis.average_width(vcat(fill([s0, s1], 5)...))

    @test abs(0.5 - avg) < 0.005
end

@testset "sigmoid" begin
    xs = collect(0:0.01:1)

    s = DataAnalysis.sigmoid(xs, [0, 1])

    @test s[1] ≈ 0.5
    @test s[101] ≈ 1 / (1 + exp(-1))
end

@testset "calc_sigmoid_curve" begin
    locations_x = collect(0:0.01:0.99)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    s = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)
    midpoint = DataAnalysis.spaced_locations[argmin(abs.(s .- 0.5))]
    @test length(s) == 1001
    @test midpoint ≈ 0.495
    @test s[1] == 0
    @test s[1000] == 1
end

@testset "calc_sigmoid_curves" begin
    locations_x = collect(0.0f0:0.01f0:0.99f0)
    locations = Vector(undef, 10)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    for i in 1:10
        locations[i] = Location.(locations_x, Ref(Float32((i - 1) * 0.1) + 0.001f0))
    end

    locations = vcat(locations...)

    hybrid_indices = vcat(fill(hybrid_indices, 10)...)

    @test length(locations) == length(hybrid_indices)

    sigmoid_curves = DataAnalysis.calc_sigmoid_curves(locations, hybrid_indices)
    s = sigmoid_curves[1]

    midpoint = DataAnalysis.spaced_locations[argmin(abs.(s .- 0.5))]
    @test all(Ref(s) .== sigmoid_curves)
    @test length(s) == 1001
    @test midpoint ≈ 0.495
    @test s[1] == 0
    @test s[1000] == 1
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