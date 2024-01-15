@testset "PopulationTrackingData" begin
    genotypes = [[0 0; 1 1], [0 1; 0 0], [0 0; 0 0], [1 1; 1 1]]
    locations_x = [0.1f0, 0.6f0, 0.2f0, 0.7f0]
    locations_y = fill(0.51f0, 4)

    loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3,], [1, 2], [3], 1:2, 1:2, 1:2, 1:2))

    tracking_data = DataAnalysis.PopulationTrackingData(genotypes, locations_x, locations_y, (1:2), 0)

    @test tracking_data.population_size == 4

    @test tracking_data.hybridity ≈ 1.5 / 4

    @test tracking_data.overlap ≈ 0

    @test tracking_data.width > 0.95
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
          [[1], [], [], [], [3], [], [], [], [], [], [2], [], [], [], [], [], [], [], [4], []]
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
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    locations_x = vcat(fill(locations_x, 20)...)
    locations_y = vcat([fill(Float32((i - 1) * 0.05) + 0.001f0, 100) for i in 1:20]...)

    hybrid_indices = vcat(fill(hybrid_indices, 20)...)

    @test length(locations_x) == length(hybrid_indices)

    sigmoid_curves = DataAnalysis.calc_sigmoid_curves(locations_x, locations_y, hybrid_indices)
    s = sigmoid_curves[1]

    midpoint = DataAnalysis.spaced_locations[argmin(abs.(s .- 0.5))]
    @test maximum(sigmoid_curves[1].-sigmoid_curves[10]) <0.01
    @test length(s) == 1001
    @test midpoint ≈ 0.495
    @test s[1] == 0
    @test s[1000] == 1
end
